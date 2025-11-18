mod matrix;
use matrix::{ABMatrix, EquationIndex, VariableIndex};

use crate::{
    be_solver::matrix::XMatrix,
    components::{Netlist, Resistor, VoltageSource},
};

/// A Backward Euler method solver for solving transient circuits.
pub struct BESolver<'n> {
    netlist: &'n mut Netlist,
}

impl<'n> BESolver<'n> {
    /// Creates a new BESolver with a given number of nodes.
    pub fn new(netlist: &'n mut Netlist) -> Self {
        Self { netlist }
    }

    /// Solves the system for the next timestep dt.
    pub fn solve(&mut self, _dt: f64) {
        // Compute the dimensionality of the matrix we are to solve.
        //
        // This is the number of nodes plus the number of voltages sources.
        //
        // This is because for each node we have a variable (node voltages) and an equation (KCL at
        // that node).
        //
        // For each voltages sources we have a variable (current through the voltage source) and an
        // equation (setting the voltage potential between the two nodes).
        let num_nodes = self.netlist.get_num_nodes();
        let num_voltage_sources = self.netlist.get_voltages_sources().len();

        let mut problem = ABMatrix::new(num_nodes, num_voltage_sources);

        self.netlist
            .get_resistors()
            .iter()
            .enumerate()
            .for_each(|(i, r)| Self::stamp_resistor(r, i, &mut problem));

        self.netlist
            .get_voltages_sources()
            .iter()
            .enumerate()
            .for_each(|(i, v)| Self::stamp_voltage_source(v, i, &mut problem));

        let (a, b) = problem.into_ab();
        let solution = XMatrix::new(a.try_inverse().unwrap() * b, num_nodes, num_voltage_sources);

        self.netlist
            .get_resistors_mut()
            .iter_mut()
            .enumerate()
            .for_each(|(i, r)| Self::update_resistor(r, i, &solution));

        self.netlist
            .get_voltages_sources_mut()
            .iter_mut()
            .enumerate()
            .for_each(|(i, v)| Self::update_voltage_source(v, i, &solution));
    }

    fn stamp_resistor(resistor: &Resistor, _resistor_index: usize, problem: &mut ABMatrix) {
        let positive_equation_index = EquationIndex::NodalEquation(resistor.get_positive_node());
        let negative_equation_index = EquationIndex::NodalEquation(resistor.get_negative_node());

        let positive_voltage_index = VariableIndex::NodalVoltage(resistor.get_positive_node());
        let negative_voltage_index = VariableIndex::NodalVoltage(resistor.get_negative_node());

        // Compute resistance conductance
        let g = 1.0 / resistor.get_resistance();

        // Current flowing out of positive node is (v_positive - v_negative) / R
        problem.coefficient_add(positive_equation_index, positive_voltage_index, g);
        problem.coefficient_add(positive_equation_index, negative_voltage_index, -g);

        // Current flowing out of positive node is -(v_positive - v_negative) / R
        problem.coefficient_add(negative_equation_index, positive_voltage_index, -g);
        problem.coefficient_add(negative_equation_index, negative_voltage_index, g);
    }

    fn stamp_voltage_source(
        voltage_source: &VoltageSource,
        voltage_source_index: usize,
        problem: &mut ABMatrix,
    ) {
        let positive_equation_index =
            EquationIndex::NodalEquation(voltage_source.get_positive_node());
        let negative_equation_index =
            EquationIndex::NodalEquation(voltage_source.get_negative_node());
        let source_equation_index = EquationIndex::VoltageSourceEquation(voltage_source_index);

        let positive_voltage_index =
            VariableIndex::NodalVoltage(voltage_source.get_positive_node());
        let negative_voltage_index =
            VariableIndex::NodalVoltage(voltage_source.get_negative_node());
        let source_current_index = VariableIndex::VoltageSourceCurrent(voltage_source_index);

        // Current flowing out of positive node is -i_source
        problem.coefficient_add(positive_equation_index, source_current_index, -1.0);
        // Current flowing out of negative node is i_source
        problem.coefficient_add(negative_equation_index, source_current_index, 1.0);

        // Source equation is v_positive - v_negative = v_source
        problem.coefficient_add(source_equation_index, positive_voltage_index, 1.0);
        problem.coefficient_add(source_equation_index, negative_voltage_index, -1.0);
        problem.result_add(source_equation_index, voltage_source.get_voltage());
    }

    fn update_resistor(resistor: &mut Resistor, _resistor_index: usize, solution: &XMatrix) {
        let positive_voltage_index = VariableIndex::NodalVoltage(resistor.get_positive_node());
        let negative_voltage_index = VariableIndex::NodalVoltage(resistor.get_negative_node());

        resistor.set_voltage(
            solution.get_variable(positive_voltage_index).unwrap()
                - solution.get_variable(negative_voltage_index).unwrap(),
        );
    }

    fn update_voltage_source(
        voltage_source: &mut VoltageSource,
        voltage_source_index: usize,
        solution: &XMatrix,
    ) {
        let source_current_index = VariableIndex::VoltageSourceCurrent(voltage_source_index);

        voltage_source.set_current(solution.get_variable(source_current_index).unwrap());
    }
}

#[cfg(test)]
mod test {
    use crate::{
        BESolver,
        components::{Netlist, Resistor, VoltageSource},
    };

    #[test]
    fn test_single_resistor() {
        let mut netlist = Netlist::new();
        netlist
            .add_voltage_source(VoltageSource::new(1, 0, 10.0))
            .add_resistor(Resistor::new(1, 0, 2.0));

        let mut solver = BESolver::new(&mut netlist);
        solver.solve(0.001);

        println!("{:?}", netlist);

        assert_eq!(netlist.get_voltages_sources()[0].get_current(), 5.0);
        assert_eq!(netlist.get_voltages_sources()[0].get_power(), 50.0);
        assert_eq!(netlist.get_resistors()[0].get_voltage(), 10.0);
        assert_eq!(netlist.get_resistors()[0].get_current(), 5.0);
        assert_eq!(netlist.get_resistors()[0].get_power(), 50.0);
    }

    #[test]
    fn test_single_resistor_source_flip() {
        let mut netlist = Netlist::new();
        netlist
            .add_voltage_source(VoltageSource::new(0, 1, 10.0))
            .add_resistor(Resistor::new(1, 0, 2.0));

        let mut solver = BESolver::new(&mut netlist);
        solver.solve(0.001);

        println!("{:?}", netlist);

        assert_eq!(netlist.get_voltages_sources()[0].get_current(), 5.0);
        assert_eq!(netlist.get_voltages_sources()[0].get_power(), 50.0);
        assert_eq!(netlist.get_resistors()[0].get_voltage(), -10.0);
        assert_eq!(netlist.get_resistors()[0].get_current(), -5.0);
        assert_eq!(netlist.get_resistors()[0].get_power(), 50.0);
    }

    #[test]
    fn test_voltage_divider() {
        let mut netlist = Netlist::new();
        netlist
            .add_voltage_source(VoltageSource::new(1, 0, 5.0))
            .add_resistor(Resistor::new(1, 2, 4.0))
            .add_resistor(Resistor::new(2, 0, 1.0));

        let mut solver = BESolver::new(&mut netlist);
        solver.solve(0.001);

        println!("{:?}", netlist);

        assert_eq!(netlist.get_voltages_sources()[0].get_current(), 1.0);
        assert_eq!(netlist.get_voltages_sources()[0].get_power(), 5.0);
        assert_eq!(netlist.get_resistors()[0].get_voltage(), 4.0);
        assert_eq!(netlist.get_resistors()[0].get_current(), 1.0);
        assert_eq!(netlist.get_resistors()[0].get_power(), 4.0);
        assert_eq!(netlist.get_resistors()[1].get_voltage(), 1.0);
        assert_eq!(netlist.get_resistors()[1].get_current(), 1.0);
        assert_eq!(netlist.get_resistors()[1].get_power(), 1.0);
    }

    // #[test]
    // fn test_complex() {
    //     let mut netlist = Netlist::new();
    //     netlist
    //         .add_voltage_source(VoltageSource::new(1, 0, 10.0))
    //         .add_resistor(Resistor::new(1, 2, 5.0))
    //         .add_resistor(Resistor::new(1, 3, 3.0))
    //         .add_resistor(Resistor::new(2, 4, 3.0))
    //         .add_voltage_source(VoltageSource::new(3, 4, 2.0))
    //         .add_resistor(Resistor::new(4, 0, 5.0))
    //         .add_resistor(Resistor::new(4, 0, 4.0));
    //
    //     let mut solver = BESolver::new(&mut netlist);
    //     solver.solve(0.001);
    //
    //     println!("{:?}", netlist);
    // }
}
