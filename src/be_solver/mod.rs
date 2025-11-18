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
        let p_eq = EquationIndex::NodalEquation(resistor.get_positive_node());
        let n_eq = EquationIndex::NodalEquation(resistor.get_negative_node());

        let p_v = VariableIndex::NodalVoltage(resistor.get_positive_node());
        let n_v = VariableIndex::NodalVoltage(resistor.get_negative_node());

        // Compute resistance conductance
        let g = 1.0 / resistor.get_resistance();

        problem.get_coefficient_mut(p_eq, p_v).map(|a| *a += g);
        problem.get_coefficient_mut(p_eq, n_v).map(|a| *a -= g);

        problem.get_coefficient_mut(n_eq, p_v).map(|a| *a += g);
        problem.get_coefficient_mut(n_eq, n_v).map(|a| *a -= g);
    }

    fn stamp_voltage_source(
        voltage_source: &VoltageSource,
        voltage_source_index: usize,
        problem: &mut ABMatrix,
    ) {
        let p_eq = EquationIndex::NodalEquation(voltage_source.get_positive_node());
        let n_eq = EquationIndex::NodalEquation(voltage_source.get_negative_node());
        let v_eq = EquationIndex::VoltageSourceEquation(voltage_source_index);

        let p_v = VariableIndex::NodalVoltage(voltage_source.get_positive_node());
        let n_v = VariableIndex::NodalVoltage(voltage_source.get_negative_node());
        let v_i = VariableIndex::VoltageSourceCurrent(voltage_source_index);

        problem.get_coefficient_mut(p_eq, v_i).map(|a| *a -= 1.0);
        problem.get_coefficient_mut(n_eq, v_i).map(|a| *a += 1.0);

        problem.get_coefficient_mut(v_eq, p_v).map(|a| *a += 1.0);
        problem.get_coefficient_mut(v_eq, n_v).map(|a| *a -= 1.0);
        problem
            .get_result_mut(v_eq)
            .map(|a| *a += voltage_source.get_voltage());
    }

    fn update_resistor(resistor: &mut Resistor, _resistor_index: usize, solution: &XMatrix) {
        let p_v = VariableIndex::NodalVoltage(resistor.get_positive_node());
        let n_v = VariableIndex::NodalVoltage(resistor.get_negative_node());

        resistor
            .set_voltage(solution.get_variable(p_v).unwrap() - solution.get_variable(n_v).unwrap());
    }

    fn update_voltage_source(
        voltage_source: &mut VoltageSource,
        voltage_source_index: usize,
        solution: &XMatrix,
    ) {
        let v_i = VariableIndex::VoltageSourceCurrent(voltage_source_index);

        voltage_source.set_current(solution.get_variable(v_i).unwrap());
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
        assert_eq!(netlist.get_resistors()[0].get_current(), -5.0);
        assert_eq!(netlist.get_resistors()[0].get_power(), 50.0);
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
