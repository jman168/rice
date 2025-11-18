mod matrix;
use matrix::{ABMatrix, EquationIndex, VariableIndex};

use crate::{
    be_solver::matrix::XMatrix,
    components::{Capacitor, CurrentSource, Netlist, Resistor, VoltageSource},
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
    pub fn solve(&mut self, dt: f64) {
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
            .for_each(|(i, c)| Self::stamp_resistor(c, i, &mut problem, dt));

        self.netlist
            .get_capacitors()
            .iter()
            .enumerate()
            .for_each(|(i, c)| Self::stamp_capacitor(c, i, &mut problem, dt));

        self.netlist
            .get_voltages_sources()
            .iter()
            .enumerate()
            .for_each(|(i, c)| Self::stamp_voltage_source(c, i, &mut problem, dt));

        self.netlist
            .get_current_sources()
            .iter()
            .enumerate()
            .for_each(|(i, c)| Self::stamp_current_source(c, i, &mut problem, dt));

        let (a, b) = problem.into_ab();
        let solution = XMatrix::new(a.try_inverse().unwrap() * b, num_nodes, num_voltage_sources);

        self.netlist
            .get_resistors_mut()
            .iter_mut()
            .enumerate()
            .for_each(|(i, c)| Self::update_resistor(c, i, &solution, dt));

        self.netlist
            .get_capacitors_mut()
            .iter_mut()
            .enumerate()
            .for_each(|(i, c)| Self::update_capacitor(c, i, &solution, dt));

        self.netlist
            .get_voltages_sources_mut()
            .iter_mut()
            .enumerate()
            .for_each(|(i, c)| Self::update_voltage_source(c, i, &solution, dt));

        self.netlist
            .get_current_sources_mut()
            .iter_mut()
            .enumerate()
            .for_each(|(i, c)| Self::update_current_source(c, i, &solution, dt));
    }

    fn stamp_resistor(
        resistor: &Resistor,
        _resistor_index: usize,
        problem: &mut ABMatrix,
        _dt: f64,
    ) {
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

    fn stamp_capacitor(
        capacitor: &Capacitor,
        _capacitor_index: usize,
        problem: &mut ABMatrix,
        dt: f64,
    ) {
        let positive_equation_index = EquationIndex::NodalEquation(capacitor.get_positive_node());
        let negative_equation_index = EquationIndex::NodalEquation(capacitor.get_negative_node());

        let positive_voltage_index = VariableIndex::NodalVoltage(capacitor.get_positive_node());
        let negative_voltage_index = VariableIndex::NodalVoltage(capacitor.get_negative_node());

        let c = capacitor.get_capacitance();

        // The differential equation describing a capacitor is i = C*dv/dt.
        // Discretizing we get i = C*(v_new - v_old)/dt.
        // Further expanding this we get i = C*(v_positive - v_negative)/dt - C*v_old/dt.

        // Current flowing out of the positive node is C*v_positive/dt - C*v_negative/dt - C*v_old/dt.
        problem.coefficient_add(positive_equation_index, positive_voltage_index, c / dt);
        problem.coefficient_add(positive_equation_index, negative_voltage_index, -c / dt);
        problem.result_add(positive_equation_index, c * capacitor.get_voltage() / dt);

        // Current flowing out of the negative node is -C*v_positive/dt + C*v_negative/dt + C*v_old/dt.
        problem.coefficient_add(negative_equation_index, positive_voltage_index, -c / dt);
        problem.coefficient_add(negative_equation_index, negative_voltage_index, c / dt);
        problem.result_add(negative_equation_index, -c * capacitor.get_voltage() / dt);
    }

    fn stamp_voltage_source(
        voltage_source: &VoltageSource,
        voltage_source_index: usize,
        problem: &mut ABMatrix,
        _dt: f64,
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

    fn stamp_current_source(
        current_source: &CurrentSource,
        _current_source_index: usize,
        problem: &mut ABMatrix,
        _dt: f64,
    ) {
        let positive_equation_index =
            EquationIndex::NodalEquation(current_source.get_positive_node());
        let negative_equation_index =
            EquationIndex::NodalEquation(current_source.get_negative_node());

        // NOTE: the signs are flipped here because they take the form of constants, not
        // coefficients.

        // Current flowing out of positive node is -i_source
        problem.result_add(positive_equation_index, current_source.get_current());
        // Current flowing out of negative node is i_source
        problem.result_add(negative_equation_index, -current_source.get_current());
    }

    fn update_resistor(
        resistor: &mut Resistor,
        _resistor_index: usize,
        solution: &XMatrix,
        _dt: f64,
    ) {
        let positive_voltage_index = VariableIndex::NodalVoltage(resistor.get_positive_node());
        let negative_voltage_index = VariableIndex::NodalVoltage(resistor.get_negative_node());

        resistor.set_voltage(
            solution.get_variable(positive_voltage_index).unwrap()
                - solution.get_variable(negative_voltage_index).unwrap(),
        );
    }

    fn update_capacitor(
        capacitor: &mut Capacitor,
        _capacitor_index: usize,
        solution: &XMatrix,
        dt: f64,
    ) {
        let positive_voltage_index = VariableIndex::NodalVoltage(capacitor.get_positive_node());
        let negative_voltage_index = VariableIndex::NodalVoltage(capacitor.get_negative_node());

        let new_voltage = solution.get_variable(positive_voltage_index).unwrap()
            - solution.get_variable(negative_voltage_index).unwrap();

        // Discretized equation is i = C*(v_new - v_old)/dt (see capacitor stamping function).

        capacitor.set_current(
            capacitor.get_capacitance() * (new_voltage - capacitor.get_voltage()) / dt,
        );

        capacitor.set_voltage(new_voltage);
    }

    fn update_voltage_source(
        voltage_source: &mut VoltageSource,
        voltage_source_index: usize,
        solution: &XMatrix,
        _dt: f64,
    ) {
        let source_current_index = VariableIndex::VoltageSourceCurrent(voltage_source_index);

        voltage_source.set_current(solution.get_variable(source_current_index).unwrap());
    }

    fn update_current_source(
        current_source: &mut CurrentSource,
        _current_source_index: usize,
        solution: &XMatrix,
        _dt: f64,
    ) {
        let positive_voltage_index =
            VariableIndex::NodalVoltage(current_source.get_positive_node());
        let negative_voltage_index =
            VariableIndex::NodalVoltage(current_source.get_negative_node());

        current_source.set_voltage(
            solution.get_variable(positive_voltage_index).unwrap()
                - solution.get_variable(negative_voltage_index).unwrap(),
        );
    }
}

#[cfg(test)]
mod test {
    use crate::{
        BESolver,
        components::{Capacitor, CurrentSource, Netlist, Resistor, VoltageSource},
    };

    use approx::assert_relative_eq;

    #[test]
    fn test_voltage_source_resistor() {
        let mut netlist = Netlist::new();
        netlist
            .add_voltage_source(VoltageSource::new(1, 0, 10.0))
            .add_resistor(Resistor::new(1, 0, 2.0));

        let mut solver = BESolver::new(&mut netlist);
        solver.solve(0.001);

        println!("{:?}", netlist);

        assert_relative_eq!(
            netlist.get_voltages_sources()[0].get_voltage(),
            10.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_voltages_sources()[0].get_current(),
            5.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_resistors()[0].get_voltage(),
            10.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_resistors()[0].get_current(),
            5.0,
            max_relative = 0.001
        );
    }

    #[test]
    fn test_two_voltage_source_resistor() {
        let mut netlist = Netlist::new();
        netlist
            .add_voltage_source(VoltageSource::new(1, 0, 10.0))
            .add_resistor(Resistor::new(1, 2, 2.0))
            .add_voltage_source(VoltageSource::new(2, 0, 5.0));

        let mut solver = BESolver::new(&mut netlist);
        solver.solve(0.001);

        println!("{:?}", netlist);

        assert_relative_eq!(
            netlist.get_voltages_sources()[0].get_voltage(),
            10.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_voltages_sources()[0].get_current(),
            2.5,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_voltages_sources()[1].get_voltage(),
            5.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_voltages_sources()[1].get_current(),
            -2.5,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_resistors()[0].get_voltage(),
            5.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_resistors()[0].get_current(),
            2.5,
            max_relative = 0.001
        );
    }

    #[test]
    fn test_voltage_source_voltage_divider() {
        let mut netlist = Netlist::new();
        netlist
            .add_voltage_source(VoltageSource::new(1, 0, 5.0))
            .add_resistor(Resistor::new(1, 2, 4.0))
            .add_resistor(Resistor::new(2, 0, 1.0));

        let mut solver = BESolver::new(&mut netlist);
        solver.solve(0.001);

        println!("{:?}", netlist);

        assert_relative_eq!(
            netlist.get_voltages_sources()[0].get_voltage(),
            5.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_voltages_sources()[0].get_current(),
            1.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_resistors()[0].get_voltage(),
            4.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_resistors()[0].get_current(),
            1.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_resistors()[1].get_voltage(),
            1.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_resistors()[1].get_current(),
            1.0,
            max_relative = 0.001
        );
    }

    #[test]
    fn test_current_source_resistor() {
        let mut netlist = Netlist::new();
        netlist
            .add_current_source(CurrentSource::new(1, 0, 5.0))
            .add_resistor(Resistor::new(1, 0, 2.0));

        let mut solver = BESolver::new(&mut netlist);
        solver.solve(0.001);

        println!("{:?}", netlist);

        assert_relative_eq!(
            netlist.get_current_sources()[0].get_voltage(),
            10.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_current_sources()[0].get_current(),
            5.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_resistors()[0].get_voltage(),
            10.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_resistors()[0].get_current(),
            5.0,
            max_relative = 0.001
        );
    }

    #[test]
    fn test_current_source_capacitor() {
        let mut netlist = Netlist::new();
        netlist
            .add_current_source(CurrentSource::new(1, 0, 1.0))
            .add_capacitor(Capacitor::new(1, 0, 0.5, 0.0));

        let mut solver = BESolver::new(&mut netlist);
        solver.solve(0.25);

        println!("{:?}", netlist);

        assert_relative_eq!(
            netlist.get_current_sources()[0].get_voltage(),
            0.50,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_current_sources()[0].get_current(),
            1.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_capacitors()[0].get_voltage(),
            0.50,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_capacitors()[0].get_current(),
            1.0,
            max_relative = 0.001
        );

        let mut solver = BESolver::new(&mut netlist);
        solver.solve(0.75);

        println!("{:?}", netlist);

        assert_relative_eq!(
            netlist.get_current_sources()[0].get_voltage(),
            2.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_current_sources()[0].get_current(),
            1.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_capacitors()[0].get_voltage(),
            2.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_capacitors()[0].get_current(),
            1.0,
            max_relative = 0.001
        );
    }

    #[test]
    fn test_rc_step() {
        let mut netlist = Netlist::new();
        netlist
            .add_voltage_source(VoltageSource::new(1, 0, 1.0))
            .add_resistor(Resistor::new(1, 2, 1000.0))
            .add_capacitor(Capacitor::new(2, 0, 0.001, 0.0));

        let mut solver = BESolver::new(&mut netlist);
        for _ in 0..1000 {
            solver.solve(0.001);
        }

        println!("{:?}", netlist);

        assert_relative_eq!(
            netlist.get_voltages_sources()[0].get_voltage(),
            1.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_voltages_sources()[0].get_current(),
            0.000367879441171,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_resistors()[0].get_voltage(),
            0.367879441171,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_resistors()[0].get_current(),
            0.000367879441171,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_capacitors()[0].get_voltage(),
            0.632120558829,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_capacitors()[0].get_current(),
            0.000367879441171,
            max_relative = 0.001
        );
    }
}
