mod matrix_view;
mod stampable;

use nalgebra::DMatrix;

use matrix_view::{ABMatrixView, XMatrixView};
use stampable::Stampable;

use crate::components::Netlist;

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
        let num_variables: usize = self
            .netlist
            .get_components()
            .iter()
            .map(|c| c.num_variables())
            .sum();

        let mut a = DMatrix::zeros(num_nodes + num_variables, num_nodes + num_variables);

        let mut b = DMatrix::zeros(num_nodes + num_variables, 1);

        self.netlist
            .get_components()
            .iter()
            .fold(num_nodes, |variables_start, c| {
                let mut view = ABMatrixView::new(
                    &mut a,
                    &mut b,
                    num_nodes,
                    c.num_variables(),
                    variables_start,
                );
                c.stamp(&mut view, dt);
                variables_start + c.num_variables()
            });

        let x = a.try_inverse().unwrap() * b;

        self.netlist
            .get_components_mut()
            .iter_mut()
            .fold(num_nodes, |variables_start, c| {
                let view = XMatrixView::new(&x, num_nodes, c.num_variables(), variables_start);
                c.update(&view, dt);
                variables_start + c.num_variables()
            });
    }
}

#[cfg(test)]
mod test {
    use crate::{
        BESolver,
        components::{Capacitor, CurrentSource, Inductor, Netlist, Resistor, VoltageSource},
    };

    use approx::assert_relative_eq;

    #[test]
    fn test_voltage_source_resistor() {
        let mut netlist = Netlist::new();
        netlist
            .add_component(VoltageSource::new(1, 0, 10.0))
            .add_component(Resistor::new(1, 0, 2.0));

        let mut solver = BESolver::new(&mut netlist);
        solver.solve(0.001);

        println!("{:?}", netlist);

        assert_relative_eq!(
            netlist.get_components()[0].get_voltage(),
            10.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[0].get_current(),
            5.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[1].get_voltage(),
            10.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[1].get_current(),
            5.0,
            max_relative = 0.001
        );
    }

    #[test]
    fn test_two_voltage_source_resistor() {
        let mut netlist = Netlist::new();
        netlist
            .add_component(VoltageSource::new(1, 0, 10.0))
            .add_component(Resistor::new(1, 2, 2.0))
            .add_component(VoltageSource::new(2, 0, 5.0));

        let mut solver = BESolver::new(&mut netlist);
        solver.solve(0.001);

        println!("{:?}", netlist);

        assert_relative_eq!(
            netlist.get_components()[0].get_voltage(),
            10.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[0].get_current(),
            2.5,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[1].get_voltage(),
            5.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[1].get_current(),
            2.5,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[2].get_voltage(),
            5.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[2].get_current(),
            -2.5,
            max_relative = 0.001
        );
    }

    #[test]
    fn test_voltage_source_voltage_divider() {
        let mut netlist = Netlist::new();
        netlist
            .add_component(VoltageSource::new(1, 0, 5.0))
            .add_component(Resistor::new(1, 2, 4.0))
            .add_component(Resistor::new(2, 0, 1.0));

        let mut solver = BESolver::new(&mut netlist);
        solver.solve(0.001);

        println!("{:?}", netlist);

        assert_relative_eq!(
            netlist.get_components()[0].get_voltage(),
            5.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[0].get_current(),
            1.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[1].get_voltage(),
            4.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[1].get_current(),
            1.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[2].get_voltage(),
            1.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[2].get_current(),
            1.0,
            max_relative = 0.001
        );
    }

    #[test]
    fn test_current_source_resistor() {
        let mut netlist = Netlist::new();
        netlist
            .add_component(CurrentSource::new(1, 0, 5.0))
            .add_component(Resistor::new(1, 0, 2.0));

        let mut solver = BESolver::new(&mut netlist);
        solver.solve(0.001);

        println!("{:?}", netlist);

        assert_relative_eq!(
            netlist.get_components()[0].get_voltage(),
            10.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[0].get_current(),
            5.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[1].get_voltage(),
            10.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[1].get_current(),
            5.0,
            max_relative = 0.001
        );
    }

    #[test]
    fn test_current_source_capacitor() {
        let mut netlist = Netlist::new();
        netlist
            .add_component(CurrentSource::new(1, 0, 1.0))
            .add_component(Capacitor::new(1, 0, 0.5, 0.0));

        let mut solver = BESolver::new(&mut netlist);
        solver.solve(0.25);

        println!("{:?}", netlist);

        assert_relative_eq!(
            netlist.get_components()[0].get_voltage(),
            0.50,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[0].get_current(),
            1.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[1].get_voltage(),
            0.50,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[1].get_current(),
            1.0,
            max_relative = 0.001
        );

        let mut solver = BESolver::new(&mut netlist);
        solver.solve(0.75);

        println!("{:?}", netlist);

        assert_relative_eq!(
            netlist.get_components()[0].get_voltage(),
            2.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[0].get_current(),
            1.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[1].get_voltage(),
            2.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[1].get_current(),
            1.0,
            max_relative = 0.001
        );
    }

    #[test]
    fn test_rc_step() {
        let mut netlist = Netlist::new();
        netlist
            .add_component(VoltageSource::new(1, 0, 1.0))
            .add_component(Resistor::new(1, 2, 1000.0))
            .add_component(Capacitor::new(2, 0, 0.001, 0.0));

        let mut solver = BESolver::new(&mut netlist);
        for _ in 0..1000 {
            solver.solve(0.001);
        }

        println!("{:?}", netlist);

        assert_relative_eq!(
            netlist.get_components()[0].get_voltage(),
            1.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[0].get_current(),
            0.000367879441171,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[1].get_voltage(),
            0.367879441171,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[1].get_current(),
            0.000367879441171,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[2].get_voltage(),
            0.632120558829,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[2].get_current(),
            0.000367879441171,
            max_relative = 0.001
        );
    }

    #[test]
    fn test_voltage_source_inductor() {
        let mut netlist = Netlist::new();
        netlist
            .add_component(VoltageSource::new(1, 0, 1.0))
            .add_component(Inductor::new(1, 0, 0.5, 0.0));

        let mut solver = BESolver::new(&mut netlist);
        solver.solve(0.25);

        println!("{:?}", netlist);

        assert_relative_eq!(
            netlist.get_components()[0].get_voltage(),
            1.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[0].get_current(),
            0.5,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[1].get_voltage(),
            1.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[1].get_current(),
            0.5,
            max_relative = 0.001
        );

        let mut solver = BESolver::new(&mut netlist);
        solver.solve(0.75);

        println!("{:?}", netlist);

        assert_relative_eq!(
            netlist.get_components()[0].get_voltage(),
            1.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[0].get_current(),
            2.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[1].get_voltage(),
            1.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[1].get_current(),
            2.0,
            max_relative = 0.001
        );
    }

    #[test]
    fn test_rl_step() {
        let mut netlist = Netlist::new();
        netlist
            .add_component(VoltageSource::new(1, 0, 1.0))
            .add_component(Resistor::new(1, 2, 0.001))
            .add_component(Inductor::new(2, 0, 0.01, 0.0));

        let mut solver = BESolver::new(&mut netlist);
        for _ in 0..1000 {
            solver.solve(0.001);
        }

        println!("{:?}", netlist);

        assert_relative_eq!(
            netlist.get_components()[0].get_voltage(),
            1.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[0].get_current(),
            95.162581964,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[1].get_voltage(),
            0.095162581964,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[1].get_current(),
            95.162581964,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[2].get_voltage(),
            0.904837418036,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_components()[2].get_current(),
            95.162581964,
            max_relative = 0.001
        );
    }
}
