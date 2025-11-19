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
        let num_voltage_sources = self.netlist.get_voltages_sources().len();

        let mut a = DMatrix::zeros(
            num_nodes + num_voltage_sources,
            num_nodes + num_voltage_sources,
        );

        let mut b = DMatrix::zeros(num_nodes + num_voltage_sources, 1);

        self.netlist.get_resistors().iter().for_each(|c| {
            c.stamp(
                &mut ABMatrixView::new(&mut a, &mut b, num_nodes, c.num_variables(), num_nodes),
                dt,
            )
        });

        self.netlist.get_capacitors().iter().for_each(|c| {
            c.stamp(
                &mut ABMatrixView::new(&mut a, &mut b, num_nodes, c.num_variables(), num_nodes),
                dt,
            )
        });

        self.netlist.get_inductors().iter().for_each(|c| {
            c.stamp(
                &mut ABMatrixView::new(&mut a, &mut b, num_nodes, c.num_variables(), num_nodes),
                dt,
            )
        });

        self.netlist
            .get_voltages_sources()
            .iter()
            .enumerate()
            .for_each(|(i, c)| {
                c.stamp(
                    &mut ABMatrixView::new(
                        &mut a,
                        &mut b,
                        num_nodes,
                        c.num_variables(),
                        num_nodes + i,
                    ),
                    dt,
                )
            });

        self.netlist.get_current_sources().iter().for_each(|c| {
            c.stamp(
                &mut ABMatrixView::new(&mut a, &mut b, num_nodes, c.num_variables(), num_nodes),
                dt,
            )
        });

        let x = a.try_inverse().unwrap() * b;

        self.netlist.get_resistors_mut().iter_mut().for_each(|c| {
            c.update(
                &XMatrixView::new(&x, num_nodes, c.num_variables(), num_nodes),
                dt,
            )
        });

        self.netlist.get_capacitors_mut().iter_mut().for_each(|c| {
            c.update(
                &XMatrixView::new(&x, num_nodes, c.num_variables(), num_nodes),
                dt,
            )
        });

        self.netlist.get_inductors_mut().iter_mut().for_each(|c| {
            c.update(
                &XMatrixView::new(&x, num_nodes, c.num_variables(), num_nodes),
                dt,
            )
        });

        self.netlist
            .get_voltages_sources_mut()
            .iter_mut()
            .enumerate()
            .for_each(|(i, c)| {
                c.update(
                    &XMatrixView::new(&x, num_nodes, c.num_variables(), num_nodes + i),
                    dt,
                )
            });

        self.netlist
            .get_current_sources_mut()
            .iter_mut()
            .for_each(|c| {
                c.update(
                    &XMatrixView::new(&x, num_nodes, c.num_variables(), num_nodes),
                    dt,
                )
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

    #[test]
    fn test_voltage_source_inductor() {
        let mut netlist = Netlist::new();
        netlist
            .add_voltage_source(VoltageSource::new(1, 0, 1.0))
            .add_inductor(Inductor::new(1, 0, 0.5, 0.0));

        let mut solver = BESolver::new(&mut netlist);
        solver.solve(0.25);

        println!("{:?}", netlist);

        assert_relative_eq!(
            netlist.get_voltages_sources()[0].get_voltage(),
            1.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_voltages_sources()[0].get_current(),
            0.5,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_inductors()[0].get_voltage(),
            1.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_inductors()[0].get_current(),
            0.5,
            max_relative = 0.001
        );

        let mut solver = BESolver::new(&mut netlist);
        solver.solve(0.75);

        println!("{:?}", netlist);

        assert_relative_eq!(
            netlist.get_voltages_sources()[0].get_voltage(),
            1.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_voltages_sources()[0].get_current(),
            2.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_inductors()[0].get_voltage(),
            1.0,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_inductors()[0].get_current(),
            2.0,
            max_relative = 0.001
        );
    }

    #[test]
    fn test_rl_step() {
        let mut netlist = Netlist::new();
        netlist
            .add_voltage_source(VoltageSource::new(1, 0, 1.0))
            .add_resistor(Resistor::new(1, 2, 0.001))
            .add_inductor(Inductor::new(2, 0, 0.01, 0.0));

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
            95.162581964,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_resistors()[0].get_voltage(),
            0.095162581964,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_resistors()[0].get_current(),
            95.162581964,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_inductors()[0].get_voltage(),
            0.904837418036,
            max_relative = 0.001
        );
        assert_relative_eq!(
            netlist.get_inductors()[0].get_current(),
            95.162581964,
            max_relative = 0.001
        );
    }
}
