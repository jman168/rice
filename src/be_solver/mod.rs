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

        // Create matrices for the Ax=B linear system of equations.
        let mut a = DMatrix::zeros(num_nodes + num_variables, num_nodes + num_variables);
        let mut b = DMatrix::zeros(num_nodes + num_variables, 1);
        let mut x = DMatrix::zeros(num_nodes + num_variables, 1);

        // Variable for keeping track of how many iterations we have completed so far.
        let mut iteration = 0;

        loop {
            // Bound the maximum number of iterations so we don't accidentally run forever.
            if iteration >= 1000 {
                panic!(
                    "failed to solve in less than {} iterations {:?}",
                    iteration, x
                );
            }

            // Reset the A and B matrices so they can be stamped in the next step.
            a.fill(0.0);
            b.fill(0.0);

            // Stamp the matrices with all the components equations.
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
                    let op_point =
                        XMatrixView::new(&x, num_nodes, c.num_variables(), variables_start);
                    c.stamp(&mut view, &op_point, dt);
                    variables_start + c.num_variables()
                });

            // Try to invert the matrix in place.
            if !a.try_inverse_mut() {
                panic!("failed to invert problem matrix!");
            }

            // Compute the new solution to the equation.
            let x_prime = &a * &b;

            // Compute the tolerance for each solved parameter (fractional difference between the new iteration and old
            // iteration).
            let mut x_t = &x_prime - &x;
            x_t.component_div_assign(&x_prime);
            // We do this instead of .abs() as nalgebra has no native in place abs operation.
            x_t.iter_mut().for_each(|x| *x = x.abs());

            // Compute the overall tolerance of the entire matrix (worse case tolerance).
            let tolerance = x_t.max();

            // Update the solution matrix for later use.
            x = x_prime;

            // If we are happy with the solution, we are done!
            if tolerance < 1e-4 {
                break;
            }
            // If we are not happy with the solution yet, try again.
            else {
                iteration += 1;
            }
        }

        // Now that we are happy with the given solution, update all the components to reflect it.
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
        components::{Capacitor, CurrentSource, Diode, Inductor, Netlist, Resistor, VoltageSource},
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

        let v: VoltageSource = netlist.get_components()[0].try_into().unwrap();
        let r: Resistor = netlist.get_components()[1].try_into().unwrap();

        assert_relative_eq!(v.get_voltage(), 10.0, max_relative = 0.001);
        assert_relative_eq!(v.get_current(), 5.0, max_relative = 0.001);
        assert_relative_eq!(r.get_voltage(), 10.0, max_relative = 0.001);
        assert_relative_eq!(r.get_current(), 5.0, max_relative = 0.001);
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

        let v1: VoltageSource = netlist.get_components()[0].try_into().unwrap();
        let r: Resistor = netlist.get_components()[1].try_into().unwrap();
        let v2: VoltageSource = netlist.get_components()[2].try_into().unwrap();

        assert_relative_eq!(v1.get_voltage(), 10.0, max_relative = 0.001);
        assert_relative_eq!(v1.get_current(), 2.5, max_relative = 0.001);
        assert_relative_eq!(r.get_voltage(), 5.0, max_relative = 0.001);
        assert_relative_eq!(r.get_current(), 2.5, max_relative = 0.001);
        assert_relative_eq!(v2.get_voltage(), 5.0, max_relative = 0.001);
        assert_relative_eq!(v2.get_current(), -2.5, max_relative = 0.001);
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

        let v: VoltageSource = netlist.get_components()[0].try_into().unwrap();
        let r1: Resistor = netlist.get_components()[1].try_into().unwrap();
        let r2: Resistor = netlist.get_components()[2].try_into().unwrap();

        assert_relative_eq!(v.get_voltage(), 5.0, max_relative = 0.001);
        assert_relative_eq!(v.get_current(), 1.0, max_relative = 0.001);
        assert_relative_eq!(r1.get_voltage(), 4.0, max_relative = 0.001);
        assert_relative_eq!(r1.get_current(), 1.0, max_relative = 0.001);
        assert_relative_eq!(r2.get_voltage(), 1.0, max_relative = 0.001);
        assert_relative_eq!(r2.get_current(), 1.0, max_relative = 0.001);
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

        let i: CurrentSource = netlist.get_components()[0].try_into().unwrap();
        let r: Resistor = netlist.get_components()[1].try_into().unwrap();

        assert_relative_eq!(i.get_voltage(), 10.0, max_relative = 0.001);
        assert_relative_eq!(i.get_current(), 5.0, max_relative = 0.001);
        assert_relative_eq!(r.get_voltage(), 10.0, max_relative = 0.001);
        assert_relative_eq!(r.get_current(), 5.0, max_relative = 0.001);
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

        let i: CurrentSource = netlist.get_components()[0].try_into().unwrap();
        let c: Capacitor = netlist.get_components()[1].try_into().unwrap();

        assert_relative_eq!(i.get_voltage(), 0.50, max_relative = 0.001);
        assert_relative_eq!(i.get_current(), 1.0, max_relative = 0.001);
        assert_relative_eq!(c.get_voltage(), 0.50, max_relative = 0.001);
        assert_relative_eq!(c.get_current(), 1.0, max_relative = 0.001);

        let mut solver = BESolver::new(&mut netlist);
        solver.solve(0.75);

        println!("{:?}", netlist);

        let i: CurrentSource = netlist.get_components()[0].try_into().unwrap();
        let c: Capacitor = netlist.get_components()[1].try_into().unwrap();

        assert_relative_eq!(i.get_voltage(), 2.0, max_relative = 0.001);
        assert_relative_eq!(i.get_current(), 1.0, max_relative = 0.001);
        assert_relative_eq!(c.get_voltage(), 2.0, max_relative = 0.001);
        assert_relative_eq!(c.get_current(), 1.0, max_relative = 0.001);
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

        let v: VoltageSource = netlist.get_components()[0].try_into().unwrap();
        let r: Resistor = netlist.get_components()[1].try_into().unwrap();
        let c: Capacitor = netlist.get_components()[2].try_into().unwrap();

        assert_relative_eq!(v.get_voltage(), 1.0, max_relative = 0.001);
        assert_relative_eq!(v.get_current(), 0.000367879441171, max_relative = 0.001);
        assert_relative_eq!(r.get_voltage(), 0.367879441171, max_relative = 0.001);
        assert_relative_eq!(r.get_current(), 0.000367879441171, max_relative = 0.001);
        assert_relative_eq!(c.get_voltage(), 0.632120558829, max_relative = 0.001);
        assert_relative_eq!(c.get_current(), 0.000367879441171, max_relative = 0.001);
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

        let v: VoltageSource = netlist.get_components()[0].try_into().unwrap();
        let l: Inductor = netlist.get_components()[1].try_into().unwrap();

        assert_relative_eq!(v.get_voltage(), 1.0, max_relative = 0.001);
        assert_relative_eq!(v.get_current(), 0.5, max_relative = 0.001);
        assert_relative_eq!(l.get_voltage(), 1.0, max_relative = 0.001);
        assert_relative_eq!(l.get_current(), 0.5, max_relative = 0.001);

        let mut solver = BESolver::new(&mut netlist);
        solver.solve(0.75);

        println!("{:?}", netlist);

        let v: VoltageSource = netlist.get_components()[0].try_into().unwrap();
        let l: Inductor = netlist.get_components()[1].try_into().unwrap();

        assert_relative_eq!(v.get_voltage(), 1.0, max_relative = 0.001);
        assert_relative_eq!(v.get_current(), 2.0, max_relative = 0.001);
        assert_relative_eq!(l.get_voltage(), 1.0, max_relative = 0.001);
        assert_relative_eq!(l.get_current(), 2.0, max_relative = 0.001);
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

        let v: VoltageSource = netlist.get_components()[0].try_into().unwrap();
        let r: Resistor = netlist.get_components()[1].try_into().unwrap();
        let l: Inductor = netlist.get_components()[2].try_into().unwrap();

        assert_relative_eq!(v.get_voltage(), 1.0, max_relative = 0.001);
        assert_relative_eq!(v.get_current(), 95.162581964, max_relative = 0.001);
        assert_relative_eq!(r.get_voltage(), 0.095162581964, max_relative = 0.001);
        assert_relative_eq!(r.get_current(), 95.162581964, max_relative = 0.001);
        assert_relative_eq!(l.get_voltage(), 0.904837418036, max_relative = 0.001);
        assert_relative_eq!(l.get_current(), 95.162581964, max_relative = 0.001);
    }

    #[test]
    fn test_r_diode() {
        let mut netlist = Netlist::new();
        netlist
            .add_component(VoltageSource::new(1, 0, 10.0))
            .add_component(Resistor::new(1, 2, 10.0))
            // Approximately a Vishay 1N4148 just for fun
            .add_component(Diode::new(2, 0, 25.0e-9, 2.0, 26.0e-3));

        let mut solver = BESolver::new(&mut netlist);
        solver.solve(0.001);

        println!("{:?}", netlist);

        let v: VoltageSource = netlist.get_components()[0].try_into().unwrap();
        let r: Resistor = netlist.get_components()[1].try_into().unwrap();
        let d: Diode = netlist.get_components()[2].try_into().unwrap();

        assert_relative_eq!(v.get_voltage(), 10.0, max_relative = 0.001);
        assert_relative_eq!(v.get_current(), 0.9094706, max_relative = 0.001);
        assert_relative_eq!(r.get_voltage(), 9.094706, max_relative = 0.001);
        assert_relative_eq!(r.get_current(), 0.9094706, max_relative = 0.001);
        assert_relative_eq!(d.get_voltage(), 0.905294, max_relative = 0.001);
        assert_relative_eq!(d.get_current(), 0.9094706, max_relative = 0.001);
    }
}
