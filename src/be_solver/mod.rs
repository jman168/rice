use nalgebra::DMatrix;

use crate::components::{Netlist, Resistor, VoltageSource};

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
        let dim = num_nodes + num_voltage_sources;

        // Create the matrix which we will later stamp and then solve.
        let mut a: DMatrix<f64> = DMatrix::zeros(dim, dim);
        // Create the B matrix which we will later use to get the solution matrix x.
        let mut b: DMatrix<f64> = DMatrix::zeros(dim, 1);

        self.netlist
            .get_resistors()
            .iter()
            .enumerate()
            .for_each(|(i, r)| Self::stamp_resistor(r, i, &mut a, &mut b, num_nodes));

        self.netlist
            .get_voltages_sources()
            .iter()
            .enumerate()
            .for_each(|(i, v)| Self::stamp_voltage_source(v, i, &mut a, &mut b, num_nodes));

        let x = a.try_inverse().unwrap() * b;

        self.netlist
            .get_resistors_mut()
            .iter_mut()
            .enumerate()
            .for_each(|(i, r)| Self::update_resistor(r, i, &x, num_nodes));

        self.netlist
            .get_voltages_sources_mut()
            .iter_mut()
            .enumerate()
            .for_each(|(i, v)| Self::update_voltage_source(v, i, &x, num_nodes));
    }

    fn stamp_resistor(
        resistor: &Resistor,
        _resistor_index: usize,
        a: &mut DMatrix<f64>,
        _b: &mut DMatrix<f64>,
        _num_nodes: usize,
    ) {
        let p = resistor.get_positive_node();
        let n = resistor.get_negative_node();
        let g = 1.0 / resistor.get_resistance();

        if p != 0 {
            *a.get_mut((p - 1, p - 1)).unwrap() -= g;
            if n != 0 {
                *a.get_mut((p - 1, n - 1)).unwrap() += g;
            }
        }

        if n != 0 {
            *a.get_mut((n - 1, n - 1)).unwrap() -= g;
            if p != 0 {
                *a.get_mut((n - 1, p - 1)).unwrap() += g;
            }
        }
    }

    fn stamp_voltage_source(
        voltage_source: &VoltageSource,
        voltage_source_index: usize,
        a: &mut DMatrix<f64>,
        b: &mut DMatrix<f64>,
        num_nodes: usize,
    ) {
        let p = voltage_source.get_positive_node();
        let n = voltage_source.get_negative_node();

        if p != 0 {
            *a.get_mut((p - 1, num_nodes + voltage_source_index))
                .unwrap() = 1.0;
            *a.get_mut((num_nodes + voltage_source_index, p - 1))
                .unwrap() = 1.0;
        }

        if n != 0 {
            *a.get_mut((n - 1, num_nodes + voltage_source_index))
                .unwrap() = -1.0;
            *a.get_mut((num_nodes + voltage_source_index, n - 1))
                .unwrap() = -1.0;
        }

        *b.get_mut((num_nodes + voltage_source_index, 0)).unwrap() = voltage_source.get_voltage();
    }

    fn update_resistor(
        resistor: &mut Resistor,
        _resistor_index: usize,
        x: &DMatrix<f64>,
        _num_nodes: usize,
    ) {
        let p = resistor.get_positive_node();
        let n = resistor.get_negative_node();

        let vp = if p != 0 {
            *x.get((p - 1, 0)).unwrap()
        } else {
            0.0
        };

        let vn = if n != 0 {
            *x.get((n - 1, 0)).unwrap()
        } else {
            0.0
        };

        resistor.set_voltage(vp - vn);
    }

    fn update_voltage_source(
        voltage_source: &mut VoltageSource,
        voltage_source_index: usize,
        x: &DMatrix<f64>,
        num_nodes: usize,
    ) {
        voltage_source.set_current(*x.get((num_nodes + voltage_source_index, 0)).unwrap());
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
