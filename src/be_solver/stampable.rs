use crate::{
    be_solver::matrix_view::{ABMatrixView, ViewEquationIndex, ViewVariableIndex, XMatrixView},
    components::{Capacitor, Component, CurrentSource, Inductor, Resistor, VoltageSource},
};

pub trait Stampable {
    /// Returns the number of additional variables this component will add to the matrix.
    fn num_variables(&self) -> usize;

    /// Stamps the coefficients of the component.
    ///
    /// If this component is non-linear (needing Newton-Raphson iteration), it may choose to
    /// compute its partials given the op_point XMatrixView. This is computed as either an initial
    /// guess, or the solution to the last solve.
    fn stamp(&self, view: &mut ABMatrixView, op_point: &XMatrixView, dt: f64);

    /// Updates the component state based on the given solution.
    fn update(&mut self, view: &XMatrixView, dt: f64);
}

impl Stampable for Resistor {
    fn num_variables(&self) -> usize {
        0
    }

    fn stamp(&self, view: &mut ABMatrixView, _op_point: &XMatrixView, _dt: f64) {
        let positive_equation_index = ViewEquationIndex::NodalEquation(self.get_positive_node());
        let negative_equation_index = ViewEquationIndex::NodalEquation(self.get_negative_node());

        let positive_voltage_index = ViewVariableIndex::NodeVoltage(self.get_positive_node());
        let negative_voltage_index = ViewVariableIndex::NodeVoltage(self.get_negative_node());

        // Compute resistor conductance
        let g = 1.0 / self.get_resistance();

        // Current flowing out of positive node is (v_positive - v_negative) / R
        view.coefficient_add(positive_equation_index, positive_voltage_index, g);
        view.coefficient_add(positive_equation_index, negative_voltage_index, -g);

        // Current flowing out of positive node is -(v_positive - v_negative) / R
        view.coefficient_add(negative_equation_index, positive_voltage_index, -g);
        view.coefficient_add(negative_equation_index, negative_voltage_index, g);
    }

    fn update(&mut self, view: &XMatrixView, _dt: f64) {
        let positive_voltage_index = ViewVariableIndex::NodeVoltage(self.get_positive_node());
        let negative_voltage_index = ViewVariableIndex::NodeVoltage(self.get_negative_node());

        self.set_voltage(
            view.get_variable(positive_voltage_index).unwrap()
                - view.get_variable(negative_voltage_index).unwrap(),
        );
    }
}

impl Stampable for Capacitor {
    fn num_variables(&self) -> usize {
        0
    }

    fn stamp(&self, view: &mut ABMatrixView, _op_point: &XMatrixView, dt: f64) {
        let positive_equation_index = ViewEquationIndex::NodalEquation(self.get_positive_node());
        let negative_equation_index = ViewEquationIndex::NodalEquation(self.get_negative_node());

        let positive_voltage_index = ViewVariableIndex::NodeVoltage(self.get_positive_node());
        let negative_voltage_index = ViewVariableIndex::NodeVoltage(self.get_negative_node());

        let c = self.get_capacitance();

        // The differential equation describing a capacitor is i = C*dv/dt.
        // Discretizing we get i = C*(v_new - v_old)/dt.
        // Further expanding this we get i = C*(v_positive - v_negative)/dt - C*v_old/dt.

        // Current flowing out of the positive node is C*v_positive/dt - C*v_negative/dt - C*v_old/dt.
        view.coefficient_add(positive_equation_index, positive_voltage_index, c / dt);
        view.coefficient_add(positive_equation_index, negative_voltage_index, -c / dt);
        view.result_add(positive_equation_index, c * self.get_voltage() / dt);

        // Current flowing out of the negative node is -C*v_positive/dt + C*v_negative/dt + C*v_old/dt.
        view.coefficient_add(negative_equation_index, positive_voltage_index, -c / dt);
        view.coefficient_add(negative_equation_index, negative_voltage_index, c / dt);
        view.result_add(negative_equation_index, -c * self.get_voltage() / dt);
    }

    fn update(&mut self, view: &XMatrixView, dt: f64) {
        let positive_voltage_index = ViewVariableIndex::NodeVoltage(self.get_positive_node());
        let negative_voltage_index = ViewVariableIndex::NodeVoltage(self.get_negative_node());

        let new_voltage = view.get_variable(positive_voltage_index).unwrap()
            - view.get_variable(negative_voltage_index).unwrap();

        // Discretized equation is i = C*(v_new - v_old)/dt (see capacitor stamping function).

        self.set_current(self.get_capacitance() * (new_voltage - self.get_voltage()) / dt);

        self.set_voltage(new_voltage);
    }
}

impl Stampable for Inductor {
    fn num_variables(&self) -> usize {
        0
    }

    fn stamp(&self, view: &mut ABMatrixView, _op_point: &XMatrixView, dt: f64) {
        let positive_equation_index = ViewEquationIndex::NodalEquation(self.get_positive_node());
        let negative_equation_index = ViewEquationIndex::NodalEquation(self.get_negative_node());

        let positive_voltage_index = ViewVariableIndex::NodeVoltage(self.get_positive_node());
        let negative_voltage_index = ViewVariableIndex::NodeVoltage(self.get_negative_node());

        let l = self.get_inductance();

        // The differential equation describing an inductor is v = L*di/dt.
        // Discretizing we get v = L*(i_new - i_old)/dt.
        // Further expanding this we get v_positve - v_negative = L*(i_new - i_old)/dt.
        // Doing some algebra to solver for i_new we get:
        // i_new = v_positive*dt/L - v_negative*dt/L + i_old.

        // Current flowing out of the positive node is v_positive*dt/L - v_negative*dt/L + i_old.
        view.coefficient_add(positive_equation_index, positive_voltage_index, dt / l);
        view.coefficient_add(positive_equation_index, negative_voltage_index, -dt / l);
        view.result_add(positive_equation_index, -self.get_current());

        // Current flowing out of the negative node is -v_positive*dt/L + v_negative*dt/L - i_old.
        view.coefficient_add(negative_equation_index, positive_voltage_index, -dt / l);
        view.coefficient_add(negative_equation_index, negative_voltage_index, dt / l);
        view.result_add(negative_equation_index, self.get_current());
    }

    fn update(&mut self, view: &XMatrixView, dt: f64) {
        let positive_voltage_index = ViewVariableIndex::NodeVoltage(self.get_positive_node());
        let negative_voltage_index = ViewVariableIndex::NodeVoltage(self.get_negative_node());

        self.set_voltage(
            view.get_variable(positive_voltage_index).unwrap()
                - view.get_variable(negative_voltage_index).unwrap(),
        );

        // Discretized equation is i_new = v_positive*dt/L - v_negative*dt/L + i_old (see inductor stamping function).

        self.set_current(self.get_voltage() * dt / self.get_inductance() + self.get_current());
    }
}

impl Stampable for VoltageSource {
    fn num_variables(&self) -> usize {
        1
    }

    fn stamp(&self, view: &mut ABMatrixView, _op_point: &XMatrixView, _dt: f64) {
        let positive_equation_index = ViewEquationIndex::NodalEquation(self.get_positive_node());
        let negative_equation_index = ViewEquationIndex::NodalEquation(self.get_negative_node());
        let specific_equation_index = ViewEquationIndex::SpecificEquation(0);

        let positive_voltage_index = ViewVariableIndex::NodeVoltage(self.get_positive_node());
        let negative_voltage_index = ViewVariableIndex::NodeVoltage(self.get_negative_node());
        let current_index = ViewVariableIndex::SpecificVariable(0);

        // Current flowing out of positive node is -i_source
        view.coefficient_add(positive_equation_index, current_index, -1.0);
        // Current flowing out of negative node is i_source
        view.coefficient_add(negative_equation_index, current_index, 1.0);

        // Source equation is v_positive - v_negative = v_source
        view.coefficient_add(specific_equation_index, positive_voltage_index, 1.0);
        view.coefficient_add(specific_equation_index, negative_voltage_index, -1.0);
        view.result_add(specific_equation_index, self.get_voltage());
    }

    fn update(&mut self, view: &XMatrixView, _dt: f64) {
        let current_index = ViewVariableIndex::SpecificVariable(0);
        self.set_current(view.get_variable(current_index).unwrap());
    }
}

impl Stampable for CurrentSource {
    fn num_variables(&self) -> usize {
        0
    }

    fn stamp(&self, view: &mut ABMatrixView, _op_point: &XMatrixView, _dt: f64) {
        let positive_equation_index = ViewEquationIndex::NodalEquation(self.get_positive_node());
        let negative_equation_index = ViewEquationIndex::NodalEquation(self.get_negative_node());

        // NOTE: the signs are flipped here because they take the form of constants, not
        // coefficients.

        // Current flowing out of positive node is -i_source
        view.result_add(positive_equation_index, self.get_current());
        // Current flowing out of negative node is i_source
        view.result_add(negative_equation_index, -self.get_current());
    }

    fn update(&mut self, view: &XMatrixView, _dt: f64) {
        let positive_voltage_index = ViewVariableIndex::NodeVoltage(self.get_positive_node());
        let negative_voltage_index = ViewVariableIndex::NodeVoltage(self.get_negative_node());

        self.set_voltage(
            view.get_variable(positive_voltage_index).unwrap()
                - view.get_variable(negative_voltage_index).unwrap(),
        );
    }
}

impl Stampable for Component {
    fn num_variables(&self) -> usize {
        match self {
            Self::Resistor(c) => c.num_variables(),
            Self::Capacitor(c) => c.num_variables(),
            Self::Inductor(c) => c.num_variables(),
            Self::VoltageSource(c) => c.num_variables(),
            Self::CurrentSource(c) => c.num_variables(),
        }
    }

    fn stamp(&self, view: &mut ABMatrixView, op_point: &XMatrixView, dt: f64) {
        match self {
            Self::Resistor(c) => c.stamp(view, op_point, dt),
            Self::Capacitor(c) => c.stamp(view, op_point, dt),
            Self::Inductor(c) => c.stamp(view, op_point, dt),
            Self::VoltageSource(c) => c.stamp(view, op_point, dt),
            Self::CurrentSource(c) => c.stamp(view, op_point, dt),
        }
    }

    fn update(&mut self, view: &XMatrixView, dt: f64) {
        match self {
            Self::Resistor(c) => c.update(view, dt),
            Self::Capacitor(c) => c.update(view, dt),
            Self::Inductor(c) => c.update(view, dt),
            Self::VoltageSource(c) => c.update(view, dt),
            Self::CurrentSource(c) => c.update(view, dt),
        }
    }
}
