use std::fmt::Debug;

use crate::components::Component;

#[derive(Clone, Copy, PartialEq)]
pub struct Diode {
    // Static variables
    positive_node: usize,
    negative_node: usize,
    reverse_saturation_current: f64,
    ideality_factor: f64,
    thermal_voltage: f64,

    // Computed variables
    voltage: f64,
}

impl Diode {
    pub fn new(
        positive_node: usize,
        negative_node: usize,
        reverse_saturation_current: f64,
        ideality_factor: f64,
        thermal_voltage: f64,
    ) -> Self {
        Self {
            positive_node,
            negative_node,
            reverse_saturation_current,
            ideality_factor,
            thermal_voltage,
            voltage: 0.0,
        }
    }

    pub fn max_node(&self) -> usize {
        self.get_positive_node().max(self.get_negative_node())
    }

    pub fn get_positive_node(&self) -> usize {
        self.positive_node
    }

    pub fn get_negative_node(&self) -> usize {
        self.negative_node
    }

    pub fn get_reverse_saturation_current(&self) -> f64 {
        self.reverse_saturation_current
    }

    pub fn get_ideality_factor(&self) -> f64 {
        self.ideality_factor
    }

    pub fn get_thermal_voltage(&self) -> f64 {
        self.thermal_voltage
    }

    pub fn get_voltage(&self) -> f64 {
        self.voltage
    }

    pub fn set_voltage(&mut self, voltage: f64) {
        self.voltage = voltage;
    }

    pub fn get_current(&self) -> f64 {
        self.get_current_with_voltage(self.voltage)
    }

    pub fn get_power(&self) -> f64 {
        self.get_voltage() * self.get_current()
    }

    pub fn get_current_with_voltage(&self, voltage: f64) -> f64 {
        self.reverse_saturation_current
            * ((voltage / (self.ideality_factor * self.thermal_voltage)).exp() - 1.0)
    }

    pub fn get_di_dv(&self, voltage: f64) -> f64 {
        (self.reverse_saturation_current / (self.ideality_factor * self.thermal_voltage))
            * (voltage / (self.ideality_factor * self.thermal_voltage)).exp()
    }
}

impl Debug for Diode {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "{{v: {}, i: {}, p: {}}}",
            self.get_voltage(),
            self.get_current(),
            self.get_power()
        )
    }
}

impl TryFrom<Component> for Diode {
    type Error = ();

    fn try_from(value: Component) -> Result<Self, Self::Error> {
        match value {
            Component::Diode(c) => Ok(c),
            _ => Err(()),
        }
    }
}
