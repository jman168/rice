use std::fmt::Debug;

use crate::components::Component;

#[derive(Clone, Copy, PartialEq)]
pub struct CurrentSource {
    // Static variables
    positive_node: usize,
    negative_node: usize,
    current: f64,

    // Computed variables
    voltage: f64,
}

impl CurrentSource {
    pub fn new(positive_node: usize, negative_node: usize, current: f64) -> Self {
        Self {
            positive_node,
            negative_node,
            current,
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

    pub fn get_current(&self) -> f64 {
        self.current
    }

    pub fn get_voltage(&self) -> f64 {
        self.voltage
    }

    pub fn set_voltage(&mut self, voltage: f64) {
        self.voltage = voltage;
    }

    pub fn get_power(&self) -> f64 {
        self.get_voltage() * self.get_current()
    }
}

impl Debug for CurrentSource {
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

impl TryFrom<Component> for CurrentSource {
    type Error = ();

    fn try_from(value: Component) -> Result<Self, Self::Error> {
        match value {
            Component::CurrentSource(c) => Ok(c),
            _ => Err(()),
        }
    }
}
