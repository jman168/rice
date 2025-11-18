use std::fmt::Debug;

pub struct Inductor {
    // Static variables
    positive_node: usize,
    negative_node: usize,
    inductance: f64,

    // State variables
    current: f64,

    // Computed variables
    voltage: f64,
}

impl Inductor {
    pub fn new(
        positive_node: usize,
        negative_node: usize,
        inductance: f64,
        initial_current: f64,
    ) -> Self {
        Self {
            positive_node,
            negative_node,
            inductance,
            current: initial_current,
            voltage: 0.0,
        }
    }

    pub fn get_positive_node(&self) -> usize {
        self.positive_node
    }

    pub fn get_negative_node(&self) -> usize {
        self.negative_node
    }

    pub fn get_inductance(&self) -> f64 {
        self.inductance
    }

    pub fn get_current(&self) -> f64 {
        self.current
    }

    pub fn set_current(&mut self, current: f64) {
        self.current = current;
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

impl Debug for Inductor {
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
