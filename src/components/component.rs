use crate::components::{Capacitor, CurrentSource, Inductor, Resistor, VoltageSource};

#[derive(Debug)]
pub enum Component {
    Resistor(Resistor),
    Capacitor(Capacitor),
    Inductor(Inductor),
    VoltageSource(VoltageSource),
    CurrentSource(CurrentSource),
}

impl Component {
    pub fn max_node(&self) -> usize {
        match self {
            Self::Resistor(c) => c.max_node(),
            Self::Capacitor(c) => c.max_node(),
            Self::Inductor(c) => c.max_node(),
            Self::VoltageSource(c) => c.max_node(),
            Self::CurrentSource(c) => c.max_node(),
        }
    }

    pub fn get_voltage(&self) -> f64 {
        match self {
            Self::Resistor(c) => c.get_voltage(),
            Self::Capacitor(c) => c.get_voltage(),
            Self::Inductor(c) => c.get_voltage(),
            Self::VoltageSource(c) => c.get_voltage(),
            Self::CurrentSource(c) => c.get_voltage(),
        }
    }

    pub fn get_current(&self) -> f64 {
        match self {
            Self::Resistor(c) => c.get_current(),
            Self::Capacitor(c) => c.get_current(),
            Self::Inductor(c) => c.get_current(),
            Self::VoltageSource(c) => c.get_current(),
            Self::CurrentSource(c) => c.get_current(),
        }
    }

    pub fn get_power(&self) -> f64 {
        match self {
            Self::Resistor(c) => c.get_power(),
            Self::Capacitor(c) => c.get_power(),
            Self::Inductor(c) => c.get_power(),
            Self::VoltageSource(c) => c.get_power(),
            Self::CurrentSource(c) => c.get_power(),
        }
    }
}

impl From<Resistor> for Component {
    fn from(value: Resistor) -> Self {
        Self::Resistor(value)
    }
}

impl From<Capacitor> for Component {
    fn from(value: Capacitor) -> Self {
        Self::Capacitor(value)
    }
}

impl From<Inductor> for Component {
    fn from(value: Inductor) -> Self {
        Self::Inductor(value)
    }
}

impl From<VoltageSource> for Component {
    fn from(value: VoltageSource) -> Self {
        Self::VoltageSource(value)
    }
}

impl From<CurrentSource> for Component {
    fn from(value: CurrentSource) -> Self {
        Self::CurrentSource(value)
    }
}
