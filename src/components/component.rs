use crate::components::{Capacitor, CurrentSource, Diode, Inductor, Resistor, VoltageSource};

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Component {
    Resistor(Resistor),
    Capacitor(Capacitor),
    Inductor(Inductor),
    Diode(Diode),
    VoltageSource(VoltageSource),
    CurrentSource(CurrentSource),
}

impl Component {
    pub fn max_node(&self) -> usize {
        match self {
            Self::Resistor(c) => c.max_node(),
            Self::Capacitor(c) => c.max_node(),
            Self::Inductor(c) => c.max_node(),
            Self::Diode(c) => c.max_node(),
            Self::VoltageSource(c) => c.max_node(),
            Self::CurrentSource(c) => c.max_node(),
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

impl From<Diode> for Component {
    fn from(value: Diode) -> Self {
        Self::Diode(value)
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
