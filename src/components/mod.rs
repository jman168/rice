use std::fmt::Debug;

mod resistor;
pub use resistor::Resistor;

mod voltage_source;
pub use voltage_source::VoltageSource;

/// A struct that represents a collection of electrical components connected to each other.
pub struct Netlist {
    resistors: Vec<Resistor>,
    voltage_sources: Vec<VoltageSource>,
}

impl Netlist {
    /// Creates a new empty netlist.
    pub fn new() -> Self {
        Self {
            resistors: Vec::new(),
            voltage_sources: Vec::new(),
        }
    }

    /// Adds a single resistor to the netlist.
    pub fn add_resistor(&mut self, resistor: Resistor) -> &mut Self {
        self.resistors.push(resistor);
        self
    }

    /// Adds multiple resistors to the netlist.
    pub fn add_resistors(&mut self, resistors: impl Iterator<Item = Resistor>) -> &mut Self {
        self.resistors.extend(resistors);
        self
    }

    /// Gets all the resistors in the netlist in the order they were added.
    pub fn get_resistors(&self) -> &Vec<Resistor> {
        &self.resistors
    }

    /// Gets mutatable reference to all the resistors in the netlist in the order they were added.
    pub fn get_resistors_mut(&mut self) -> &mut Vec<Resistor> {
        &mut self.resistors
    }

    /// Adds a single voltage source to the netlist.
    pub fn add_voltage_source(&mut self, voltage_source: VoltageSource) -> &mut Self {
        self.voltage_sources.push(voltage_source);
        self
    }

    /// Adds multiple voltage sources to the netlist.
    pub fn add_voltage_sources(
        &mut self,
        voltage_sources: impl Iterator<Item = VoltageSource>,
    ) -> &mut Self {
        self.voltage_sources.extend(voltage_sources);
        self
    }

    /// Gets all the voltage sources in the netlist in the order they were added.
    pub fn get_voltages_sources(&self) -> &Vec<VoltageSource> {
        &self.voltage_sources
    }

    /// Gets mutatable references to all the voltage sources in the netlist in the order they were
    /// added.
    pub fn get_voltages_sources_mut(&mut self) -> &mut Vec<VoltageSource> {
        &mut self.voltage_sources
    }

    /// Gets the total number of nodes in the netlist (maximum node number referenced by all
    /// components.
    pub fn get_num_nodes(&self) -> usize {
        let max_r_node = self
            .resistors
            .iter()
            .map(|r| r.get_positive_node().max(r.get_negative_node()))
            .max()
            .unwrap_or(0);
        let max_v_node = self
            .voltage_sources
            .iter()
            .map(|v| v.get_positive_node().max(v.get_negative_node()))
            .max()
            .unwrap_or(0);

        max_r_node.max(max_v_node)
    }
}

impl Debug for Netlist {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(f, "Netlist:\n")?;
        write!(f, "\tResistors: [\n")?;
        for r in &self.resistors {
            write!(f, "\t\t{:?},\n", r)?;
        }
        write!(f, "\t]\n")?;

        write!(f, "\tVoltage Sources: [\n")?;
        for v in &self.voltage_sources {
            write!(f, "\t\t{:?},\n", v)?;
        }
        write!(f, "\t]\n")
    }
}

#[cfg(test)]
mod test {
    use super::*;

    #[test]
    fn test_get_num_nodes() {
        let mut netlist = Netlist::new();
        netlist
            .add_resistor(Resistor::new(1, 2, 1.0))
            .add_voltage_source(VoltageSource::new(3, 4, 1.0));
        assert_eq!(netlist.get_num_nodes(), 4);

        let mut netlist = Netlist::new();
        netlist
            .add_voltage_source(VoltageSource::new(1, 2, 1.0))
            .add_resistor(Resistor::new(3, 4, 1.0));
        assert_eq!(netlist.get_num_nodes(), 4);
    }
}
