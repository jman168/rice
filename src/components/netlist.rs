use crate::components::Component;

#[derive(Debug)]
pub struct Netlist {
    components: Vec<Component>,
}

impl Netlist {
    pub fn new() -> Self {
        Self {
            components: Vec::new(),
        }
    }

    /// Adds a single component to the netlist.
    pub fn add_component(&mut self, component: impl Into<Component>) -> &mut Self {
        self.components.push(component.into());
        self
    }

    /// Adds multiple components to the netlist.
    pub fn add_components(
        &mut self,
        components: impl Iterator<Item = impl Into<Component>>,
    ) -> &mut Self {
        self.components.extend(components.map(|c| c.into()));
        self
    }

    /// Gets all the components in the netlist in the order they were added.
    pub fn get_components(&self) -> &Vec<Component> {
        &self.components
    }

    /// Gets mutatable references to all the components in the netlist in the order they were
    /// added.
    pub fn get_components_mut(&mut self) -> &mut Vec<Component> {
        &mut self.components
    }

    pub fn get_num_nodes(&self) -> usize {
        self.components
            .iter()
            .map(|c| c.max_node())
            .max()
            .unwrap_or(0)
    }
}

impl Default for Netlist {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::components::{Resistor, VoltageSource};

    #[test]
    fn test_get_num_nodes() {
        let mut netlist = Netlist::new();
        netlist
            .add_component(Resistor::new(1, 2, 1.0))
            .add_component(VoltageSource::new(3, 4, 1.0));
        assert_eq!(netlist.get_num_nodes(), 4);

        let mut netlist = Netlist::new();
        netlist
            .add_component(VoltageSource::new(1, 2, 1.0))
            .add_component(Resistor::new(3, 4, 1.0));
        assert_eq!(netlist.get_num_nodes(), 4);
    }
}
