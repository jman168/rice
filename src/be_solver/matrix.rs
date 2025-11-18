use nalgebra::DMatrix;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum EquationIndex {
    NodalEquation(usize),
    VoltageSourceEquation(usize),
}

impl EquationIndex {
    fn into_row(self, num_nodes: usize, num_voltage_sources: usize) -> Option<usize> {
        match self {
            EquationIndex::NodalEquation(idx) => {
                if idx == 0 {
                    return None;
                }

                if idx > num_nodes {
                    return None;
                }

                Some(idx - 1)
            }
            EquationIndex::VoltageSourceEquation(idx) => {
                if idx >= num_voltage_sources {
                    return None;
                }

                Some(num_nodes + idx)
            }
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VariableIndex {
    NodalVoltage(usize),
    VoltageSourceCurrent(usize),
}

impl VariableIndex {
    fn into_col(self, num_nodes: usize, num_voltage_sources: usize) -> Option<usize> {
        match self {
            VariableIndex::NodalVoltage(idx) => {
                if idx == 0 {
                    return None;
                }

                if idx > num_nodes {
                    return None;
                }

                Some(idx - 1)
            }
            VariableIndex::VoltageSourceCurrent(idx) => {
                if idx >= num_voltage_sources {
                    return None;
                }

                Some(num_nodes + idx)
            }
        }
    }
}

pub struct XMatrix {
    x: DMatrix<f64>,
    num_nodes: usize,
    num_voltage_sources: usize,
}

impl XMatrix {
    pub fn new(x: DMatrix<f64>, num_nodes: usize, num_voltage_sources: usize) -> Self {
        Self {
            x,
            num_nodes,
            num_voltage_sources,
        }
    }

    pub fn get_variable(&self, variable_index: VariableIndex) -> Option<f64> {
        match variable_index {
            VariableIndex::NodalVoltage(0) => Some(0.0),
            _ => self
                .x
                .get((
                    variable_index.into_col(self.num_nodes, self.num_voltage_sources)?,
                    0,
                ))
                .cloned(),
        }
    }
}

pub struct ABMatrix {
    a: DMatrix<f64>,
    b: DMatrix<f64>,
    num_nodes: usize,
    num_voltage_sources: usize,
}

impl ABMatrix {
    pub fn new(num_nodes: usize, num_voltage_sources: usize) -> Self {
        Self {
            a: DMatrix::zeros(
                num_nodes + num_voltage_sources,
                num_nodes + num_voltage_sources,
            ),
            b: DMatrix::zeros(num_nodes + num_voltage_sources, 1),
            num_nodes,
            num_voltage_sources,
        }
    }

    pub fn get_coefficient_mut(
        &mut self,
        equation_index: EquationIndex,
        variable_index: VariableIndex,
    ) -> Option<&mut f64> {
        self.a.get_mut((
            equation_index.into_row(self.num_nodes, self.num_voltage_sources)?,
            variable_index.into_col(self.num_nodes, self.num_voltage_sources)?,
        ))
    }

    pub fn get_result_mut(&mut self, equation_index: EquationIndex) -> Option<&mut f64> {
        self.b.get_mut((
            equation_index.into_row(self.num_nodes, self.num_voltage_sources)?,
            0,
        ))
    }

    pub fn into_ab(self) -> (DMatrix<f64>, DMatrix<f64>) {
        (self.a, self.b)
    }
}
