use nalgebra::DMatrix;

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ViewEquationIndex {
    NodalEquation(usize),
    SpecificEquation(usize),
}

impl ViewEquationIndex {
    fn into_global_index(
        self,
        num_nodes: usize,
        num_variables: usize,
        variables_start: usize,
    ) -> Option<usize> {
        match self {
            Self::NodalEquation(idx) => {
                if idx == 0 {
                    return None;
                }

                if idx > num_nodes {
                    return None;
                }

                Some(idx - 1)
            }
            Self::SpecificEquation(idx) => {
                if idx >= num_variables {
                    return None;
                }

                Some(variables_start + idx)
            }
        }
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ViewVariableIndex {
    NodeVoltage(usize),
    SpecificVariable(usize),
}

impl ViewVariableIndex {
    fn into_global_index(
        self,
        num_nodes: usize,
        num_variables: usize,
        variables_start: usize,
    ) -> Option<usize> {
        match self {
            Self::NodeVoltage(idx) => {
                if idx == 0 {
                    return None;
                }

                if idx > num_nodes {
                    return None;
                }

                Some(idx - 1)
            }
            Self::SpecificVariable(idx) => {
                if idx >= num_variables {
                    return None;
                }

                Some(variables_start + idx)
            }
        }
    }
}

pub struct ABMatrixView<'a> {
    a: &'a mut DMatrix<f64>,
    b: &'a mut DMatrix<f64>,
    num_nodes: usize,
    num_variables: usize,
    variables_start: usize,
}

impl<'a> ABMatrixView<'a> {
    pub fn new(
        a: &'a mut DMatrix<f64>,
        b: &'a mut DMatrix<f64>,
        num_nodes: usize,
        num_variables: usize,
        variables_start: usize,
    ) -> Self {
        Self {
            a,
            b,
            num_nodes,
            num_variables,
            variables_start,
        }
    }

    fn get_coefficient_mut(
        &mut self,
        equation: ViewEquationIndex,
        variable: ViewVariableIndex,
    ) -> Option<&mut f64> {
        self.a.get_mut((
            equation.into_global_index(self.num_nodes, self.num_variables, self.variables_start)?,
            variable.into_global_index(self.num_nodes, self.num_variables, self.variables_start)?,
        ))
    }

    pub fn coefficient_add(
        &mut self,
        equation: ViewEquationIndex,
        variable: ViewVariableIndex,
        value: f64,
    ) {
        if let Some(a) = self.get_coefficient_mut(equation, variable) {
            *a += value;
        }
    }

    fn get_result_mut(&mut self, equation: ViewEquationIndex) -> Option<&mut f64> {
        self.b.get_mut((
            equation.into_global_index(self.num_nodes, self.num_variables, self.variables_start)?,
            0,
        ))
    }

    pub fn result_add(&mut self, equation: ViewEquationIndex, value: f64) {
        if let Some(a) = self.get_result_mut(equation) {
            *a += value;
        }
    }
}

pub struct XMatrixView<'a> {
    x: &'a DMatrix<f64>,
    num_nodes: usize,
    num_variables: usize,
    variables_start: usize,
}

impl<'a> XMatrixView<'a> {
    pub fn new(
        x: &'a DMatrix<f64>,
        num_nodes: usize,
        num_variables: usize,
        variables_start: usize,
    ) -> Self {
        Self {
            x,
            num_nodes,
            num_variables,
            variables_start,
        }
    }

    pub fn get_variable(&self, variable: ViewVariableIndex) -> Option<f64> {
        match variable {
            ViewVariableIndex::NodeVoltage(0) => Some(0.0),
            _ => self
                .x
                .get((
                    variable.into_global_index(
                        self.num_nodes,
                        self.num_variables,
                        self.variables_start,
                    )?,
                    0,
                ))
                .cloned(),
        }
    }
}
