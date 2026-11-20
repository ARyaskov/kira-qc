use crate::core::model::BaseGroup;

#[derive(Clone, Debug)]
pub struct PerBaseContentRow {
    pub group: BaseGroup,
    pub g: f64,
    pub a: f64,
    pub t: f64,
    pub c: f64,
}
