use crate::core::model::BaseGroup;

#[derive(Clone, Debug)]
pub struct PerBaseNRow {
    pub group: BaseGroup,
    pub n_percent: f64,
}
