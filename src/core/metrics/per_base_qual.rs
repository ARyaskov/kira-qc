use crate::core::model::BaseGroup;

#[derive(Clone, Debug)]
pub struct PerBaseQualRow {
    pub group: BaseGroup,
    pub mean: f64,
    pub median: u8,
    pub lower_quartile: u8,
    pub upper_quartile: u8,
    pub p10: u8,
    pub p90: u8,
}
