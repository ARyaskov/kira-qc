//! Optional `KIRA_STATS=1` timing output, shared by the CLI and the engine.

use std::time::{Duration, Instant};

#[derive(Clone, Copy, Debug)]
pub struct StatsLog {
    enabled: bool,
}

impl StatsLog {
    pub fn from_env() -> Self {
        Self {
            enabled: matches!(std::env::var("KIRA_STATS").as_deref(), Ok("1")),
        }
    }

    pub fn enabled(&self) -> bool {
        self.enabled
    }

    /// Report the time elapsed since `start` under `name`.
    pub fn stage(&self, name: &str, start: Instant) {
        if self.enabled {
            eprintln!(
                "KIRA_STATS stage={} time={}",
                name,
                fmt_dur(start.elapsed())
            );
        }
    }

    pub fn line(&self, args: std::fmt::Arguments<'_>) {
        if self.enabled {
            eprintln!("KIRA_STATS {}", args);
        }
    }
}

pub fn fmt_dur(d: Duration) -> String {
    if d.as_secs_f64() < 1.0 {
        format!("{}ms", d.as_millis())
    } else {
        format!("{:.3}s", d.as_secs_f64())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn durations_switch_units_at_one_second() {
        assert_eq!(fmt_dur(Duration::from_millis(250)), "250ms");
        assert_eq!(fmt_dur(Duration::from_millis(1500)), "1.500s");
    }
}
