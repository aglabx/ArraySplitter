//! Tunable parameters of the `--method autocorr` decomposition pipeline.
//!
//! The four levers exposed here cover the entire autocorr method:
//!   * `excess_floor`        — minimum autocorrelation excess for `find_period`
//!     and `find_period_refined` to accept a candidate period.
//!   * `recursion_termination` — autocorrelation threshold below which the
//!     recursive HOR descent stops (a sub-monomer is treated as a base leaf).
//!   * `multiplet_factor`    — ratio `len / period` at which an anchor-bounded
//!     segment becomes a multiplet and gets split.
//!   * `period_finder`       — which period-detection function runs at the
//!     top level vs inside the recursive descent.
//!
//! Defaults reproduce the historical hard-coded behaviour exactly, so
//! `arraysplitter` with no new flags is byte-identical to the prior binary.

/// Which period-detection function to use.
///
/// `ModeDependent` keeps the historical split:
///   * top level (`decompose_array_autocorr`) uses `find_period` (raw)
///   * recursive descent (`recursive_hor`) uses `find_period_refined`
///
/// `Raw` forces `find_period` at every level; `Refined` forces
/// `find_period_refined` everywhere.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PeriodFinder {
    ModeDependent,
    Raw,
    Refined,
}

impl PeriodFinder {
    /// Parse the CLI string form (`mode-dependent` | `raw` | `refined`).
    pub fn from_cli_str(s: &str) -> Result<Self, String> {
        match s.to_ascii_lowercase().as_str() {
            "mode-dependent" | "mode_dependent" | "modedependent" => Ok(Self::ModeDependent),
            "raw" => Ok(Self::Raw),
            "refined" => Ok(Self::Refined),
            other => Err(format!(
                "invalid --period-finder value: '{}' (expected mode-dependent | raw | refined)",
                other
            )),
        }
    }

    pub fn as_cli_str(&self) -> &'static str {
        match self {
            Self::ModeDependent => "mode-dependent",
            Self::Raw => "raw",
            Self::Refined => "refined",
        }
    }
}

/// Tunable parameters of the autocorr pipeline. See module docs for which
/// downstream call site each field drives.
#[derive(Debug, Clone, Copy)]
pub struct AutocorrParams {
    pub excess_floor: f64,
    pub recursion_termination: f64,
    pub multiplet_factor: f64,
    pub period_finder: PeriodFinder,
}

impl Default for AutocorrParams {
    fn default() -> Self {
        Self {
            excess_floor: 0.05,
            recursion_termination: 0.5,
            multiplet_factor: 1.5,
            period_finder: PeriodFinder::ModeDependent,
        }
    }
}

impl AutocorrParams {
    /// True iff every field equals its historical default. The writer uses
    /// this to keep the `method` column in `summary.tsv` as plain `"autocorr"`
    /// (preserving regression-test golden MD5) instead of the verbose tagged
    /// form when no sweep has actually overridden anything.
    pub fn is_default(&self) -> bool {
        let d = Self::default();
        self.excess_floor == d.excess_floor
            && self.recursion_termination == d.recursion_termination
            && self.multiplet_factor == d.multiplet_factor
            && self.period_finder == d.period_finder
    }

    /// Render the `method` column label used in `summary.tsv`. Returns the
    /// bare base name when at defaults; otherwise appends the param values
    /// in the `base:excess=E;recurse=R;mult=M;finder=F` shape.
    pub fn method_label(&self, base: &str) -> String {
        if self.is_default() {
            base.to_string()
        } else {
            format!(
                "{}:excess={};recurse={};mult={};finder={}",
                base,
                self.excess_floor,
                self.recursion_termination,
                self.multiplet_factor,
                self.period_finder.as_cli_str(),
            )
        }
    }
}
