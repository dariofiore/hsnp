//use crate::String;

/// The error type for `HSNP`.
#[derive(Debug)]
pub enum Error {
    /// The constraints < num inputs>
    TooFewConstraints,
    /// The index is too large for the universal public parameters.
    IndexTooLarge,

}

impl core::fmt::Display for Error {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Error::TooFewConstraints => write!(
                f,
                "Too few constraints",
            ),
            Error::IndexTooLarge => write!(
                f,
                "Index too large",
            ),
        }
    }
}

impl ark_std::error::Error for Error {}
