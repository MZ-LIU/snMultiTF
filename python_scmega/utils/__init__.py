"""
Utilities Module

Core utility functions used throughout the PyMEGA package.

Key components:
- helpers.py: Core utility functions (getQuantiles, centerRollMean equivalent)
  - get_quantiles: Quantile calculation for trajectory analysis
  - center_rolling_mean: Centered rolling mean for data smoothing
  - safe_divide: Safe division with zero handling
  - match_cells: Find matching cells between datasets
  - filter_chromosomes: Filter peaks by chromosome
  - parse_peak_name: Parse peak coordinates
  - format_peak_name: Format peak names

References:
- Original R code: R/utils.R and helper functions
- ArchR package utilities
"""

from .helpers import (
    get_quantiles,
    center_rolling_mean,
    safe_divide,
    match_cells,
    filter_chromosomes,
    parse_peak_name,
    format_peak_name
)

__all__ = [
    "get_quantiles",
    "center_rolling_mean",
    "safe_divide",
    "match_cells",
    "filter_chromosomes",
    "parse_peak_name",
    "format_peak_name"
]
