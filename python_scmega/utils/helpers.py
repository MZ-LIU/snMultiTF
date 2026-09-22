"""
Helper Functions Module

Core utility functions adapted from the original R scMEGA package.
These functions are used throughout the analysis pipeline, particularly
in trajectory analysis and data processing.

Key functions:
- get_quantiles: Quantile calculation (getQuantiles equivalent)
- center_rolling_mean: Centered rolling mean (centerRollMean equivalent)

References:
- Original R code: R/utils.R
- ArchR package hidden utilities: https://github.com/GreenleafLab/ArchR/blob/master/R/HiddenUtils.R
"""

import numpy as np
from typing import Union, Optional
from scipy import stats


def get_quantiles(v: Union[np.ndarray, list], 
                 length: Optional[int] = None) -> np.ndarray:
    """
    Calculate quantiles for a vector - R-compatible implementation.
    
    This function replicates the getQuantiles function from the original R code:
    ```R
    getQuantiles <- function(v = NULL, len = length(v)){
      if(length(v) < len){
        v2 <- rep(0, len)
        v2[seq_along(v)] <- v
      }else{
        v2 <- v
      }
      p <- trunc(rank(v2))/length(v2)  # R's rank() uses method="average" by default
      if(length(v) < len){
        p <- p[seq_along(v)]
      }
      return(p)
    }
    ```
    
    This function is crucial for trajectory analysis, particularly in the
    AddTrajectory function for pseudotime calculation.
    
    **R compatibility note**: R's rank() function uses method="average" by default,
    which assigns the average rank to tied values. This differs from method="min"
    or method="ordinal".
    
    Args:
        v: Input vector/array
        length: Target length (default: length of v)
        
    Returns:
        Array of quantile values between 0 and 1
        
    Example:
        ```python
        import numpy as np
        values = np.array([1, 5, 3, 8, 2])
        quantiles = get_quantiles(values)
        # Returns quantile ranks normalized to [0, 1]
        ```
    """
    v = np.asarray(v)
    
    if length is None:
        length = len(v)
    
    # Handle case where input is shorter than target length
    if len(v) < length:
        v2 = np.zeros(length)
        v2[:len(v)] = v
    else:
        v2 = v.copy()
    
    # Calculate truncated ranks and normalize
    # IMPORTANT: Use 'average' method to match R's rank() default behavior
    # This is critical for handling tied values correctly
    ranks = stats.rankdata(v2, method='average')  # R-compatible: average rank for ties
    p = np.floor(ranks) / len(v2)
    
    # Return only original length if input was padded
    if len(v) < length:
        p = p[:len(v)]
    
    return p


def center_rolling_mean(v: Union[np.ndarray, list], 
                       k: int) -> np.ndarray:
    """
    Calculate centered rolling mean with edge handling - OPTIMIZED R-compatible version.
    
    This function replicates the centerRollMean function from the original R code:
    ```R
    centerRollMean <- function(v = NULL, k = NULL){
      o1 <- data.table::frollmean(v, k, align = "right", na.rm = FALSE)
      if(k%%2==0){
        o2 <- c(rep(o1[k], floor(k/2)-1), o1[-seq_len(k-1)], rep(o1[length(o1)], floor(k/2)))
      }else if(k%%2==1){
        o2 <- c(rep(o1[k], floor(k/2)), o1[-seq_len(k-1)], rep(o1[length(o1)], floor(k/2)))
      }else{
        stop("Error!")
      }
      o2
    }
    ```
    
    This function is used in trajectory analysis for smoothing data along
    pseudotime trajectories (GetTrajectory function).
    
    **Performance optimization**: Uses pandas.rolling() instead of Python loops
    for 10-100x speedup on large arrays, matching R's data.table::frollmean performance.
    
    Args:
        v: Input vector/array
        k: Window size for rolling mean
        
    Returns:
        Array with centered rolling mean values
        
    Raises:
        ValueError: If k is not a positive integer
        
    Example:
        ```python
        import numpy as np
        data = np.array([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
        smoothed = center_rolling_mean(data, k=3)
        # Returns smoothed trajectory data
        ```
    """
    import pandas as pd
    
    v = np.asarray(v)
    
    if not isinstance(k, int) or k <= 0:
        raise ValueError("k must be a positive integer")
    
    if k > len(v):
        raise ValueError("Window size k cannot be larger than array length")
    
    # OPTIMIZED: Use pandas rolling mean (vectorized, ~10-100x faster than Python loop)
    # This replicates R's data.table::frollmean behavior with align="right"
    rolling_mean = pd.Series(v).rolling(window=k, min_periods=k).mean().values
    
    # Handle edge cases based on window size parity (following R logic exactly)
    if k % 2 == 0:
        # Even window size
        left_pad = k // 2 - 1
        right_pad = k // 2
    elif k % 2 == 1:
        # Odd window size  
        left_pad = k // 2
        right_pad = k // 2
    else:
        raise ValueError("Error!")  # This matches the R code error condition
    
    # Create output array with edge handling (following R logic)
    result = np.full(len(v), np.nan)
    
    # Fill the main part (excluding edges) - equivalent to o1[-seq_len(k-1)]
    valid_start = k - 1
    result[left_pad:len(v)-right_pad] = rolling_mean[valid_start:len(v)]
    
    # Fill left edge with first valid value - equivalent to rep(o1[k], ...)
    if not np.isnan(rolling_mean[k-1]):
        result[:left_pad] = rolling_mean[k-1]
    
    # Fill right edge with last valid value - equivalent to rep(o1[length(o1)], ...)
    if not np.isnan(rolling_mean[-1]):
        result[len(v)-right_pad:] = rolling_mean[-1]
    
    return result


def safe_divide(numerator: Union[np.ndarray, float], 
                denominator: Union[np.ndarray, float],
                fill_value: float = 0.0) -> Union[np.ndarray, float]:
    """
    Perform safe division with handling of zero denominators.
    
    Args:
        numerator: Numerator values
        denominator: Denominator values  
        fill_value: Value to use when denominator is zero
        
    Returns:
        Division result with safe handling of zero denominators
    """
    if isinstance(denominator, (int, float)):
        if denominator == 0:
            return fill_value
        return numerator / denominator
    
    # Handle array case
    numerator = np.asarray(numerator)
    denominator = np.asarray(denominator)
    
    result = np.full_like(numerator, fill_value, dtype=float)
    mask = denominator != 0
    result[mask] = numerator[mask] / denominator[mask]
    
    return result


def match_cells(cells1: list, cells2: list) -> tuple:
    """
    Find matching cells between two cell lists.
    
    Equivalent to R's intersect() function for finding shared cells
    between RNA and ATAC data.
    
    Args:
        cells1: First list of cell barcodes
        cells2: Second list of cell barcodes
        
    Returns:
        Tuple of (shared_cells, indices_in_cells1, indices_in_cells2)
    """
    cells1_set = set(cells1)
    cells2_set = set(cells2)
    
    shared_cells = sorted(list(cells1_set.intersection(cells2_set)))
    
    # Get indices for shared cells in original lists
    indices1 = [cells1.index(cell) for cell in shared_cells if cell in cells1]
    indices2 = [cells2.index(cell) for cell in shared_cells if cell in cells2]
    
    return shared_cells, indices1, indices2


def filter_chromosomes(peak_names: list, 
                      chr_pattern: str = "chr",
                      exclude_patterns: Optional[list] = None) -> np.ndarray:
    """
    Filter peaks by chromosome names.
    
    Equivalent to R code: grep("chr", rownames(atac_counts))
    
    Args:
        peak_names: List of peak names (e.g., "chr1:1000-2000")
        chr_pattern: Pattern to match for valid chromosomes
        exclude_patterns: Patterns to exclude (e.g., ["chrM", "chrY"])
        
    Returns:
        Boolean mask for valid peaks
    """
    if exclude_patterns is None:
        exclude_patterns = []
    
    mask = np.array([
        peak.startswith(chr_pattern) and 
        not any(exclude in peak for exclude in exclude_patterns)
        for peak in peak_names
    ])
    
    return mask


def parse_peak_name(peak_name: str) -> tuple:
    """
    Parse peak name into chromosome, start, and end coordinates.
    
    Args:
        peak_name: Peak name in format "chr1:1000-2000"
        
    Returns:
        Tuple of (chromosome, start, end)
        
    Raises:
        ValueError: If peak name format is invalid
    """
    try:
        chr_part, coord_part = peak_name.split(":")
        start_str, end_str = coord_part.split("-")
        return chr_part, int(start_str), int(end_str)
    except (ValueError, IndexError):
        raise ValueError(f"Invalid peak name format: {peak_name}. Expected format: chr:start-end")


def format_peak_name(chromosome: str, start: int, end: int) -> str:
    """
    Format chromosome coordinates into peak name.
    
    Args:
        chromosome: Chromosome name (e.g., "chr1")
        start: Start coordinate
        end: End coordinate
        
    Returns:
        Peak name in format "chr1:1000-2000"
    """
    return f"{chromosome}:{start}-{end}"



#for saving dictionaries
def save_stuff(stuff,path):
    u"""for saving dictionaries, but probably works with lists and other pickleable objects"""
    import pickle
    with open(path+u'.pickle', u'wb') as handle:
        pickle.dump(stuff, handle, protocol=pickle.HIGHEST_PROTOCOL)


def load_stuff(path,encoding='ASCII'):
    """for loading object saved using 'save_stuff'.
    I had to use encoding='bytes' to load in python3 
    certain data pickled in python2."""
    import pickle
    with open(path, u'rb') as handle:
        return pickle.load(handle,encoding=encoding)
    
def now():
    """Generate a filename-friendly string for the current date and time"""
    return datetime.datetime.now().strftime('%y%m%d_%Hh%M')

def oset(a_list):
    """Given a list, return an ordered set (remove duplication but maintain order)"""
    seen = set()
    seen_add = seen.add
    return [x for x in a_list if not (x in seen or seen_add(x))]
