"""
Trajectory Analysis Module

This module provides comprehensive tools for single-cell trajectory analysis,
including trajectory inference, pseudotime analysis,  and trajectory validation.

Key Components:
- trajectory_inference: Trajectory inference algorithms (AddTrajectory equivalent)
- pseudotime_analysis: Pseudotime analysis and gene dynamics (GetTrajectory equivalent)
- trajectory_validation: Trajectory quality assessment and validation

Main Functions:
- add_trajectory: Infer trajectories using various methods
- get_trajectory_data: Extract trajectory data for analysis
- validate_trajectory_quality: Comprehensive trajectory validation

Example Usage:
```python
import python_scmega as pymega

# Trajectory inference (R-compatible API)
multiome = pymega.add_trajectory(
    multiome,
    trajectory=["0", "1", "2"],  # R: trajectory = c("0", "1", "2")
    group_by="leiden",           # R: group.by = "leiden"
    reduction="pca",
    dims=list(range(30)),
    pre_filter_quantile=0.9,
    post_filter_quantile=0.9,
    dof=250,
    spar=1.0,
    name="Trajectory"
)

# Extract trajectory data (R GetTrajectory equivalent)
trajectory_result = pymega.get_trajectory_data(
    multiome,
    trajectory_name="Trajectory",
    assay="RNA",
    group_every=1,
    log2_norm=True,
    scale_to=10000,
    smooth_window=11
)


# Validate trajectory quality
validation = pymega.validate_pseudotime_ordering(
    multiome,
    trajectory_name="Trajectory",
    known_markers={
        "early": ["CD34", "KIT"],
        "late": ["CD14", "CD16"]
    }
)
```
"""

from .trajectory_inference import (
    add_trajectory,
    add_trajectory_archr
)

from .pseudotime_analysis import (
    get_trajectory_data
)



__all__ = [
    # Trajectory inference
    'add_trajectory',
    'add_trajectory_archr',
    
    # Pseudotime analysis
    'get_trajectory_data',
]
