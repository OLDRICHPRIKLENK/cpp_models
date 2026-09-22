"""
downstream_analysis package
===========================
Python downstream processing and 2D planar slicing for 3D Paul Trap VTU simulations.
"""

from .slice_extractor import (
    FemSliceData,
    GridSliceData,
    load_vtu,
    extract_fem_slice,
    extract_grid_slice,
)

from .visualizer import (
    plot_fem_potential,
    plot_grid_potential,
    plot_grid_electric_field,
    plot_combined_dashboard,
    plot_1d_cuts,
    save_figure,
)

__all__ = [
    "FemSliceData",
    "GridSliceData",
    "load_vtu",
    "extract_fem_slice",
    "extract_grid_slice",
    "plot_fem_potential",
    "plot_grid_potential",
    "plot_grid_electric_field",
    "plot_combined_dashboard",
    "plot_1d_cuts",
    "save_figure",
]
