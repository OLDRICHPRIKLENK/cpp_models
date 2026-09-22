"""
slice_extractor.py
==================
Extraction and slicing utilities for Paul trap VTU datasets.
Supports both unstructured 3D FEM tetrahedral meshes and structured regular grids.
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Tuple, Optional, Union
import numpy as np
import pyvista as pv


@dataclass
class FemSliceData:
    """Container for sliced 2D FEM data."""
    z_val: float
    x: np.ndarray             # 1D array of x coordinates in meters
    y: np.ndarray             # 1D array of y coordinates in meters
    triangles: np.ndarray     # (N, 3) integer indices defining Delaunay/cut triangles
    potential: np.ndarray     # 1D array of scalar potential in Volts
    x_bounds: Tuple[float, float]
    y_bounds: Tuple[float, float]

    @property
    def x_mm(self) -> np.ndarray:
        """x coordinates in millimeters."""
        return self.x * 1e3

    @property
    def y_mm(self) -> np.ndarray:
        """y coordinates in millimeters."""
        return self.y * 1e3


@dataclass
class GridSliceData:
    """Container for sliced 2D regular grid data."""
    z_val: float
    x_1d: np.ndarray          # 1D sorted unique x coordinates in meters
    y_1d: np.ndarray          # 1D sorted unique y coordinates in meters
    X: np.ndarray             # 2D meshgrid of x coordinates in meters
    Y: np.ndarray             # 2D meshgrid of y coordinates in meters
    potential: np.ndarray     # 2D array of scalar potential in Volts
    Ex: np.ndarray            # 2D array of x-component electric field in V/m
    Ey: np.ndarray            # 2D array of y-component electric field in V/m
    Ez: np.ndarray            # 2D array of z-component electric field in V/m
    E_mag: np.ndarray         # 2D array of 3D electric field magnitude in V/m
    E_inplane: np.ndarray     # 2D array of in-plane electric field magnitude in V/m
    x_bounds: Tuple[float, float]
    y_bounds: Tuple[float, float]

    @property
    def x_1d_mm(self) -> np.ndarray:
        """1D x coordinates in millimeters."""
        return self.x_1d * 1e3

    @property
    def y_1d_mm(self) -> np.ndarray:
        """1D y coordinates in millimeters."""
        return self.y_1d * 1e3

    @property
    def X_mm(self) -> np.ndarray:
        """2D x coordinates in millimeters."""
        return self.X * 1e3

    @property
    def Y_mm(self) -> np.ndarray:
        """2D y coordinates in millimeters."""
        return self.Y * 1e3

    @property
    def E_mag_kV_m(self) -> np.ndarray:
        """Electric field magnitude in kV/m."""
        return self.E_mag / 1e3

    @property
    def E_inplane_kV_m(self) -> np.ndarray:
        """In-plane electric field magnitude in kV/m."""
        return self.E_inplane / 1e3


def resolve_file_path(filepath: Union[str, Path]) -> Path:
    """
    Resolve file path, checking current directory, parent directory, and project root.
    """
    p = Path(filepath)
    if p.exists():
        return p.resolve()

    # Try relative to parent directories
    for parent in [Path.cwd(), Path(__file__).resolve().parent, Path(__file__).resolve().parent.parent]:
        candidate = (parent / p).resolve()
        if candidate.exists():
            return candidate

    raise FileNotFoundError(f"VTU file not found: {filepath}. Searched in cwd and project root.")


def load_vtu(filepath: Union[str, Path]) -> pv.UnstructuredGrid:
    """Load a VTU file using PyVista with clear error reporting."""
    resolved = resolve_file_path(filepath)
    try:
        mesh = pv.read(str(resolved))
        return mesh
    except Exception as e:
        raise RuntimeError(f"Failed to read VTU file from {resolved}: {e}") from e


def extract_fem_slice(vtu_source: Union[str, Path, pv.UnstructuredGrid], z_val: float = 0.0) -> FemSliceData:
    """
    Extract a 2D planar slice at z = z_val from an unstructured FEM VTU mesh.

    Parameters:
    -----------
    vtu_source : str, Path, or pv.UnstructuredGrid
        File path to fem_solution.vtu or an already loaded PyVista mesh.
    z_val : float
        z-coordinate of the planar slice in meters (default: 0.0).

    Returns:
    --------
    FemSliceData containing coordinates, triangulation indices, and potential.
    """
    if isinstance(vtu_source, (str, Path)):
        mesh = load_vtu(vtu_source)
    else:
        mesh = vtu_source

    if "Potential" not in mesh.point_data:
        available = list(mesh.point_data.keys())
        raise KeyError(f"'Potential' scalar array not found in FEM mesh. Available arrays: {available}")

    # Slice at normal z
    plane_slice = mesh.slice(normal="z", origin=(0.0, 0.0, z_val))
    if plane_slice.n_points == 0:
        z_bounds = mesh.bounds[4:6]
        raise ValueError(
            f"Planar cut at z={z_val:.4e} m returned 0 points. Mesh z bounds are [{z_bounds[0]:.4e}, {z_bounds[1]:.4e}] m."
        )

    # Triangulate polygon faces for clean Matplotlib Triangulation
    tri_slice = plane_slice.triangulate()
    faces = tri_slice.faces.reshape(-1, 4)
    triangles = faces[:, 1:4]

    pts = tri_slice.points
    x = pts[:, 0]
    y = pts[:, 1]
    potential = tri_slice.point_data["Potential"]

    return FemSliceData(
        z_val=z_val,
        x=x,
        y=y,
        triangles=triangles,
        potential=potential,
        x_bounds=(float(x.min()), float(x.max())),
        y_bounds=(float(y.min()), float(y.max()))
    )


def extract_grid_slice(vtu_source: Union[str, Path, pv.UnstructuredGrid], z_val: float = 0.0) -> GridSliceData:
    """
    Extract a 2D planar slice at z = z_val from a regular/interpolated grid VTU.

    Parameters:
    -----------
    vtu_source : str, Path, or pv.UnstructuredGrid
        File path to interpolated_grid.vtu or an already loaded PyVista mesh.
    z_val : float
        z-coordinate of the planar slice in meters (default: 0.0).

    Returns:
    --------
    GridSliceData containing 2D structured arrays of potential and electric field.
    """
    if isinstance(vtu_source, (str, Path)):
        mesh = load_vtu(vtu_source)
    else:
        mesh = vtu_source

    if "Potential" not in mesh.point_data:
        raise KeyError(f"'Potential' scalar array not found in grid mesh. Available: {list(mesh.point_data.keys())}")
    if "ElectricField" not in mesh.point_data:
        raise KeyError(f"'ElectricField' vector array not found in grid mesh. Available: {list(mesh.point_data.keys())}")

    plane_slice = mesh.slice(normal="z", origin=(0.0, 0.0, z_val))
    if plane_slice.n_points == 0:
        z_bounds = mesh.bounds[4:6]
        raise ValueError(
            f"Planar cut at z={z_val:.4e} m returned 0 points. Mesh z bounds are [{z_bounds[0]:.4e}, {z_bounds[1]:.4e}] m."
        )

    pts = plane_slice.points
    pot = plane_slice.point_data["Potential"]
    ef = plane_slice.point_data["ElectricField"]

    # Extract unique coordinates along x and y with precision rounding (nanometer tolerance)
    xs_unique = np.sort(np.unique(np.round(pts[:, 0], 9)))
    ys_unique = np.sort(np.unique(np.round(pts[:, 1], 9)))
    nx = len(xs_unique)
    ny = len(ys_unique)

    # Map scattered slice points into 2D structured grid
    x_idx = np.searchsorted(xs_unique, np.round(pts[:, 0], 9))
    y_idx = np.searchsorted(ys_unique, np.round(pts[:, 1], 9))

    V_grid = np.zeros((ny, nx), dtype=float)
    Ex_grid = np.zeros((ny, nx), dtype=float)
    Ey_grid = np.zeros((ny, nx), dtype=float)
    Ez_grid = np.zeros((ny, nx), dtype=float)

    V_grid[y_idx, x_idx] = pot
    Ex_grid[y_idx, x_idx] = ef[:, 0]
    Ey_grid[y_idx, x_idx] = ef[:, 1]
    Ez_grid[y_idx, x_idx] = ef[:, 2]

    X, Y = np.meshgrid(xs_unique, ys_unique)
    E_mag = np.sqrt(Ex_grid**2 + Ey_grid**2 + Ez_grid**2)
    E_inplane = np.hypot(Ex_grid, Ey_grid)

    return GridSliceData(
        z_val=z_val,
        x_1d=xs_unique,
        y_1d=ys_unique,
        X=X,
        Y=Y,
        potential=V_grid,
        Ex=Ex_grid,
        Ey=Ey_grid,
        Ez=Ez_grid,
        E_mag=E_mag,
        E_inplane=E_inplane,
        x_bounds=(float(xs_unique.min()), float(xs_unique.max())),
        y_bounds=(float(ys_unique.min()), float(ys_unique.max()))
    )
