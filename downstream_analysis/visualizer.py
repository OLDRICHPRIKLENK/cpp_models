"""
visualizer.py
=============
Publication-quality Matplotlib visualization routines for Paul trap 2D slices.
Visualizes Electric Potential and Electric Field distributions.
"""

from pathlib import Path
from typing import Optional, Union, Tuple
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.tri as mtri
import numpy as np

from .slice_extractor import FemSliceData, GridSliceData


def setup_matplotlib_style():
    """Apply clean, modern publication styling to Matplotlib."""
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.size": 11,
        "axes.labelsize": 12,
        "axes.titlesize": 13,
        "xtick.labelsize": 10,
        "ytick.labelsize": 10,
        "legend.fontsize": 10,
        "figure.titlesize": 15,
        "axes.linewidth": 1.0,
        "grid.linewidth": 0.5,
        "grid.alpha": 0.35,
    })


def plot_fem_potential(
    fem_data: FemSliceData,
    ax: Optional[plt.Axes] = None,
    levels: int = 50,
    cmap: str = "coolwarm",
    show_colorbar: bool = True,
    show_contours: bool = True,
    show_mesh: bool = False,
    electrode_radius_mm: Optional[float] = 5.0,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot 2D Electric Potential from unstructured FEM slice.

    Parameters:
    -----------
    fem_data : FemSliceData
        FEM slice data container.
    ax : Optional[plt.Axes]
        Matplotlib Axes to plot onto. If None, a new figure is created.
    levels : int
        Number of contour color levels.
    cmap : str
        Colormap name.
    show_colorbar : bool
        Whether to add a colorbar.
    show_contours : bool
        Whether to overlay labeled equipotential curves.
    show_mesh : bool
        Whether to overlay mesh triangle edges.
    electrode_radius_mm : Optional[float]
        Radius in mm to draw reference electrode boundary circle.

    Returns:
    --------
    (fig, ax) tuple.
    """
    setup_matplotlib_style()
    if ax is None:
        fig, ax = plt.subplots(figsize=(6.5, 5.5))
    else:
        fig = ax.figure

    triang = mtri.Triangulation(fem_data.x_mm, fem_data.y_mm, triangles=fem_data.triangles)

    # Filled potential contour
    cf = ax.tricontourf(triang, fem_data.potential, levels=levels, cmap=cmap)

    # Equipotential lines with labels
    if show_contours:
        c_lines = ax.tricontour(triang, fem_data.potential, levels=14, colors="k", linewidths=0.6, alpha=0.55)
        ax.clabel(c_lines, inline=True, fontsize=8, fmt="%.0f V")

    # Optional mesh wireframe
    if show_mesh:
        ax.triplot(triang, color="gray", linewidth=0.25, alpha=0.4)

    # Reference electrode cylinder boundary
    if electrode_radius_mm is not None:
        theta = np.linspace(0, 2 * np.pi, 200)
        ax.plot(
            electrode_radius_mm * np.cos(theta),
            electrode_radius_mm * np.sin(theta),
            "k--",
            linewidth=1.0,
            alpha=0.6,
            label=f"Electrode Surface (r₀ = {electrode_radius_mm:.1f} mm)"
        )

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_title(f"FEM Electric Potential (z = {fem_data.z_val*1e3:.2f} mm)")

    if show_colorbar:
        cbar = fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("Potential $\\Phi$ [V]")

    return fig, ax


def plot_grid_potential(
    grid_data: GridSliceData,
    ax: Optional[plt.Axes] = None,
    levels: int = 50,
    cmap: str = "coolwarm",
    show_colorbar: bool = True,
    show_contours: bool = True,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot 2D Electric Potential from interpolated Cartesian grid slice.
    """
    setup_matplotlib_style()
    if ax is None:
        fig, ax = plt.subplots(figsize=(6.5, 5.5))
    else:
        fig = ax.figure

    cf = ax.contourf(grid_data.x_1d_mm, grid_data.y_1d_mm, grid_data.potential, levels=levels, cmap=cmap)

    if show_contours:
        c_lines = ax.contour(
            grid_data.x_1d_mm,
            grid_data.y_1d_mm,
            grid_data.potential,
            levels=14,
            colors="k",
            linewidths=0.6,
            alpha=0.55
        )
        ax.clabel(c_lines, inline=True, fontsize=8, fmt="%.1f V")

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_title(f"Interpolated Potential (z = {grid_data.z_val*1e3:.2f} mm)")

    if show_colorbar:
        cbar = fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("Potential $\\Phi$ [V]")

    return fig, ax


def plot_grid_electric_field(
    grid_data: GridSliceData,
    ax: Optional[plt.Axes] = None,
    style: str = "streamlines",
    levels: int = 50,
    cmap: str = "viridis",
    show_colorbar: bool = True,
    overlay_equipotentials: bool = False,
    quiver_subsample: int = 2,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot 2D Electric Field magnitude and vectors from interpolated grid slice.

    Parameters:
    -----------
    grid_data : GridSliceData
        Interpolated grid slice data container.
    ax : Optional[plt.Axes]
        Target axes.
    style : str
        'streamlines', 'quiver', or 'both'.
    levels : int
        Number of contour levels for magnitude heatmap.
    cmap : str
        Colormap name.
    show_colorbar : bool
        Whether to add colorbar.
    overlay_equipotentials : bool
        Whether to overlay equipotential contours to highlight E perpendicular to V.
    quiver_subsample : int
        Subsampling step for quiver arrows to avoid clutter.
    """
    setup_matplotlib_style()
    if ax is None:
        fig, ax = plt.subplots(figsize=(6.5, 5.5))
    else:
        fig = ax.figure

    # In-plane magnitude in kV/m
    cf = ax.contourf(grid_data.x_1d_mm, grid_data.y_1d_mm, grid_data.E_inplane_kV_m, levels=levels, cmap=cmap)

    # Optional equipotential lines (dashed) to demonstrate orthogonality
    if overlay_equipotentials:
        ax.contour(
            grid_data.x_1d_mm,
            grid_data.y_1d_mm,
            grid_data.potential,
            levels=10,
            colors="cyan",
            linestyles="dashed",
            linewidths=0.7,
            alpha=0.5
        )

    # Electric field vectors
    if style in ("streamlines", "both"):
        ax.streamplot(
            grid_data.x_1d_mm,
            grid_data.y_1d_mm,
            grid_data.Ex,
            grid_data.Ey,
            color="white",
            linewidth=0.9,
            density=1.1,
            arrowsize=1.0
        )

    if style in ("quiver", "both"):
        sub = quiver_subsample
        X_sub = grid_data.X_mm[::sub, ::sub]
        Y_sub = grid_data.Y_mm[::sub, ::sub]
        Ex_sub = grid_data.Ex[::sub, ::sub]
        Ey_sub = grid_data.Ey[::sub, ::sub]
        arrow_color = "red" if style == "both" else "white"
        ax.quiver(
            X_sub,
            Y_sub,
            Ex_sub,
            Ey_sub,
            color=arrow_color,
            scale=None,
            pivot="middle",
            alpha=0.85
        )

    ax.set_aspect("equal", adjustable="box")
    ax.set_xlabel("x [mm]")
    ax.set_ylabel("y [mm]")
    ax.set_title(f"Electric Field $\\vec{{E}} = -\\nabla V$ (z = {grid_data.z_val*1e3:.2f} mm)")

    if show_colorbar:
        cbar = fig.colorbar(cf, ax=ax, fraction=0.046, pad=0.04)
        cbar.set_label("Electric Field $|\\vec{E}_{xy}|$ [kV/m]")

    return fig, ax


def plot_1d_cuts(
    grid_data: GridSliceData,
    fem_data: Optional[FemSliceData] = None,
    ax: Optional[plt.Axes] = None,
) -> Tuple[plt.Figure, plt.Axes]:
    """
    Plot 1D cross-sectional potential profiles along x-axis (y=0) and y-axis (x=0).
    Demonstrates the saddle potential and harmonic trap confinement.
    """
    setup_matplotlib_style()
    if ax is None:
        fig, ax = plt.subplots(figsize=(6.5, 5.5))
    else:
        fig = ax.figure

    # Find y closest to 0
    iy_zero = int(np.argmin(np.abs(grid_data.y_1d)))
    ix_zero = int(np.argmin(np.abs(grid_data.x_1d)))

    x_line = grid_data.x_1d_mm
    y_line = grid_data.y_1d_mm
    Vx_line = grid_data.potential[iy_zero, :]
    Vy_line = grid_data.potential[:, ix_zero]

    ax.plot(x_line, Vx_line, "b-", linewidth=2.0, label=f"V(x, y=0) [Confining Axis]")
    ax.plot(y_line, Vy_line, "r--", linewidth=2.0, label=f"V(x=0, y) [Anti-confining Axis]")

    ax.axhline(0, color="gray", linestyle=":", alpha=0.5)
    ax.axvline(0, color="gray", linestyle=":", alpha=0.5)
    ax.set_xlabel("Coordinate [mm]")
    ax.set_ylabel("Potential [V]")
    ax.set_title(f"1D Central Potential Well Profiles (z = {grid_data.z_val*1e3:.2f} mm)")
    ax.legend(frameon=True, loc="upper center")
    ax.grid(True, linestyle="--", alpha=0.5)

    return fig, ax


def plot_combined_dashboard(
    fem_data: FemSliceData,
    grid_data: GridSliceData,
    efield_style: str = "streamlines",
    save_path: Optional[Union[str, Path]] = None,
    dpi: int = 300,
) -> plt.Figure:
    """
    Create a 2x2 comprehensive downstream analysis dashboard:
    1. FEM Potential (Entire physical domain)
    2. Interpolated Grid Potential (Central trap zone)
    3. Electric Field Magnitude & Field Lines (Streamlines / Quiver)
    4. 1D Line Cut Profiles (Potential curvature along x and y)
    """
    
    setup_matplotlib_style()
    fig, axs = plt.subplots(2, 2, figsize=(14, 12))

    # Panel (0, 0): FEM Potential
    plot_fem_potential(fem_data, ax=axs[0, 0], show_colorbar=True)
    axs[0, 0].set_title(f"(a) Physics Engine Potential: Full Domain (z = {fem_data.z_val*1e3:.2f} mm)", fontweight="bold")

    # Panel (0, 1): Grid Potential
    plot_grid_potential(grid_data, ax=axs[0, 1], show_colorbar=True)
    axs[0, 1].set_title(f"(b) Interpolated Potential: Central Zone (z = {grid_data.z_val*1e3:.2f} mm)", fontweight="bold")

    # Panel (1, 0): Grid Electric Field
    plot_grid_electric_field(grid_data, ax=axs[1, 0], style=efield_style, show_colorbar=True, overlay_equipotentials=True)
    axs[1, 0].set_title(f"(c) Electric Field: Streamlines & |E| (z = {grid_data.z_val*1e3:.2f} mm)", fontweight="bold")

    # Panel (1, 1): 1D Line Cuts
    plot_1d_cuts(grid_data, fem_data=fem_data, ax=axs[1, 1])
    axs[1, 1].set_title(f"(d) 1D Central Potential Well Profiles", fontweight="bold")

    plt.suptitle(
        f"3D Paul Trap Quadrupole Slice Analysis (z = {grid_data.z_val*1e3:.2f} mm)",
        fontsize=16,
        fontweight="bold",
        y=0.98
    )
    plt.tight_layout(rect=[0, 0, 1, 0.96])

    if save_path is not None:
        save_figure(fig, save_path, dpi=dpi)

    return fig


def save_figure(fig: plt.Figure, filepath: Union[str, Path], dpi: int = 300):
    """Save a Matplotlib figure ensuring parent directories exist."""
    p = Path(filepath)
    p.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(p), dpi=dpi, bbox_inches="tight")
    print(f"[Visualizer] Saved figure to: {p.resolve()}")
