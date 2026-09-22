#!/usr/bin/env python3
"""
slice_plotter.py
================
CLI tool for extracting and displaying 2D mid-plane (z = 0) slices of 
FEM potential and interpolated grid electric field from Paul trap VTU simulations.

Usage:
------
# Basic run (extracts at z=0 and saves figures to downstream_analysis/plots/):
python slice_plotter.py

# Interactive display:
python slice_plotter.py --show

# Custom slice position and output directory:
python slice_plotter.py --z-val 0.0 --output-dir ./plots --format pdf
"""

import argparse
import sys
from pathlib import Path

# Add project directory / parent directory to path so imports work seamlessly
script_dir = Path(__file__).resolve().parent
if str(script_dir.parent) not in sys.path:
    sys.path.insert(0, str(script_dir.parent))
if str(script_dir) not in sys.path:
    sys.path.insert(0, str(script_dir))

import matplotlib

from downstream_analysis.slice_extractor import (
    extract_fem_slice,
    extract_grid_slice,
    resolve_file_path,
)
from downstream_analysis.visualizer import (
    plot_fem_potential,
    plot_grid_potential,
    plot_grid_electric_field,
    plot_combined_dashboard,
    plot_1d_cuts,
    save_figure,
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Extract and plot 2D planar slices (z=0) of FEM potential and interpolated electric field."
    )
    parser.add_argument(
        "--fem-vtu",
        type=str,
        default="output/vtk/fem_solution.vtu",
        help="Path to the FEM solution VTU file (default: output/vtk/fem_solution.vtu)"
    )
    parser.add_argument(
        "--grid-vtu",
        type=str,
        default="output/vtk/interpolated_grid.vtu",
        help="Path to the interpolated grid VTU file (default: output/vtk/interpolated_grid.vtu)"
    )
    parser.add_argument(
        "--z-val",
        type=float,
        default=0.0,
        help="Z coordinate for planar slice in meters (default: 0.0)"
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default=str(script_dir / "plots"),
        help="Directory to save generated plot files (default: downstream_analysis/plots)"
    )
    parser.add_argument(
        "--mode",
        type=str,
        choices=["all", "dashboard", "fem", "grid-potential", "grid-efield", "cuts"],
        default="all",
        help="Which plots to generate (default: all)"
    )
    parser.add_argument(
        "--efield-style",
        type=str,
        choices=["streamlines", "quiver", "both"],
        default="streamlines",
        help="Vector visualization style for electric field (default: streamlines)"
    )
    parser.add_argument(
        "--cmap-potential",
        type=str,
        default="coolwarm",
        help="Colormap for electric potential (default: coolwarm)"
    )
    parser.add_argument(
        "--cmap-efield",
        type=str,
        default="viridis",
        help="Colormap for electric field magnitude (default: viridis)"
    )
    parser.add_argument(
        "--format",
        type=str,
        choices=["png", "pdf", "svg"],
        default="png",
        help="File format for exported plots (default: png)"
    )
    parser.add_argument(
        "--dpi",
        type=int,
        default=300,
        help="Resolution for exported raster plots (default: 300)"
    )
    parser.add_argument(
        "--show-mesh",
        action="store_true",
        help="Overlay FEM mesh wireframe edges on FEM slice plot"
    )
    parser.add_argument(
        "--show",
        action="store_true",
        help="Display figures interactively using Matplotlib GUI window"
    )
    parser.add_argument(
        "--no-save",
        action="store_true",
        help="Do not save plots to disk (useful if only using --show)"
    )
    return parser.parse_args()


def main():
    args = parse_args()

    # Set non-interactive backend if not displaying GUI
    if not args.show:
        matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out_dir = Path(args.output_dir)
    fmt = args.format
    z_str = f"z{int(round(args.z_val*1e3))}mm" if args.z_val != 0.0 else "z0"

    print("==================================================================")
    print("  Paul Trap Downstream Analysis: 2D Planar Field Slicing Tool     ")
    print("==================================================================")
    print(f"Slice plane:          z = {args.z_val:.4e} m ({args.z_val*1e3:.2f} mm)")
    print(f"FEM VTU path:         {args.fem_vtu}")
    print(f"Grid VTU path:        {args.grid_vtu}")
    print(f"Output directory:     {out_dir}")
    print(f"Plot mode:            {args.mode}")
    print(f"Field vector style:   {args.efield_style}")
    print("------------------------------------------------------------------")

    # Step 1: Extract FEM Slice
    print(f"[1/3] Extracting FEM slice at z = {args.z_val:.4e} m...")
    fem_data = extract_fem_slice(args.fem_vtu, z_val=args.z_val)
    print(f"      FEM Points: {len(fem_data.x)}, Triangles: {len(fem_data.triangles)}")
    print(f"      Potential range: [{fem_data.potential.min():.2f}, {fem_data.potential.max():.2f}] V")
    print(f"      Spatial bounds: x in [{fem_data.x_mm.min():.2f}, {fem_data.x_mm.max():.2f}] mm, "
          f"y in [{fem_data.y_mm.min():.2f}, {fem_data.y_mm.max():.2f}] mm")

    # Step 2: Extract Grid Slice
    print(f"[2/3] Extracting Regular Grid slice at z = {args.z_val:.4e} m...")
    grid_data = extract_grid_slice(args.grid_vtu, z_val=args.z_val)
    print(f"      Grid resolution: {len(grid_data.x_1d)} x {len(grid_data.y_1d)} points")
    print(f"      Potential range: [{grid_data.potential.min():.2f}, {grid_data.potential.max():.2f}] V")
    print(f"      Electric field |E| max: {grid_data.E_mag_kV_m.max():.2f} kV/m ({grid_data.E_mag.max():.2e} V/m)")
    print(f"      Spatial bounds: x in [{grid_data.x_1d_mm.min():.2f}, {grid_data.x_1d_mm.max():.2f}] mm, "
          f"y in [{grid_data.y_1d_mm.min():.2f}, {grid_data.y_1d_mm.max():.2f}] mm")

    # Step 3: Render and Save Plots
    print("[3/3] Generating visualizations...")
    saved_files = []

    # Combined Dashboard
    if args.mode in ("all", "dashboard"):
        print("      Generating combined dashboard...")
        fig_dash = plot_combined_dashboard(
            fem_data,
            grid_data,
            efield_style=args.efield_style,
        )
        if not args.no_save:
            dash_path = out_dir / f"slice_dashboard_{z_str}.{fmt}"
            save_figure(fig_dash, dash_path, dpi=args.dpi)
            saved_files.append(dash_path)

    # Standalone FEM Potential
    if args.mode in ("all", "fem"):
        print("      Generating FEM potential plot...")
        fig_fem, _ = plot_fem_potential(
            fem_data,
            cmap=args.cmap_potential,
            show_mesh=args.show_mesh,
        )
        if not args.no_save:
            fem_path = out_dir / f"fem_potential_{z_str}.{fmt}"
            save_figure(fig_fem, fem_path, dpi=args.dpi)
            saved_files.append(fem_path)

    # Standalone Grid Potential
    if args.mode in ("all", "grid-potential"):
        print("      Generating Interpolated Grid potential plot...")
        fig_gpot, _ = plot_grid_potential(
            grid_data,
            cmap=args.cmap_potential,
        )
        if not args.no_save:
            gpot_path = out_dir / f"grid_potential_{z_str}.{fmt}"
            save_figure(fig_gpot, gpot_path, dpi=args.dpi)
            saved_files.append(gpot_path)

    # Standalone Grid Electric Field
    if args.mode in ("all", "grid-efield"):
        print("      Generating Interpolated Grid electric field plot...")
        fig_ef, _ = plot_grid_electric_field(
            grid_data,
            style=args.efield_style,
            cmap=args.cmap_efield,
            overlay_equipotentials=True,
        )
        if not args.no_save:
            ef_path = out_dir / f"grid_electric_field_{z_str}.{fmt}"
            save_figure(fig_ef, ef_path, dpi=args.dpi)
            saved_files.append(ef_path)

    # Standalone 1D Profiles
    if args.mode in ("all", "cuts"):
        print("      Generating 1D potential profile cuts...")
        fig_cuts, _ = plot_1d_cuts(grid_data, fem_data=fem_data)
        if not args.no_save:
            cuts_path = out_dir / f"potential_cuts_{z_str}.{fmt}"
            save_figure(fig_cuts, cuts_path, dpi=args.dpi)
            saved_files.append(cuts_path)

    print("------------------------------------------------------------------")
    if saved_files:
        print(f"Successfully generated and saved {len(saved_files)} plot(s):")
        for f in saved_files:
            print(f"  -> {f.resolve()}")
    else:
        print("Plots generated (saving skipped via --no-save).")

    if args.show:
        print("Displaying interactive plot window(s)... (Close windows to exit)")
        plt.show()

    print("==================================================================")


if __name__ == "__main__":
    main()
