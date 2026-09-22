"""
test_downstream.py
==================
Unit and regression test suite for downstream_analysis package.
Verifies VTU slice extraction, array shapes, physics bounds, and plot rendering.
"""

import sys
import unittest
from pathlib import Path
import tempfile
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

pkg_dir = Path(__file__).resolve().parent
repo_dir = pkg_dir.parent
if str(repo_dir) not in sys.path:
    sys.path.insert(0, str(repo_dir))
if str(pkg_dir) not in sys.path:
    sys.path.insert(0, str(pkg_dir))

from downstream_analysis.slice_extractor import (
    extract_fem_slice,
    extract_grid_slice,
    resolve_file_path,
    FemSliceData,
    GridSliceData,
)
from downstream_analysis.visualizer import (
    plot_fem_potential,
    plot_grid_potential,
    plot_grid_electric_field,
    plot_combined_dashboard,
    plot_1d_cuts,
    save_figure,
)


class TestSliceExtractor(unittest.TestCase):

    def setUp(self):
        self.fem_vtu = resolve_file_path("output/vtk/fem_solution.vtu")
        self.grid_vtu = resolve_file_path("output/vtk/interpolated_grid.vtu")

    def test_fem_slice_extraction(self):
        """Verify extraction from unstructured FEM tetrahedral mesh."""
        fem_data = extract_fem_slice(self.fem_vtu, z_val=0.0)
        self.assertIsInstance(fem_data, FemSliceData)
        self.assertEqual(fem_data.z_val, 0.0)
        self.assertGreater(len(fem_data.x), 500)
        self.assertEqual(len(fem_data.x), len(fem_data.y))
        self.assertEqual(len(fem_data.x), len(fem_data.potential))
        self.assertEqual(fem_data.triangles.shape[1], 3)

        # Check physical bounds (cylinder radius r0 = 5 mm = 0.005 m)
        self.assertAlmostEqual(fem_data.x_bounds[0], -0.005, delta=2e-4)
        self.assertAlmostEqual(fem_data.x_bounds[1], 0.005, delta=2e-4)
        self.assertAlmostEqual(fem_data.y_bounds[0], -0.005, delta=2e-4)
        self.assertAlmostEqual(fem_data.y_bounds[1], 0.005, delta=2e-4)

        # Voltage bounds (-100 V to +100 V applied on radial electrodes)
        self.assertLessEqual(fem_data.potential.min(), -90.0)
        self.assertGreaterEqual(fem_data.potential.max(), 90.0)

    def test_grid_slice_extraction(self):
        """Verify extraction and restructuring from interpolated Cartesian grid."""
        grid_data = extract_grid_slice(self.grid_vtu, z_val=0.0)
        self.assertIsInstance(grid_data, GridSliceData)
        self.assertEqual(grid_data.z_val, 0.0)

        # 21x21 resolution
        self.assertEqual(len(grid_data.x_1d), 21)
        self.assertEqual(len(grid_data.y_1d), 21)
        self.assertEqual(grid_data.potential.shape, (21, 21))
        self.assertEqual(grid_data.Ex.shape, (21, 21))
        self.assertEqual(grid_data.Ey.shape, (21, 21))
        self.assertEqual(grid_data.Ez.shape, (21, 21))
        self.assertEqual(grid_data.E_mag.shape, (21, 21))

        # Check grid domain bounds (-2 mm to +2 mm)
        self.assertAlmostEqual(grid_data.x_bounds[0], -0.002, delta=1e-5)
        self.assertAlmostEqual(grid_data.x_bounds[1], 0.002, delta=1e-5)

        # Electric field physics checks:
        # At trap center (0, 0), |E| should approach minimum
        mid_idx = 10
        center_emag = grid_data.E_mag[mid_idx, mid_idx]
        max_emag = grid_data.E_mag.max()
        self.assertLess(center_emag, 0.01 * max_emag)  # Center field < 1% of edge field


class TestVisualizer(unittest.TestCase):

    def setUp(self):
        self.fem_vtu = resolve_file_path("output/vtk/fem_solution.vtu")
        self.grid_vtu = resolve_file_path("output/vtk/interpolated_grid.vtu")
        self.fem_data = extract_fem_slice(self.fem_vtu, z_val=0.0)
        self.grid_data = extract_grid_slice(self.grid_vtu, z_val=0.0)

    def test_plot_fem_potential(self):
        fig, ax = plot_fem_potential(self.fem_data, show_mesh=True)
        self.assertIsInstance(fig, plt.Figure)
        plt.close(fig)

    def test_plot_grid_potential(self):
        fig, ax = plot_grid_potential(self.grid_data)
        self.assertIsInstance(fig, plt.Figure)
        plt.close(fig)

    def test_plot_grid_electric_field_styles(self):
        for style in ["streamlines", "quiver", "both"]:
            fig, ax = plot_grid_electric_field(self.grid_data, style=style, overlay_equipotentials=True)
            self.assertIsInstance(fig, plt.Figure)
            plt.close(fig)

    def test_plot_1d_cuts(self):
        fig, ax = plot_1d_cuts(self.grid_data, fem_data=self.fem_data)
        self.assertIsInstance(fig, plt.Figure)
        plt.close(fig)

    def test_plot_combined_dashboard_and_save(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out_file = Path(tmpdir) / "test_dash.png"
            fig = plot_combined_dashboard(self.fem_data, self.grid_data, save_path=out_file, dpi=100)
            self.assertTrue(out_file.exists())
            self.assertGreater(out_file.stat().st_size, 1000)
            plt.close(fig)


if __name__ == "__main__":
    unittest.main()
