# Downstream Analysis: 2D Planar Slicing of Paul Trap Fields

This package provides post-processing and visualization tools for the 3D Paul trap field solver (`paul_trap_sim`). It extracts planar 2D slices (by default at the trap center, $z = 0$) from the simulation output VTU files (`fem_solution.vtu` and `interpolated_grid.vtu`), producing publication-quality Matplotlib figures of the **electric potential** and **electric field**.

---

## 1. Directory Structure

```text
downstream_analysis/
├── .venv/                   # Isolated Python virtual environment
├── requirements.txt         # Pinned Python package dependencies
├── slice_extractor.py       # VTU loading and planar slice extraction
├── visualizer.py            # Matplotlib 2D plotting functions
├── slice_plotter.py         # Command-line driver script
├── test_downstream.py       # Unit and regression test suite
├── plots/                   # Default export directory for figures
│   ├── slice_dashboard_z0.png
│   ├── fem_potential_z0.png
│   ├── grid_potential_z0.png
│   ├── grid_electric_field_z0.png
│   └── potential_cuts_z0.png
└── README.md                # Documentation (this file)
```

---

## 2. Environment Setup

The `.venv` virtual environment is scoped solely inside `downstream_analysis/`:

```bash
# Navigate to the downstream analysis directory
cd downstream_analysis

# Activate the scoped virtual environment
source .venv/bin/activate

# (Optional) Verify or reinstall dependencies
pip install -r requirements.txt
```

---

## 3. Quick Start & CLI Usage

### Basic Run
Extracts the slice at $z = 0$, processes both the FEM and interpolated grid solutions, and exports 5 publication-ready plots into `downstream_analysis/plots/`:

```bash
python slice_plotter.py
```

### Interactive Display
To inspect plots interactively in a window rather than just saving to disk:

```bash
python slice_plotter.py --show
```

### Advanced Options

| Flag | Type | Default | Description |
| :--- | :--- | :--- | :--- |
| `--fem-vtu` | `str` | `output/vtk/fem_solution.vtu` | Path to FEM VTU file |
| `--grid-vtu` | `str` | `output/vtk/interpolated_grid.vtu` | Path to regular grid VTU file |
| `--z-val` | `float` | `0.0` | Planar slice height along $z$ (meters) |
| `--output-dir` | `str` | `downstream_analysis/plots` | Target directory for generated figures |
| `--mode` | `choice` | `all` | `all`, `dashboard`, `fem`, `grid-potential`, `grid-efield`, `cuts` |
| `--efield-style` | `choice` | `streamlines` | `streamlines`, `quiver`, or `both` |
| `--cmap-potential`| `str` | `coolwarm` | Colormap for electric potential |
| `--cmap-efield` | `str` | `viridis` | Colormap for electric field magnitude |
| `--format` | `choice` | `png` | Output format: `png`, `pdf`, `svg` |
| `--dpi` | `int` | `300` | Resolution for exported raster figures |
| `--show-mesh` | `flag` | `False` | Overlay FEM triangular wireframe edges |
| `--show` | `flag` | `False` | Display interactive Matplotlib window |
| `--no-save` | `flag` | `False` | Disable file export (useful with `--show`) |

#### Examples:
```bash
# Generate PDF figures with quiver vector arrows
python slice_plotter.py --efield-style quiver --format pdf

# Slice at z = 1 mm with custom colormaps
python slice_plotter.py --z-val 0.001 --cmap-potential RdBu_r --cmap-efield magma

# Overlay the FEM triangular mesh
python slice_plotter.py --mode fem --show-mesh
```

---

## 4. Python API Usage

You can also import the modules into custom Python scripts or Jupyter notebooks:

```python
from downstream_analysis.slice_extractor import extract_fem_slice, extract_grid_slice
from downstream_analysis.visualizer import plot_combined_dashboard, plot_grid_electric_field
import matplotlib.pyplot as plt

# Extract slice at z = 0
fem_data = extract_fem_slice("output/vtk/fem_solution.vtu", z_val=0.0)
grid_data = extract_grid_slice("output/vtk/interpolated_grid.vtu", z_val=0.0)

# Create 4-panel dashboard
fig = plot_combined_dashboard(fem_data, grid_data, efield_style="streamlines")
plt.show()
```

---

## 5. Physical Interpretation of Plots

For standard voltages ($V_1 = +100\,\text{V}, V_2 = -100\,\text{V}, V_3 = +100\,\text{V}, V_4 = -100\,\text{V}, V_5 = +50\,\text{V}, V_6 = +50\,\text{V}$):

1. **FEM Potential (Full Domain):**
   - Shows the entire cross section of the cylindrical apparatus ($r \le r_0 = 5\,\text{mm}$).
   - Highlights the 4 radial quadrant boundary conditions ($+100\,\text{V}$ on $x$-axis electrodes, $-100\,\text{V}$ on $y$-axis electrodes) and the central smooth field transition.
2. **Interpolated Grid Potential (Central Zone):**
   - Zooms into the central trapping region ($-2\,\text{mm} \le x, y \le 2\,\text{mm}$).
   - Features labeled equipotential hyperbolas with saddle point at the origin.
3. **Electric Field ($\vec{E} = -\nabla V$):**
   - Background colormap indicates field magnitude $|\vec{E}_{xy}|$ in $\text{kV/m}$.
   - Streamlines point inward along the confining $x$-axis and outward along the anti-confining $y$-axis toward the negative electrodes.
   - Field magnitude reaches zero ($|\vec{E}| \approx 0$) at the trap center (harmonic quadrupole null point).
4. **1D Central Potential Well Profiles:**
   - Displays $V(x, y=0)$ (parabolic trapping well, positive curvature) and $V(x=0, y)$ (anti-confining barrier, negative curvature), clearly illustrating $\nabla^2 V = \frac{\partial^2 V}{\partial x^2} + \frac{\partial^2 V}{\partial y^2} = 0$.

---

## 6. Running Tests

```bash
cd downstream_analysis
./.venv/bin/python test_downstream.py
```
