# Tutorial: first runs with PercoFlow

This tutorial uses the two synthetic inputs of this folder. They are test inputs, not data. It takes about five minutes on
one core. Every number shown below was produced with PercoFlow 1.0.0 in the environment of `requirements-lock.txt`, and
`tests/test_regression.py` checks the main ones.

| file | what it is |
|---|---|
| `make_examples.py` | Regenerates the two inputs deterministically (fixed seed). |
| `synthetic_rock.png` | Image of 200 × 240 px. The pores are pure green; the solid is random discs of dark grey covering 40 % of the image. One small pore is enclosed in a solid block on purpose, so that the percolation filter has something to remove. |
| `obstacles.msh` | Gmsh triangle mesh of a 4 × 1 channel with three circular obstacles. Its physical groups are FLUID, INLET (x = 0), OUTLET (x = 4) and WALLS. |
| `config_default.yaml` | Every option with its default value and a one-line explanation. |

## 1. An image, flow from left to right

```bash
percoflow image examples/synthetic_rock.png -o runs/rock --level 3
```

```
k_D = 0.305888   (Q = 0.302659, dP = 0.99934, L = 202, W = 200)
OpenFOAM cross-check: full domain, median |FV/FE - 1| = 1.395%
```

How to read it:
* **The window.** The image is 240 px wide. A band of 8 % (19 px) is removed at the inlet side and at the outlet side,
  which leaves a sample of length L = 202 px and width W = 200 px.
* **The flow problem.** It is pressure-driven, with p_in = 1 and μ = 1. ΔP = 0.99934 is measured between the probe lines,
  one pixel inside each reservoir, so it is slightly below p_in. Q is the integrated flux.
* **The permeability.** k_D = μ Q L / (W ΔP) = 0.3059 px². With a pixel size of a metres, k = 0.3059 a² m², which
  `percoflow.descriptors.to_si` converts to m², darcy and millidarcy.
* **The cross-check.** OpenFOAM solved the same domain on the pixel staircase. Its result differs from the smooth-contour
  finite-element solution by 1.4 %; this is the size of the wall-representation effect for this image.

Open `runs/rock/result.json` to see more:

| key | value | meaning |
|---|---|---|
| `percolation.clusters_isolated` | 2 | two pore clusters (including the enclosed one) were removed before meshing |
| `percolation.connected_fraction` | 0.964 | 96.4 % of the pore area of the window carries flow |
| `critical_throat_width_px` | 2.8 | every left-to-right path crosses a throat at most 2.8 px wide |
| `solver.n_triangles` | 5880 | mesh level 3 |
| `solver.residual` | ~4e-15 | sparse LU solves the linear system to round-off |
| `solver.energy_identity_relerr` | ~8e-14 | viscous dissipation equals p_in Q, as it must |
| `descriptors.A50` | 0.105 | 10.5 % of the pore area carries half of the flow |

The files of the run directory are described in [docs/INPUTS_OUTPUTS.md](../docs/INPUTS_OUTPUTS.md).
`runs/rock/fields.png` shows the pressure and the speed, and `runs/rock/fields.vtu` opens in ParaView.

## 2. The same image, flow from top to bottom

```bash
percoflow image examples/synthetic_rock.png -o runs/rock_tb --direction tb --level 3
```

```
k_D = 6.24304   (Q = 8.5746, dP = 0.961425, L = 168, W = 240)
OpenFOAM cross-check: full domain, median |FV/FE - 1| = 0.460%
```

The permeability is about 20 times larger than from left to right. The critical throat explains why: it is 15.6 px wide
from top to bottom but 2.8 px from left to right, and the conductance of a throat grows with the cube of its width. Run both
principal directions for every sample, and report the anisotropy ratio max(k)/min(k).

## 3. A mesh you made yourself

```bash
percoflow mesh examples/obstacles.msh -o runs/obstacles
```

```
k_D = 0.0111712   (Q = 0.00279264, dP = 0.999941, L = 4.0, W = 1.0)
OpenFOAM cross-check: FV/FE - 1 = +1.841%
```

Without reservoirs the pressures are taken on the INLET and OUTLET curves, and L and W come from the bounding box (or from
`--length`, `--width`). OpenFOAM solves the same triangles extruded to prisms. Its cell-centred scheme is less accurate than
Taylor–Hood elements on the same mesh, so a difference of a few per cent here is finite-volume discretization error, and it
shrinks when the mesh is refined. `make_examples.py` shows how to name the physical groups with the Gmsh Python API.

## 4. Check a run later

```bash
percoflow verify runs/rock        # {"ok": true, "mismatch": [], "n_files": ...}
```

The command recomputes the SHA-256 of every output file listed in `runs/rock/manifest.json`. The manifest also records the
configuration, the hash of the input image, the versions of PercoFlow, Python, the libraries and OpenFOAM, the git commit,
the host and the times.

## 5. Change options

```bash
percoflow levels                                                  # mesh presets 1 (coarse) ... 5 (fine)
percoflow image my.png -o runs/my --level 4 --tolerance 0.125     # finer mesh, tighter contours
percoflow image my.png -o runs/my --config examples/config_default.yaml
```

A convergence study repeats a run at levels 1–5 and with the contour tolerance halved. Report k_D at the level where it
stops changing within the accuracy you need.
