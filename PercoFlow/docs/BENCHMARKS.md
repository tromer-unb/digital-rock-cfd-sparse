# Verification and validation suite

```bash
percoflow benchmark -o runs/bench          # quick suite, about 15 min on one core (OpenFOAM included)
percoflow benchmark -o runs/bench --full   # finer meshes and more solid fractions (hours)
```

The command writes `runs/bench/benchmark_report.json` with every number, the criterion and PASS/FAIL, and prints
`STATUS: PASS` only if every criterion passes. The OpenFOAM benchmarks are always part of the suite. With `--no-openfoam`
the status becomes `INCOMPLETE`, never `PASS`.

**Verification** asks whether the equations are solved correctly; it compares against exact or analytic solutions.
**Validation** of the flow solution compares against an independent code that shares no discretization and no software with
PercoFlow. The criteria below are written in each module's docstring and were fixed **before** any result was looked at.

## The five benchmarks

| id | what is tested | reference | criterion | module |
|---|---|---|---|---|
| B1 | plane Poiseuille flow through the **whole image pipeline** (channel images h = 20 and 40 px, L/h = 20, mesh levels 1, 3, 5): reservoirs, contours, Gmsh, Stokes, probe lines | exact value of the k_D definition, (h²/12) L_s/(L_s+2) | Taylor–Hood relative error < 1e-5 | `benchmarks/poiseuille.py` |
| B2 | periodic square array of cylinders, finite elements, analytic geometry (no pixels) | Hasimoto / Sangani–Acrivos series | ≤ 0.5 % for solid fraction c ≤ 0.15 | `benchmarks/cylinder_array.py` |
| B3 | OpenFOAM pixel route and its post-processing, on an open channel of 16 × 80 px with reservoirs | exact cell-centred value (h²+2)/12 · L_s/(L_s+2) | relative error < 1e-6 | `benchmarks/__init__.py` |
| B4 | OpenFOAM on a structured O-grid cylinder cell (c = 0.1, 0.2; 1/h = 25, 50, 100) against the finite-element value of B2 | PercoFlow FE value (Taylor–Hood, h = 0.01) | FV converged within 1 % and FV − FE within 2 % | `benchmarks/openfoam/cylinder_cell.py` |
| B5 | OpenFOAM pixel FV against contour FE on a synthetic porous window (160 × 160 px, seed 7) | none: the two represent the walls differently | both solvers converge and checkMesh reports OK; the difference is reported | `benchmarks/__init__.py` |

Why each one:
* **B1** exercises every step of the image route at once. Poiseuille flow lies inside the Taylor–Hood spaces, so any error in
  assembly, boundary conditions, flux integration or pressure probing shows up as a departure from round-off.
* **B2** tests the solver on a geometry that does not come from pixels, against an analytic series.
* **B3** checks that the OpenFOAM case and the probe-line post-processing reproduce an exact discrete solution.
* **B4** is the validation: the same physical problem solved by an independent code.
* **B5** shows on a porous geometry what the cross-check of every image run reports.

## Results of this release

The results below are for PercoFlow 1.0.0 with OpenFOAM v1912 and the quick suite. The complete report is
[reference_results/benchmark_report.json](reference_results/benchmark_report.json).

| id | result | status |
|---|---|---|
| B1 | max Taylor–Hood error 4.95e-7. Mini lies below the exact value in every case (between −19.8 % and −0.6 %), as expected for an element that cannot represent the parabola. | PASS |
| B2 | max 0.339 % for c ≤ 0.15. The departure grows above c = 0.15, where the truncated series loses accuracy. | PASS |
| B3 | relative error −3.9e-13 (round-off) | PASS |
| B4 | FV change between the two finest meshes 0.21 %. FV − FE = −0.111 % (c = 0.1) and −0.125 % (c = 0.2). | PASS |
| B5 | both converge, checkMesh OK. FV − FE = +1.65 % (pixel staircase against smooth contour). | PASS |

The quick suite took 15 minutes on one thread (2026-11-06).

## Supporting tests

`pytest -q` runs 19 tests, described in [../tests/README.md](../tests/README.md):
- configuration round trip and typo rejection;
- band cropping, percolation classification, critical throat width, contour nesting, reservoirs;
- Poiseuille exactness and the energy identity;
- flow-direction invariance (a top-to-bottom run equals the left-to-right run of the rotated image);
- the MINRES true-residual gate (a solution above `accept_residual` must raise);
- the face-connectivity rule;
- end-to-end image and mesh runs with the OpenFOAM cross-check, and the recording of a disabled cross-check;
- manifest verification;
- the regression values of the example inputs.
