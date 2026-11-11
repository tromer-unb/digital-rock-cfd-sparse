# The paper explained

A didactic walk-through of what the article does, for a reader who knows some fluid mechanics but not digital rock physics.
Equations are given in their simplest form; the full derivations are in the Supplementary Notes of the paper. Each section
ends with the place in this repository where the corresponding code or data can be found.

## 1. Permeability and Darcy's law

At the scale of a rock sample, slow flow obeys Darcy's law: the mean (superficial) velocity U is proportional to the pressure
gradient, U = (k/μ) ΔP/L. The constant k, the **permeability**, has units of area and depends only on the geometry of the pore
space. Digital rock physics computes k by solving the flow equations directly inside an image of the pores.

## 2. The flow problem

Pore-scale flow in rocks is slow and viscous, so the inertial terms vanish and the Navier–Stokes equations reduce to the
**Stokes equations**: −μ∇²u + ∇p = 0 and ∇·u = 0, with u = 0 on the pore walls. The problem is linear: twice the pressure
difference gives twice the flow, so k does not depend on the pressure we impose.

## 3. From pixels to a computational domain

1. **Segmentation** (done before this work): each pixel is pore (green) or solid.
2. **Cropping**: bands of 8 % of the image side are removed where the flow enters and leaves, to define clean inlet and outlet
   planes.
3. **Percolation filter**: pore clusters are labelled; only those touching both planes are kept. *Why it matters*: a closed
   cluster lets the pressure take any constant value with zero velocity, so the matrix of the discrete problem becomes
   singular; a cluster open only to the inlet cannot accept an imposed inflow at all. Sample R1 has no connected path and is
   rejected; R2 percolates only vertically (χ_c = 53.8 % of its pore area is connected); R3 in both directions.
4. **Reservoirs**: open strips (10 % of the sample length) are added before and after the sample. The pressure is imposed on
   their outer faces, so no artificial velocity profile is forced into the pores.
5. **Contours and mesh**: the pixel boundary becomes polygons (marching squares), simplified within 0.25 px, and Gmsh fills
   them with triangles that are small near walls and in narrow throats.

*In this repository:* `PercoFlow/src/percoflow/imaging.py`, `percolation.py`, `geometry.py`, `meshing.py`;
cluster census in `data/results/benchmarks/c3_wellposedness.json`.

## 4. Solving

Velocity and pressure are approximated with **mixed finite elements** that are mathematically stable for this saddle-point
problem: Taylor–Hood (quadratic velocity, linear pressure) as the accurate reference and Mini (linear + bubble) as a cheaper
screening option. The resulting sparse linear system (up to 2.4 million unknowns for a full 4096² image) is solved by sparse
LU in one to two minutes, or by MINRES with less memory. A subtle point: MINRES's own stopping test can accept a solution
whose true residual is 10⁻³, so the workflow checks the true residual before accepting a solution.

*In this repository:* `PercoFlow/src/percoflow/stokes.py`; MINRES study in `data/results/benchmarks/c4_minres_p1.json`.

## 5. Measuring k

k_D = μ Q L_s / (W ΔP), where Q is the computed flux, L_s the sample length, W the full width of the window, and ΔP the
difference of **flux-weighted** pressures one pixel inside each reservoir. Flux weighting is not arbitrary: Q·ΔP then equals
the power dissipated by viscosity in the sample, so k_D is the permeability of the uniform Darcy medium that dissipates the
same power at the same flux.

*In this repository:* `PercoFlow/src/percoflow/permeability.py`; the stored value of every run is the key `kD_iface_flux`
of its JSON file in `data/results/fix_round/`.

## 6. Is the solver right? Verification and validation

* **Poiseuille channel**: the exact parabolic solution lies inside the Taylor–Hood space, so the computed k must equal the
  exact value on any mesh. It does, to 5×10⁻⁷ — this checks the whole chain (mesh, assembly, boundary conditions, flux,
  pressure probes) at once.
* **Cylinder array**: k agrees with the classical Hasimoto series to 0.34 % where the series is valid.
* **OpenFOAM**: an independent finite-volume code, on the same cylinder cells, agrees to 0.09–0.11 %. Agreement between two
  codes that share nothing validates the solution, not just the implementation.

*In this repository:* `percoflow benchmark` (B1–B5) and `data/results/benchmarks/`; OpenFOAM cases and logs in
`data/results/fix_round/FR-G_openfoam_cell_v2_alienware/`.

## 7. The main result: definitions matter more than discretization

With the reference solver the mesh error is below 0.03 %. Four modelling definitions change k far more:

* **inlet condition** — a uniform inlet velocity (common in CFD set-ups) changed k by up to ×2.2 versus an imposed pressure;
* **contour tolerance** — simplification widens throats; 1.5 px overestimated k by up to 11 %, 2 px by 30 %;
* **window** — moving the cropping band by one per cent of the side changed k of R3-H by 79 %, because a constriction left the
  window; R2-V stops percolating below a 4 % band;
* **wall representation** — a pixel-staircase finite-volume mesh and the smooth-contour finite-element mesh of the *same*
  pixels differ by a median of 1.8 % (up to 32 %). Neither is wrong: they are two different domains built from one image.

*In this repository:* `data/results/audit_evidence/indep_sens/` (inlet), `fix_round/FR-C_geometry/` (tolerance),
`fix_round/FR-E2_window_scan_tol0.25/` (window), `fix_round/FR-I_subwindow_pairs_tile512/` (wall representation).

## 8. Reading the flow fields

* **R3-H** (k = 18.28 px²): two pore bodies joined by one narrow passage; almost all the pressure drop happens there and the
  speed peaks at 210 times the mean Darcy velocity. The critical throat is 26.9 px wide.
* **R2-V** (k = 60.90 px²): two bodies in series; the drop splits between the inlet openings (~27 %) and an elongated neck
  (~45 %).
* **R3-V** (k = 28 781 px²): a broad corridor (critical throat 792 px); the maximum speed is only 5.7 times the mean.

Porosity, connected fraction and hydraulic diameter of R3-H and R3-V are nearly identical, yet their permeabilities differ by
more than three orders of magnitude: scalar geometric descriptors cannot replace a resolved flow calculation.

*In this repository:* `data/results/fix_round/FR-H_fields_P1/` (fields and descriptors), `PercoFlow/src/percoflow/descriptors.py`.

## 9. Limits

Flow is 2D (throats closed in a section may be open in 3D); lengths are in pixels, so k is in px² (multiply by the square of
the pixel size for m²); the windows are not representative elementary volumes. A Stokes computation also cannot "confirm" a small Reynolds
number by itself: Re = ρ k D_h G/μ² depends on the pressure gradient and pixel size of the application.
