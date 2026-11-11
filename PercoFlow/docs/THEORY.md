# Theory: from a segmented image to a Darcy permeability

This page explains every step of PercoFlow with the equations it implements and the reason for each choice. Each section
names the module that implements it, so that the equations can be checked against the code line by line.

## 1. Physical model (`stokes.py`)

Steady, incompressible, single-phase creeping flow of a Newtonian fluid (dynamic viscosity μ) in the 2D fluid domain Ω:

```
-μ ∇²u + ∇p = 0,   ∇·u = 0        in Ω
u = 0                              on Γ_w   (pore walls and lateral edges)
μ ∂u/∂n − p n = −p_in n            on Γ_in  (pressure-driven inlet)
μ ∂u/∂n − p n = 0                  on Γ_out (traction-free outlet)
```

The open-boundary conditions are the natural ("do-nothing") conditions of the Laplacian form of the viscous term. For the
uniform parallel flow inside an open reservoir they coincide with an imposed pressure. The problem is linear: doubling `p_in`
doubles `u`, so the permeability does not depend on `p_in` (default 1) and **the computation does not fix a Reynolds number**
(Section 7).

## 2. Why the percolation filter is required (`percolation.py`)

Connected-component labelling (8-neighbour) of the window classifies each pore cluster as *percolating* (touches both
planes), *dead end* (one plane) or *isolated* (none). Only percolating clusters are meshed:

* an isolated cluster admits `u = 0, p = const` as a non-zero solution of the homogeneous problem ⇒ one kernel vector of
  the saddle-point matrix per isolated cluster ⇒ singular system;
* under a velocity inlet a dead end receives inflow with nowhere to discharge it ⇒ `∫ ∇·u = −U₀|Γ_i| ≠ 0` ⇒ no solution.

The connected pore fraction `χ_c = A_kept / A_pore` and the critical throat width `w_c` (largest wall distance for which the
pixels farther than it still percolate, ×2, found by bisection) are reported with every run.

## 3. Reservoirs and the permeability definition (`geometry.py`, `permeability.py`)

The sample (length `L_s`, transverse width `W`) is padded with fully open reservoirs of width `b = round(0.1 L_s)` on the
inlet and outlet sides. The pressure is imposed on their outer faces and measured **one pixel inside** them:

```
k_D = μ Q L_s / (W ΔP),   ΔP = P̄_in − P̄_out,   P̄ = ∫_S p u_s ds / ∫_S u_s ds
```

* `Q` integrated flux (not a nominal one), `W` full transverse width (not the open length of the inlet, which would give an
  interstitial quantity larger by 1/(open fraction)), `S` probe line sampled every half pixel.
* Flux weighting makes `Q ΔP` equal to the mechanical power delivered between the probe lines (up to a small viscous flux
  through them), so `k_D` is the permeability of the homogeneous Darcy medium that dissipates the same power.
* The discrete solution satisfies the global energy identity `Φ = μ∫|∇u|² = p_in Q` exactly; it is stored with every
  solution (`energy_identity_relerr`, ~1e-14) as an algebraic check.
* Measuring `ΔP` at the outer faces would add the access resistance of the reservoirs. This share is reported as
  `reservoir_pressure_share`; it grows with the permeability of the window.

Why not a uniform inlet velocity: it forces a plug profile into every pore opening and creates singularities where the inlet
meets a wall. The result then depends on how the inlet openings are resolved, and therefore on the contour tolerance. The
velocity inlet is kept (`solver.inlet: velocity`) only to compare with workflows that use it.

## 4. Geometry and mesh (`geometry.py`, `meshing.py`)

Marching squares at level ½ (pores fully connected) gives polygons through the pixel edges; Douglas–Peucker simplification
with tolerance ε moves no wall by more than ε. A throat of width w has a conductance proportional to w³, so moving both of
its walls by ε changes it by up to ~6ε/w. **Keep ε well below one pixel (default 0.25 px)** and check, for new kinds of
images, that `k_D` does not change when ε is halved. Contours with fewer than 8 stored points (7 distinct vertices: the
first point is repeated to close the polygon) or an area below 20 px² are dropped. Nested polygons become plane surfaces with holes
(even–odd depth). Gmsh (Frontal–Delaunay, 5 smoothing steps) meshes them with the size field

```
h(d) = h_min (1 − σ) + h_max σ,   σ = 1 / (1 + exp(−(12 r − 6))),   r = clip((d − d_min)/(d_max − d_min), 0, 1)
```

with `d` the distance to the nearest wall: small elements in throats, large ones in pore bodies. Five presets
(`percoflow levels`) support convergence studies.

## 5. Discretization (`stokes.py`)

Weak form: find `(u, p)`, `u = 0` on `Γ_w`, such that for all `(v, q)`

```
μ(∇u, ∇v) − (p, ∇·v) = −p_in ⟨v·n, 1⟩_Γin,     (q, ∇·u) = 0
```

| element | velocity | pressure | use |
|---|---|---|---|
| Taylor–Hood (`p2`) | continuous P2 | continuous P1 | reference; contains Poiseuille flow exactly |
| Mini (`mini`) | P1 + cubic bubble | continuous P1 | screening; ~1.3× fewer unknowns; converges from below |

Both satisfy the inf–sup (LBB) condition with a mesh-independent constant. Assembly with scikit-fem, quadrature exact for
degree 4. Wall unknowns are removed by index slicing, which keeps the reduced matrix symmetric.

## 6. Linear solvers (`stokes.py`)

* **Sparse LU** (default): algebraic error at round-off. A full 4096² image (≈2.4 M unknowns, Taylor–Hood) solves in
  ~1–2 min and ~15 GB.
* **MINRES** with `P = diag(diag(A), diag(M_p)/μ)`: less memory, slower. MINRES stops on the backward error of the
  *preconditioned* system, which can be orders of magnitude below the true residual. `percoflow` always computes
  `‖Kx − b‖/‖b‖` and **rejects** the solution above `solver.accept_residual` (`SolverNotConverged`).

## 7. Descriptors and the Reynolds number (`descriptors.py`)

Evaluated on the sample triangles (reservoirs excluded); speeds normalized by `U = Q/W`.

* hydraulic partition: active `u_T ≥ 0.01 u_max`, dead `u_T ≤ 0.001 u_max`, intermediate otherwise (area shares);
* transport localization `A_50`, `A_90`: smallest area fractions carrying 50 % / 90 % of `Σ u_T|T|`, triangles ranked by speed;
* hydraulic diameter `D_h = 4A_f/P_w`, domain porosity, velocity amplification `u_max/U`.

**Reynolds number.** With a mean gradient `G`, `U = k_D G/μ`, so `Re = ρ U D_h/μ = ρ k_D D_h G/μ²`. In pixel units with
pixel size `a`: `Re = ρ k̂ D̂ a³ G/μ²`. Creeping flow requires `G ≪ μ²/(ρ k̂ D̂ a³)`;
`percoflow.descriptors.reynolds_bound(k, D_h, a)` returns that gradient. A Stokes solution can never "confirm" Re ≪ 1 by
itself: the pixel size and the driving gradient of the application are needed.

## 8. Independent cross-check (`benchmarks/openfoam/`)

See [OPENFOAM.md](OPENFOAM.md). The finite-element result is accompanied by a finite-volume solution of the same domain:
pixel hexahedra for images (keeps the staircase), extruded prisms for user meshes. On a geometry without pixels (the
cylinder cell, benchmark B4) the two codes agree to about 0.1 %. On images they differ by the representation of the walls
(smooth contour versus pixel staircase), and PercoFlow reports that difference instead of hiding it.
