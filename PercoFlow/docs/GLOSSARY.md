# Glossary

Short definitions of the terms used in the code and in the documentation, in alphabetical order.

| term | meaning | where |
|---|---|---|
| **A50, A90** | Smallest fraction of the pore area that carries 50 % (90 %) of the flow, with triangles ranked from fastest to slowest. A50 = 0.5 means uniform transport; small values mean flow concentrated in a few channels. | `descriptors.py` |
| **band** | Strip removed at the inlet and at the outlet side of the image before the analysis (`domain.band_fraction`, default 8 % of the side). | `imaging.py` |
| **χ_c (connected fraction)** | Pore area kept by the percolation filter divided by the pore area of the window. | `percolation.py` |
| **contour tolerance ε** | Douglas–Peucker tolerance: no wall of the simplified polygon is farther than ε pixels from the pixel boundary. | `geometry.py` |
| **critical throat width w_c** | Width of the widest constriction that every inlet–outlet path must cross (twice the largest wall distance r for which the pixels farther than r still percolate). | `percolation.py` |
| **cross-check** | Independent OpenFOAM finite-volume solution of the same domain, run by default with every computation. | `benchmarks/openfoam/` |
| **dead end** | Pore cluster that touches only one of the two planes; it cannot carry a net flow. | `percolation.py` |
| **D_h (hydraulic diameter)** | 4 × fluid area / wetted wall length of the sample. | `descriptors.py` |
| **energy identity** | For the pressure-driven problem the viscous dissipation equals p_in Q exactly, also after discretization; its relative error is stored as a check of the assembly. | `stokes.py` |
| **face-connected** | Two pixels sharing an edge. Finite volumes need it, because pixels that touch only at a corner have no common face. | `benchmarks/openfoam/pixel.py` |
| **flux-weighted pressure** | Mean pressure on a line weighted by the normal velocity, P̄ = ∫p u ds / ∫u ds. With it, Q ΔP equals the mechanical power between the two lines. | `permeability.py` |
| **isolated cluster** | Pore cluster that touches neither plane. It would make the Stokes system singular, so it is removed. | `percolation.py` |
| **k_D** | Darcy permeability, μ Q L_s / (W ΔP), in px² for images. Multiply by the square of the pixel size for m². | `permeability.py` |
| **L_s, W** | Sample length along the flow (reservoirs excluded) and full transverse width of the window. | `permeability.py` |
| **manifest** | `manifest.json` of a run directory: configuration, hashes, versions, times, notes. | `provenance.py` |
| **mesh level** | Preset 1 (coarse) to 5 (fine) of the four parameters of the mesh-size field. | `config.py` |
| **Mini element** | P1 velocity enriched with a cubic bubble, P1 pressure. Cheaper, approaches k from below; used for screening. | `stokes.py` |
| **MINRES** | Iterative solver for symmetric indefinite systems; less memory than LU. PercoFlow always checks its true residual. | `stokes.py` |
| **percolating cluster** | Pore cluster that joins the inlet and the outlet planes (8-neighbour connectivity). Only these are meshed. | `percolation.py` |
| **probe line** | Line one pixel inside the inlet or the outlet reservoir where the pressure entering k_D is measured. | `permeability.py` |
| **reservoir** | Fully open strip of width b = 0.1 L_s added before the inlet and after the outlet. The pressure is imposed on its outer face, away from the sample. | `geometry.py` |
| **Taylor–Hood element** | P2 velocity, P1 pressure. The reference element; it contains Poiseuille flow exactly. | `stokes.py` |
| **true residual** | ‖K x − b‖ / ‖b‖ of the unpreconditioned system, computed after every solve. | `stokes.py` |
