"""PercoFlow: percolation-aware flow solver for the permeability of two-dimensional porous images and meshes.

The name: *Perco* for percolation (only the pore clusters that join the inlet and the outlet are kept, before any mesh is
made) and *Flow* for the steady Stokes flow solved on them with stable mixed finite elements and checked with OpenFOAM.

What happens to an image, step by step (module in brackets; docs/THEORY.md gives the equations and the reasons):

    image ──► boolean pore mask ......................................... [imaging.py]
          ──► remove the inlet and outlet bands ......................... [imaging.py]
          ──► keep only the percolating pore clusters ................... [percolation.py]
          ──► add open reservoirs before the inlet and after the outlet . [geometry.py]
          ──► contours (marching squares + Douglas-Peucker) ............. [geometry.py]
          ──► triangle mesh graded towards the walls (Gmsh) ............. [meshing.py]
          ──► Stokes flow, Taylor-Hood or Mini elements, LU or MINRES ... [stokes.py]
          ──► Darcy permeability k_D with flux-weighted pressures ....... [permeability.py]
          ──► pore-scale descriptors .................................... [descriptors.py]
          ──► OpenFOAM finite-volume solution of the same domain ........ [benchmarks/openfoam/]
          ──► result.json and manifest.json (hashes, versions, options) . [pipeline.py, provenance.py]

A user mesh (Gmsh .msh with physical groups FLUID, INLET, OUTLET, WALLS) enters directly at the Stokes step.

Public entry points:
    percoflow.pipeline.run_image(...)    image -> run directory
    percoflow.pipeline.run_mesh(...)     mesh  -> run directory
    percoflow.benchmarks.run_all(...)    verification and validation suite (OpenFOAM included)
    command line: ``percoflow --help``
"""

__version__ = "1.0.0"

#: Unit convention. All lengths are in the units of the input (pixels for images). With the defaults
#: mu = 1 and p_in = 1 the permeability is returned in (length unit)^2, e.g. px^2; multiply by the squared
#: pixel size to obtain m^2.
UNITS_NOTE = "lengths in input units (px for images); k_D in (unit)^2; multiply by (pixel size)^2 for SI"
