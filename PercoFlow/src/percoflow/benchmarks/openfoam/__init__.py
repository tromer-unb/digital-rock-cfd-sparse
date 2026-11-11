"""OpenFOAM cross-check: an independent finite-volume solution of the same Stokes problem.

Why: two finite-element solvers built on the same library can only verify an implementation. A code that shares neither
the discretization nor any software with the workflow validates the flow solution itself (docs/OPENFOAM.md). The
cross-check is therefore ON by default in every run and in the benchmark suite.

Modules:
    environment    find OpenFOAM, run its utilities, record the version
    polymesh       write a one-layer extruded polyMesh from 2D cells (pixels or triangles), with patch names
    case           case dictionaries (simpleFoam, SIMPLEC, laminar), run, read fields, flux-weighted pressures
    pixel          cross-check of an image domain on a pixel mesh (one hexahedron per pore pixel)
    triangles      cross-check of a user triangle mesh (one prism per triangle)
    cylinder_cell  periodic cylinder array, FV (structured O-grid) versus FE (Taylor-Hood) versus Hasimoto series
"""
from .environment import OpenFOAMNotFound, find_openfoam  # noqa: F401
