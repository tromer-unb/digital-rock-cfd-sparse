"""Run configuration.

Every option that can change a number is a field of one of the dataclasses below, and the complete configuration is
written to ``manifest.json`` with each run (see ``provenance.py``), so a result can always be traced to the options that
produced it. The defaults define the reference protocol of PercoFlow:

* pressure-driven flow between open reservoirs of width b = 0.1 L_s,
* Darcy permeability k_D = mu Q L_s / (W dP) with flux-weighted pressures one pixel inside each reservoir,
* contour tolerance 0.25 px, boundary band 8 % of the image side, Taylor-Hood elements, sparse LU,
* mesh level 2 (change ``mesh.level`` to 1-5, see ``MESH_LEVELS``).

A configuration can be written as YAML or JSON with the same nesting, e.g.::

    domain:  {direction: lr, band_fraction: 0.08, reservoir_fraction: 0.10}
    mesh:    {level: 3, contour_tolerance: 0.25}
    solver:  {element: p2, method: lu}
    openfoam: {enabled: true}
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field, fields, is_dataclass
from pathlib import Path
from typing import Optional

#: Mesh-size presets (h_min, h_max, d_min, d_max) in pixels, from coarse (1) to fine (5), for mesh-convergence studies.
#: Equation of the size field: h(d) = h_min (1 - s) + h_max s,  s = sigmoid of the wall distance d between d_min and d_max.
MESH_LEVELS = {
    1: (12.0, 70.0, 12.0, 140.0),
    2: (8.0, 50.0, 8.0, 100.0),
    3: (4.0, 30.0, 4.0, 60.0),
    4: (2.0, 15.0, 2.0, 30.0),
    5: (1.5, 12.0, 1.5, 24.0),
}


@dataclass
class ImageOptions:
    """How an image becomes a boolean pore mask (True = pore)."""
    #: "green": pore = pure-green pixels (G >= green_min and G exceeds R and B by green_delta), the usual colour-labelled segmentation;
    #: "threshold": grey value above the mid-range (or above ``threshold`` if given); "nonzero": any non-zero pixel.
    pore_rule: str = "green"
    green_min: int = 80
    green_delta: int = 30
    threshold: Optional[float] = None
    #: swap pore and solid after the rule has been applied (for images where pores are dark, for instance)
    invert: bool = False
    #: pore objects smaller than this many pixels (8-connected) are removed; 0 disables
    min_object_size: int = 20


@dataclass
class DomainOptions:
    """Which part of the mask is analysed and how flow is driven through it."""
    #: "lr": flow from left to right (x); "tb": flow from top to bottom (-y)
    direction: str = "lr"
    #: bands of this fraction of the image side are removed at the inlet and outlet sides (default 0.08). The bands keep
    #: artificial connections along the image border out of the analysis and give every sample identical planes.
    band_fraction: float = 0.08
    #: absolute band width in pixels; overrides band_fraction when given
    band_px: Optional[int] = None
    #: "all": keep every cluster that joins both planes; "largest": keep only the largest one
    keep_percolating: str = "all"
    #: reservoir width b as a fraction of the sample length L_s (0 = no reservoirs, only for velocity-driven runs)
    reservoir_fraction: float = 0.10


@dataclass
class MeshOptions:
    """Contour extraction and Gmsh meshing."""
    #: preset 1..5 from MESH_LEVELS; ignored when all four explicit sizes below are given
    level: int = 2
    h_min: Optional[float] = None
    h_max: Optional[float] = None
    d_min: Optional[float] = None
    d_max: Optional[float] = None
    #: Douglas-Peucker tolerance in pixels. Keep it well below one pixel: moving a wall by eps changes the conductance of a
    #: throat of width w by up to ~6 eps / w (docs/THEORY.md, Section 4). Check that k does not change when eps is halved.
    contour_tolerance: float = 0.25
    #: contours with fewer points or a smaller area (px^2) are discarded. A closed contour stores its first vertex twice,
    #: so the default 8 points means 7 distinct vertices.
    min_contour_points: int = 8
    min_contour_area: float = 20.0
    #: Gmsh 2D algorithm (6 = Frontal-Delaunay) and number of Laplacian smoothing steps
    algorithm: int = 6
    smoothing: int = 5
    #: polygon edges within this distance of an outer face are labelled INLET/OUTLET, all others WALLS
    face_tolerance: float = 0.75

    def sizes(self):
        """Return (h_min, h_max, d_min, d_max), from the explicit values or from the level preset."""
        explicit = (self.h_min, self.h_max, self.d_min, self.d_max)
        if all(v is not None for v in explicit):
            return tuple(float(v) for v in explicit)
        if self.level not in MESH_LEVELS:
            raise ValueError(f"mesh.level must be one of {sorted(MESH_LEVELS)}, got {self.level}")
        return MESH_LEVELS[self.level]


@dataclass
class SolverOptions:
    """Discretization, linear solver and physical constants."""
    #: "p2": Taylor-Hood P2/P1 (reference); "mini": P1-bubble/P1 (screening, fewer unknowns, underestimates k by ~0.3-2.5 %)
    element: str = "p2"
    #: "lu": sparse direct factorization; "minres": MINRES with a block-diagonal preconditioner (less memory)
    method: str = "lu"
    #: "pressure": traction -p_in n at the inlet, traction-free outlet (reference protocol);
    #: "velocity": uniform inlet velocity u0 (kept for comparisons only: it forces a plug profile into every pore opening)
    inlet: str = "pressure"
    p_in: float = 1.0
    u0: float = 1e-4
    mu: float = 1.0
    rho: float = 1.0
    #: MINRES settings. MINRES stops on a backward error of the preconditioned system, so the TRUE relative residual is
    #: always computed afterwards and the solution is rejected when it exceeds ``accept_residual``.
    minres_tol: float = 1e-7
    minres_maxiter: int = 20000
    accept_residual: float = 1e-3
    #: abort before assembly when the number of unknowns exceeds this value (protects shared machines)
    max_dofs: int = 6_000_000


@dataclass
class DescriptorOptions:
    """Thresholds of the hydraulic partition (descriptors.py)."""
    #: triangles with speed >= alpha * max speed are "active", <= beta * max speed are "dead" (defaults 0.01 and 0.001)
    alpha: float = 0.01
    beta: float = 0.001


@dataclass
class OpenFOAMOptions:
    """Finite-volume cross-check with OpenFOAM. On by default: every run is accompanied by an independent solution of the
    same domain with a code that shares neither discretization nor software with the finite-element workflow."""
    enabled: bool = True
    #: shell script that puts simpleFoam on PATH (OpenFOAM's etc/bashrc or a custom env file). None: auto-detect.
    env_script: Optional[str] = None
    #: kinematic viscosity of the FV run; large so that the run is in the creeping regime (Re << 1) at p_in = 1
    nu: float = 1.0e4
    #: images: the cross-check runs on the full domain when it has at most this many pore pixels, otherwise on
    #: square sub-windows of ``subwindow`` pixels (up to ``max_subwindows``) cut from the percolating window
    max_cells: int = 250_000
    subwindow: int = 512
    max_subwindows: int = 2
    #: MPI ranks for simpleFoam (1 = serial)
    nprocs: int = 1
    max_iterations: int = 200_000
    #: keep the OpenFOAM case directories (they can be large)
    keep_cases: bool = False


@dataclass
class OutputOptions:
    """Which files a run writes besides result.json and manifest.json (always written)."""
    save_mesh: bool = True
    save_fields: bool = True      # fields.npz (vertices, triangles, velocity, pressure)
    save_vtk: bool = True         # fields.vtu for ParaView
    figures: bool = True


@dataclass
class RunConfig:
    """The complete configuration of one run: one section per step of the workflow."""
    image: ImageOptions = field(default_factory=ImageOptions)
    domain: DomainOptions = field(default_factory=DomainOptions)
    mesh: MeshOptions = field(default_factory=MeshOptions)
    solver: SolverOptions = field(default_factory=SolverOptions)
    descriptors: DescriptorOptions = field(default_factory=DescriptorOptions)
    openfoam: OpenFOAMOptions = field(default_factory=OpenFOAMOptions)
    output: OutputOptions = field(default_factory=OutputOptions)

    # ------------------------------------------------------------------ (de)serialization
    def to_dict(self) -> dict:
        """Nested plain dict (written to manifest.json)."""
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict) -> "RunConfig":
        """Build a configuration from a nested dict; unknown keys raise an error so that typos cannot pass silently."""
        cfg = cls()
        for section, values in (data or {}).items():
            if not hasattr(cfg, section):
                raise KeyError(f"unknown configuration section '{section}'")
            target = getattr(cfg, section)
            names = {f.name for f in fields(target)}
            for key, value in (values or {}).items():
                if key not in names:
                    raise KeyError(f"unknown option '{section}.{key}'")
                setattr(target, key, value)
        cfg.validate()
        return cfg

    @classmethod
    def load(cls, path) -> "RunConfig":
        """Read a YAML (.yaml/.yml) or JSON file with the same nesting as ``to_dict``."""
        path = Path(path)
        text = path.read_text()
        if path.suffix.lower() in (".yaml", ".yml"):
            import yaml
            data = yaml.safe_load(text)
        else:
            data = json.loads(text)
        return cls.from_dict(data)

    def validate(self) -> None:
        """Reject inconsistent options before any computation starts (raises ValueError)."""
        checks = [
            (self.domain.direction in ("lr", "tb"), "domain.direction must be 'lr' or 'tb'"),
            (self.domain.keep_percolating in ("all", "largest"), "domain.keep_percolating must be 'all' or 'largest'"),
            (0.0 <= self.domain.reservoir_fraction <= 1.0, "domain.reservoir_fraction must lie in [0, 1]"),
            (self.image.pore_rule in ("green", "threshold", "nonzero"), "image.pore_rule must be green, threshold or nonzero"),
            (self.solver.element in ("p2", "mini"), "solver.element must be 'p2' or 'mini'"),
            (self.solver.method in ("lu", "minres"), "solver.method must be 'lu' or 'minres'"),
            (self.solver.inlet in ("pressure", "velocity"), "solver.inlet must be 'pressure' or 'velocity'"),
            (self.mesh.contour_tolerance >= 0, "mesh.contour_tolerance must be >= 0"),
            (not (self.solver.inlet == "pressure" and self.domain.reservoir_fraction == 0),
             "the pressure-driven definition needs reservoirs (domain.reservoir_fraction > 0)"),
        ]
        for ok, msg in checks:
            if not ok:
                raise ValueError(msg)
        self.mesh.sizes()


def dataclass_to_dict(obj):
    """Dict of a dataclass instance; any other object is returned unchanged."""
    return asdict(obj) if is_dataclass(obj) else obj
