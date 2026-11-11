"""High-level runs. Each run writes a self-describing directory:

    <out>/
      manifest.json             provenance: config, input/output hashes, versions, host, times (provenance.py)
      result.json               permeability, percolation report, descriptors, solver diagnostics, cross-check summary
      mask_selected.png         the analysed (percolating) domain                         [image runs]
      mesh/mesh.msh             the finite-element mesh (Gmsh 4.1)                        [if output.save_mesh]
      fields.npz / fields.vtu   vertex velocity and pressure                              [if output.save_fields / save_vtk]
      fields.png                pressure and speed maps                                   [if output.figures]
      openfoam/openfoam_crosscheck.json  finite-volume cross-check                        [if openfoam.enabled]

Use ``percoflow.provenance.verify_manifest(out)`` (CLI: ``percoflow verify <out>``) to check the files later.
"""
from __future__ import annotations

import json
import shutil
from pathlib import Path

import numpy as np

from . import descriptors, imaging, meshing, percolation, permeability, stokes
from .config import RunConfig
from .core import solve_window
from .provenance import Manifest


def _json(path, data):
    """Write ``data`` as indented JSON; NumPy scalars become Python numbers, anything else unknown becomes a string."""
    Path(path).write_text(json.dumps(data, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o)))


def _save_fields(sol, out: Path, cfg: RunConfig, shift=(0.0, 0.0)):
    """Write the vertex velocity and pressure as fields.npz (NumPy) and fields.vtu (ParaView), as configured."""
    ux, uy = sol.vertex_velocity()
    p = sol.vertex_pressure()
    pts = sol.mesh.p.copy()
    if cfg.output.save_fields:
        np.savez_compressed(out / "fields.npz", points=pts, triangles=sol.mesh.t, ux=ux, uy=uy, p=p,
                            note="vertex values of the Taylor-Hood/Mini solution; coordinates in the padded frame (px)")
    if cfg.output.save_vtk:
        import meshio
        meshio.write(out / "fields.vtu", meshio.Mesh(np.vstack([pts, np.zeros(pts.shape[1])]).T, [("triangle", sol.mesh.t.T)],
                                                     point_data={"velocity": np.column_stack([ux, uy, 0 * ux]), "pressure": p,
                                                                 "speed": np.hypot(ux, uy)}))


def _figure(sol, out: Path, title: str):
    """Quick-look figure fields.png: pressure (linear scale) and speed (log scale)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from matplotlib.tri import Triangulation
    ux, uy = sol.vertex_velocity()
    sp = np.hypot(ux, uy)
    tri = Triangulation(sol.mesh.p[0], sol.mesh.p[1], sol.mesh.t.T)
    fig, ax = plt.subplots(1, 2, figsize=(11, 5))
    a = ax[0].tripcolor(tri, sol.vertex_pressure(), shading="gouraud", cmap="viridis", rasterized=True)
    fig.colorbar(a, ax=ax[0], label="pressure p")
    pos = sp[sp > 0]
    b = ax[1].tripcolor(tri, np.clip(sp, pos.min() if pos.size else 1e-30, None), shading="gouraud", cmap="inferno",
                        norm=LogNorm(vmin=max(sp.max() * 1e-4, 1e-30), vmax=sp.max()), rasterized=True)
    fig.colorbar(b, ax=ax[1], label="speed |u|")
    for x, t in zip(ax, ("pressure", "speed (log)")):
        x.set_aspect("equal")
        x.set_title(t)
        x.set_xlabel("x [px]")
        x.set_ylabel("y [px]")
    fig.suptitle(title)
    fig.tight_layout()
    fig.savefig(out / "fields.png", dpi=150)
    plt.close(fig)


def run_image(image, out, cfg: RunConfig | None = None) -> dict:
    """Image -> pore mask -> window -> percolating domain -> mesh -> Stokes -> k_D -> descriptors (+ OpenFOAM cross-check)."""
    cfg = cfg or RunConfig()
    cfg.validate()
    out = Path(out)
    with Manifest(out, "image", cfg.to_dict(), inputs=[image]) as man:      # manifest.json is written even on failure
        mask = imaging.clean_mask(imaging.load_mask(image, cfg.image), cfg.image)          # image -> pore mask
        window, crop = imaging.crop_window(mask, cfg.domain)                               # remove the two bands
        selected, perc = percolation.analyse(window, cfg.domain.direction, cfg.domain.keep_percolating)
        result = {"input": str(image), "mask": imaging.describe_mask(mask), "crop": crop, "percolation": perc,
                  "direction": cfg.domain.direction, "units": "px (lengths), px^2 (permeability)"}
        from PIL import Image
        Image.fromarray((selected * 255).astype(np.uint8)).save(out / "mask_selected.png")
        if not perc["percolates"]:
            result["status"] = "not_percolating"
            man.note("no cluster joins the inlet and outlet planes: no flow problem was formulated")
            _json(out / "result.json", result)
            return result
        result["critical_throat_width_px"] = percolation.critical_throat_width(selected, cfg.domain.direction)
        fe = solve_window(selected, cfg.domain.direction, out / "mesh", cfg.mesh, cfg.solver, cfg.domain.reservoir_fraction)
        sol = fe["solution"]
        result.update(status="ok", permeability=fe["darcy"], solver=sol.info, mesh=fe["mesh"], timings=fe["timings"])
        if cfg.solver.inlet == "pressure":
            result["descriptors"] = descriptors.compute(sol, cfg.domain.direction, fe["darcy"], fe["b"], fe["L_s"],
                                                        cfg.descriptors.alpha, cfg.descriptors.beta)
        _save_fields(sol, out, cfg)
        if cfg.output.figures:
            _figure(sol, out, f"{Path(image).name}  k_D = {fe['darcy']['k_D']:.5g} px^2")
        if not cfg.output.save_mesh:
            shutil.rmtree(out / "mesh", ignore_errors=True)
        result["openfoam_crosscheck"] = _openfoam_image(selected, cfg, out, man)
        _json(out / "result.json", result)
    return result


def _openfoam_image(selected, cfg: RunConfig, out: Path, man: Manifest) -> dict:
    """Run (or record the explicit disabling of) the OpenFOAM cross-check of an image run."""
    if not cfg.openfoam.enabled:
        man.note("OpenFOAM cross-check DISABLED by the user (--no-openfoam)")
        return {"status": "disabled_by_user"}
    from .benchmarks.openfoam.pixel import crosscheck_image
    rep = crosscheck_image(selected, cfg.domain.direction, out / "openfoam", cfg.mesh, cfg.solver, cfg.openfoam,
                           cfg.domain.reservoir_fraction)
    man.data["openfoam"] = rep.get("openfoam")
    return {k: v for k, v in rep.items() if k != "pairs"} | {"pairs": [{k: p[k] for k in ("tag", "k_fe_contour", "k_fv_pixel",
                                                                                         "rel_diff_fv_vs_fe", "flag")} for p in rep.get("pairs", [])]}


def run_mesh(mesh_path, out, cfg: RunConfig | None = None, length=None, width=None) -> dict:
    """User mesh (groups INLET, OUTLET, WALLS) -> Stokes -> k_D from boundary pressures (+ OpenFOAM prism cross-check)."""
    cfg = cfg or RunConfig()
    cfg.validate()
    out = Path(out)
    with Manifest(out, "mesh", cfg.to_dict(), inputs=[mesh_path]) as man:
        mesh = meshing.load_mesh(mesh_path)
        sol = stokes.solve(mesh, cfg.domain.direction, cfg.solver)
        darcy = permeability.darcy_from_boundaries(sol, cfg.domain.direction, length, width, cfg.solver.mu)
        result = {"input": str(mesh_path), "status": "ok", "mesh": meshing.mesh_summary(mesh), "permeability": darcy,
                  "solver": sol.info, "direction": cfg.domain.direction}
        result["descriptors"] = descriptors.compute(sol, cfg.domain.direction, darcy, 0, darcy["L_s"],
                                                    cfg.descriptors.alpha, cfg.descriptors.beta)
        _save_fields(sol, out, cfg)
        if cfg.output.figures:
            _figure(sol, out, f"{Path(mesh_path).name}  k = {darcy['k_D']:.5g}")
        if cfg.openfoam.enabled:
            from .benchmarks.openfoam.triangles import crosscheck_mesh
            rep = crosscheck_mesh(mesh, darcy, out / "openfoam", cfg.openfoam, cfg.solver.p_in)
            man.data["openfoam"] = rep.get("openfoam")
            result["openfoam_crosscheck"] = rep
        else:
            man.note("OpenFOAM cross-check DISABLED by the user (--no-openfoam)")
            result["openfoam_crosscheck"] = {"status": "disabled_by_user"}
        _json(out / "result.json", result)
    return result
