#!/usr/bin/env python3
"""Generate the two small example inputs of the tutorial (README, "Quick start"):

    examples/synthetic_rock.png   200 x 240 px segmented image, pores in pure green (the default pore rule),
                                  solid discs on a fluid background, with an isolated pore on purpose (removed by the
                                  percolation filter)
    examples/obstacles.msh        Gmsh triangle mesh of a 2D channel with three circular obstacles, physical groups
                                  FLUID, INLET (x = 0), OUTLET (x = 4), WALLS (channel walls and obstacles)

Both are SYNTHETIC test inputs, not data. Deterministic (fixed seed). Run:  python examples/make_examples.py
"""
from pathlib import Path

import gmsh
import numpy as np
from PIL import Image

HERE = Path(__file__).resolve().parent


def synthetic_rock(path, ny=200, nx=240, seed=3):
    """Random solid discs (radius 6-16 px) until 40 % of the image is solid, plus one isolated pore on purpose."""
    rng = np.random.default_rng(seed)
    yy, xx = np.mgrid[0:ny, 0:nx]
    solid = np.zeros((ny, nx), bool)
    while solid.mean() < 0.40:
        cy, cx, r = rng.uniform(0, ny), rng.uniform(0, nx), rng.uniform(6, 16)
        solid |= (yy - cy) ** 2 + (xx - cx) ** 2 <= r * r
    # an isolated pore inside a solid block (the percolation filter must remove it)
    solid[90:130, 100:140] = True
    solid[100:120, 110:130] = False
    rgb = np.zeros((ny, nx, 3), np.uint8)
    rgb[~solid] = (0, 255, 0)            # pores: pure green
    rgb[solid] = (40, 40, 40)            # solid: dark grey
    Image.fromarray(rgb).save(path)


def obstacles_mesh(path, h=0.04):
    """4 x 1 channel minus three discs, meshed with Gmsh; boundary curves sorted into INLET, OUTLET and WALLS."""
    gmsh.initialize()
    gmsh.option.setNumber("General.Terminal", 0)
    try:
        occ = gmsh.model.occ
        gmsh.model.add("obstacles")
        box = occ.addRectangle(0, 0, 0, 4, 1)
        discs = [occ.addDisk(x, y, 0, r, r) for x, y, r in ((1.0, 0.45, 0.2), (2.0, 0.6, 0.25), (3.0, 0.4, 0.18))]
        fluid, _ = occ.cut([(2, box)], [(2, d) for d in discs])
        occ.synchronize()
        surf = [t for _, t in fluid]
        inlet, outlet, walls = [], [], []
        for _, c in gmsh.model.getBoundary(fluid, oriented=False):
            x0, y0, _, x1, y1, _ = gmsh.model.getBoundingBox(1, c)
            if x1 < 1e-6:            # OCC bounding boxes carry a small tolerance
                inlet.append(c)
            elif x0 > 4 - 1e-6:
                outlet.append(c)
            else:
                walls.append(c)
        gmsh.model.addPhysicalGroup(2, surf, name="FLUID")
        gmsh.model.addPhysicalGroup(1, inlet, name="INLET")
        gmsh.model.addPhysicalGroup(1, outlet, name="OUTLET")
        gmsh.model.addPhysicalGroup(1, walls, name="WALLS")
        gmsh.option.setNumber("Mesh.MeshSizeMax", h)
        gmsh.option.setNumber("Mesh.MshFileVersion", 4.1)
        gmsh.model.mesh.generate(2)
        gmsh.write(str(path))
    finally:
        gmsh.finalize()


if __name__ == "__main__":
    synthetic_rock(HERE / "synthetic_rock.png")
    obstacles_mesh(HERE / "obstacles.msh")
    print("wrote", HERE / "synthetic_rock.png", "and", HERE / "obstacles.msh")
