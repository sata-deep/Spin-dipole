#!/usr/bin/env python3
"""
Compute the atomic spin-dipole term Tz for a non-collinear (SOC) DFT calculation.

Expected inputs in the working directory
    Mx.cube  My.cube  Mz.cube   – spin-density components (μB / Å³)
    POSCAR                     – atomic positions (Å)  [ASE readable]
Adapt RWIGS{} to your PAW radii.
"""
import sys
import numpy as np
from pymatgen.io.common import VolumetricData          # Cube loader
from ase.io import read                                # POSCAR reader

# ----------------------------------------------------------------------
# 0.  Filenames
# ----------------------------------------------------------------------
cube_files = {"x": "Mx.cube", "y": "My.cube", "z": "Mz.cube"}
poscar     = "POSCAR"

# ----------------------------------------------------------------------
# 1.  Read the three spin-density cubes  -------------------------------
#     (All share the same grid and lattice.)
# ----------------------------------------------------------------------
mx = VolumetricData.from_cube(cube_files["x"])
my = VolumetricData.from_cube(cube_files["y"])
mz = VolumetricData.from_cube(cube_files["z"])

for c in (my, mz):
    assert c.dim == mx.dim,  "CUBE grids differ!"
    assert np.allclose(c.structure.lattice.matrix,
                       mx.structure.lattice.matrix), "Lattices differ!"

m_grid = np.stack([mx.data["total"],
                   my.data["total"],
                   mz.data["total"]], axis=-1)       # (Nx,Ny,Nz,3)

lattice = mx.structure.lattice.matrix               # 3×3 in Å
nx, ny, nz = mx.dim
axes = np.vstack((lattice[0] / nx,
                  lattice[1] / ny,
                  lattice[2] / nz))                 # voxel vectors
vol_per_voxel = abs(np.linalg.det(axes))            # Å³

# ----------------------------------------------------------------------
# 2.  Atomic structure (sphere centres)  ------------------------------
# ----------------------------------------------------------------------
atoms     = read(poscar)
symbols   = atoms.get_chemical_symbols()
positions = atoms.get_positions()                   # Å

# This is an example for a system with Pd and Co only. Please make your directory
# according to your POSCAR file.

RWIGS = {                                           # radii in Å
    "Pd": 2.710,
    "Co": 2.100,
}
default_rwigs = 2.0

# ----------------------------------------------------------------------
# 3.  Pre-compute voxel Cartesian coordinates  ------------------------
# ----------------------------------------------------------------------
ix, iy, iz = np.indices(mx.dim)
r_cart = (ix[..., None] * axes[0] +
          iy[..., None] * axes[1] +
          iz[..., None] * axes[2])                  # (Nx,Ny,Nz,3)

# ----------------------------------------------------------------------
# 4.  Integrate Tz inside each PAW sphere  ----------------------------
# ----------------------------------------------------------------------
tz_cell = 0.0
np.seterr(divide="ignore", invalid="ignore")        # silent 0/0 inside nuclei

for idx, (sym, pos) in enumerate(zip(symbols, positions), start=1):
    dvec = r_cart - pos                             # r − R_a
    r2   = np.einsum("...i,...i->...", dvec, dvec)  # |dvec|²
    mask = r2 < (RWIGS.get(sym, default_rwigs))**2

    # unit radial vector  r̂
    rhat = np.zeros_like(dvec)
    np.divide(dvec, np.sqrt(r2)[..., None], out=rhat, where=r2[..., None] > 0)

    # r̂ · m
    r_dot_m = np.einsum("...i,...i->...", rhat, m_grid)

    # integrand ½ [m_z − 3 r̂_z (r̂·m)]
    tz_vox = 0.5 * (m_grid[..., 2] - 3.0 * rhat[..., 2] * r_dot_m)

    tz = np.sum(tz_vox[mask]) * vol_per_voxel
    print(f"{idx:3d} {sym:2s}   Tz = {tz: .5e} μB")
    tz_cell += tz

print(f"\nCell-total Tz  = {tz_cell: .5e} μB")
print(f"Cell-total 7 Tz = {7*tz_cell: .5e} μB")
