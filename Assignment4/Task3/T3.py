from gpaw import GPAW, PW, Davidson, FermiDirac
from gpaw import restart
from gpaw.mpi import size, rank
from ase.build import fcc111
from ase.build import bulk
from ase.io import Trajectory
from ase.optimize import BFGS
import os

Elems     = ['Au', 'Pt', 'Rh']
Lattice_a = [4.18462, 3.97949, 3.85641]

for Elem, a, in zip(Elems, Lattice_a):
    print(f'Element - {Elem}.')
    # Set up the bulk copper crystal structure
    atoms = bulk(f'{Elem}', 'fcc', a = a)

    # Set up the surface slab
    slab = fcc111(f'{Elem}', size = (3, 3, 3), a = a, vacuum = 6.0)

    # Set up the calculator
    ecut = 450.0  # eV
    kpts = (4, 4, 1)
    calc = GPAW(mode = PW(ecut), xc = 'PBE', kpts = kpts, txt = f'T3_{Elem}.txt')

    # Set up the calculator for bulk
    atoms.set_calculator(calc)

    # Relax the surface
    slab.set_calculator(calc)
    dyn = BFGS(slab, trajectory = f'{Elem}_relax.traj')
    dyn.run(fmax = 0.05)

    # Calculate the surface energy
    energy_bulk    = atoms.get_potential_energy() / len(atoms)
    energy_slab    = slab.get_potential_energy()
    natoms_slab    = len(slab)
    natoms_bulk    = len(atoms)
    surface_area   = 2 * (slab.cell[0, 0] * slab.cell[1, 1])
    surface_energy = (energy_slab - natoms_slab / natoms_bulk * energy_bulk) / surface_area

    print("Surface energy: %f J/m^2" % surface_energy)
