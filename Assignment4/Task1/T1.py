import numpy as np
import matplotlib.pyplot as plt

# 1st Part - Computing:
from ase import Atoms
from ase.io.trajectory import Trajectory
from gpaw import GPAW, PW

# 2nd Part - Plotting:
from ase.io import read
from ase.units import kJ
from ase.eos import EquationOfState


# Task 1:
# 1st Part - Computing:
a = 4.0   # Approximate lattice constant
b = a / 2 # FCC lattice

for Elem in ['Au', 'Pt', 'Rh']:
    calc = GPAW(xc   = 'PBE', 
            mode = PW(450),
            kpts = [12, 12, 12],
            txt  = f'T1_{Elem}.txt')
            
    Elem_Atoms = Atoms(f'{Elem}',
               cell = [(0, b, b), (b, 0, b), (b, b, 0)], # FCC lattice
               pbc  = 1, # Periodicity in all directions. I think that's what this means. Check latter. Bulk!
               calculator = calc)
    cell = Elem_Atoms.get_cell()
    traj = Trajectory(f'T1_{Elem}.traj', 'w')
    for x in np.linspace(0.80, 1.20, 40):
        Elem_Atoms.set_cell(cell * x, scale_atoms=True)
        Elem_Atoms.get_potential_energy()
        traj.write(Elem_Atoms)

    # 2nd Part - Plotting:
    configs = read(f'T1_{Elem}.traj@0:40')  # Read all 40 configurations.
    # Extract volumes and energies:
    volumes   = [Conf.get_volume() for Conf in configs]
    energies  = [Conf.get_potential_energy() for Conf in configs]
    eos       = EquationOfState(volumes, energies)
    v0, e0, B = eos.fit()
    print(B / kJ * 1.0e24, 'GPa')
    eos.plot(f'{Elem}-eos.png')
    plt.cla()