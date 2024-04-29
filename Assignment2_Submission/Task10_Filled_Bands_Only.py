from gpaw import GPAW, FD, FermiDirac
from ase.optimize import GPMin
from ase.io import write, read

import numpy as np
import os as os

# Writing to Cuve files requires these files:
from ase.units import Bohr
from gpaw import restart

for j in [6, 7, 8]:
    filename = f'Na{j}'

    # *** Second part *** #
    # Saving to cuve file:
    # Load binary file and get calculator
    atoms, calc = restart(f'{filename}/{filename}.gpw')

    # Loop over all wavefunctions (wf's) and write their cube files:
    nbands = calc.get_number_of_bands()
    for band in range(nbands):
        n = calc.get_occupation_numbers()[band]
        if(n > 0.5):
            print(f'[band {band}] Occupied band? Yes. Occupation = {n}.')
            wf = calc.get_pseudo_wave_function(band=band)
            fname = f'{filename}/{filename}_{band}.cube'
            print('Writing wf', band, 'to file', fname)
            write(fname, atoms, data=wf * Bohr**1.5)
        else:
            print(f'[band {band}] Occupied band? No. Occupation = {n}.')
    print(f'Finishing Na{j}...\n')

# This didn't work as well as I'd liked it to.
# filename = 'Na8'
# atoms, calc = restart(f'{filename}/{filename}.gpw')
# wf = np.array([calc.get_pseudo_wave_function(band=j) for j in range(4)])
# fname = f'{filename}/{filename}_Testing.cube'
# write(fname, atoms, data=wf * Bohr**1.5)