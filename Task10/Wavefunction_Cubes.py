from gpaw import GPAW, FD, FermiDirac
from ase.optimize import GPMin
from ase.io import write, read

import numpy as np
import os as os

# Writing to Cuve files requires these files:
from ase.units import Bohr
from gpaw import restart

for n in [6, 7, 8]:
    filename = f'Na{n}'

    if not os.path.exists(filename):
        os.makedirs(filename)

    nanocluster = read(f'../Task5/{filename}.xyz')

    nanocluster.set_cell([16, 16, 16]) # 16 x 16 x 16 [Ångström] cell.
    nanocluster.center() # The cell must be centered in order to prevent atoms from lying too close to the boundary,
                             # as the boundary conditions are zero by default.

    calc = GPAW(nbands=10,
                h=0.25,
                txt=f'{filename}/out.txt',
                xc = 'LDA', # 'LDA' is the default.
                occupations=FermiDirac(0.05),
                setups={'Na': '1'},
                mode=FD(nn=3))
    nanocluster.calc = calc

    relax = GPMin(nanocluster, trajectory=f'{filename}/relax.traj', logfile=f'{filename}/relax.log') # Can also use 'BFGS', 'QuasiNewton', etc.
    relax.run(fmax=0.02, steps=100)

    E_Na6 = nanocluster.get_potential_energy()

    print(f'Final potential energy = {E_Na6} eV.')

    # Wrtie to '.gpw' file:
    calc.write(f'{filename}/{filename}.gpw', mode = 'all')

    # *** Second part *** #
    # Saving to cuve file:
    # Load binary file and get calculator
    atoms, calc = restart(f'{filename}/{filename}.gpw')

    # Loop over all wavefunctions (wf's) and write their cube files:
    nbands = calc.get_number_of_bands()
    for band in range(nbands):
        wf = calc.get_pseudo_wave_function(band=band)
        fname = f'{filename}/{filename}_{band}.cube'
        print('Writing wf', band, 'to file', fname)
        write(fname, atoms, data=wf * Bohr**1.5)