from ase import Atoms
from ase.visualize import view
from gpaw import GPAW, PW, FD, FermiDirac
from ase.optimize import GPMin

from ase.db import connect
from ase.io import write, read

import numpy as np
import os as os

for filename in ['Na6', 'Na6_2nd_Lowest_E']:

    if not os.path.exists(filename):
        os.makedirs(filename)

    nanocluster = read(f'../{filename}.xyz')

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

    calc.write(f'{filename}/{filename}.gpw', mode = 'all')