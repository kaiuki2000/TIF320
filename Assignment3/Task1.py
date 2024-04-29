from ase.io import read
from gpaw import GPAW
from ase.md.npt import NPT
from ase.units import fs, kB
from ase.io import Trajectory
from datetime import datetime

startTime = datetime.now()

atoms = read('Task1_Snapshot_12000.xyz')

calc = GPAW(mode     = 'lcao',
            xc       = 'PBE',
            basis    = 'dzp',
            symmetry = {'point_group': False}, # Turn off point-group symmetry
            # setups={'Na': '1'},              # Do we need this?
            charge   = 1,                      # Remember, we have an Na+ ion!
            txt      = 'out.txt')              # Redirects calculator output to this file!
atoms.set_calculator(calc)

dyn = NPT(atoms,
          temperature_K = 350,
          timestep = 0.5*fs,               # We chose 0.5 fs as an adequate timestep for water. Check explanation in report.
          ttime = 20*fs,                   # Characteristic timescale for the thermostat.
          pfactor = None,                  # To disable the barostat, we include this.
          externalstress = 0,              # We don't use the barostat, but this needs to be set anyways!
          logfile = 'Task1_MD_Output.log') # Outputs temperature (and more) to file at each timestep.
trajectory = Trajectory('Task1_Dynamics.traj', 'w', atoms)

dyn.attach(trajectory.write, interval = 1) # Write the current positions, etc. to file each timestep.
dyn.run(4000)                              # Run 4000 steps of MD simulation. Corresponds to 2 ps, with our time-step.

# Printing run-time:
print(datetime.now() - startTime)