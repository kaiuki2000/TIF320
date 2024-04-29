import time
start_time = time.time()

from random import random
from ase.io import write
from ase.db import connect
from ase.optimize import GPMin
from gpaw import GPAW, FermiDirac

from ase.ga.data import DataConnection
from ase.ga.population import Population
from ase.ga.standard_comparators import InteratomicDistanceComparator
from ase.ga.cutandsplicepairing import CutAndSplicePairing
from ase.ga.utilities import closest_distances_generator
from ase.ga.utilities import get_all_atom_types
from ase.ga.offspring_creator import OperationSelector
from ase.ga.standardmutations import MirrorMutation
from ase.ga.standardmutations import RattleMutation
from ase.ga.standardmutations import PermutationMutation
import argparse

# Initialize the different components of the GA
da = DataConnection('gadb.db')

write('my_all_candidates.traj', da.get_all_relaxed_candidates())