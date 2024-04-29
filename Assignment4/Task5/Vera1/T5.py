from ase.build import molecule
from gpaw import GPAW, PW
from ase.vibrations import Vibrations
from ase.thermochemistry import IdealGasThermo

# Definitions:
a = 12
h = 0.2
ecut = 450

for name in ['O2', 'CO']:
    sigma = 1 # Symmetry number - default = 1
    molecule_Obj = molecule(f'{name}')
    molecule_Obj.set_cell((a, a, a))
    molecule_Obj.center()
    calc = GPAW(h = h, mode = PW(ecut), xc = 'PBE', txt = f'T5_{name}.txt', symmetry = {'point_group': False})
    if(name == 'O2'): calc.set(hund=True); sigma = 2 # Spin-polarized calculations!
    molecule_Obj.calc = calc

    # Set up vibration calculations using previously optimized molecule
    vib = Vibrations(molecule_Obj, name = f'{name}_vib')
    vib.run()
    print(vib.summary())
    vib_energies = vib.get_energies()

    # Save vibration modes to review if imaginary frequency present
    for mode in range(3*len(molecule_Obj.numbers)):
      vib.write_mode(mode)

    potentialenergy = molecule_Obj.get_potential_energy()

    # Perform ideal gas calculations assuming ideal gas
    thermo = IdealGasThermo(vib_energies = vib_energies,
                            potentialenergy = potentialenergy,
                            atoms=molecule_Obj,
                            geometry='linear',       # adjust this
                            symmetrynumber=sigma,    # adjust this
                            spin=1                   # adjust this
                            )
    # Symmetry number of CO?
    # Calculate respective corrections at temperature and pressure specified
    G = thermo.get_entropy(temperature=300, pressure=100_000)
    print(f'[{name}]: Entropy = {G}.') # Units?
    # O2's value does seem to match the experimental value. Convert eV to kJ/mol to see this.
    # CO's also seems to be reasonably close.