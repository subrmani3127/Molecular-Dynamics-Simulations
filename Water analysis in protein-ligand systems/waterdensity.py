import MDAnalysis as mda
from MDAnalysis.analysis.density import DensityAnalysis
import numpy as np

u = mda.Universe('PR_Protein_bgc_434_center.gro', 'PR_Protein_bgc_434_center.xtc')

print("Identifying ligand...")
lig_atoms = u.select_atoms('resname bgc')
print(f"Ligand atoms: {len(lig_atoms)}")

print("Finding SOL water molecules...")
all_waters = u.select_atoms('resname SOL and name OW')
print(f"SOL oxygens: {len(all_waters)}")

waters = u.select_atoms('resname SOL and name OW and around 15 resname bgc', updating=True)
print(f"Waters around bgc: {len(waters)}")

print("Precomputing ligand COM trajectory...")
lig_coms = []
for ts in u.trajectory:
    lig_coms.append(lig_atoms.center_of_mass())
lig_coms = np.array(lig_coms)
u.trajectory[0]  # reset to frame 0

print("Density Analysis started .....")

D = DensityAnalysis(waters,
                   delta=2.0).run()

print("✓ Density computed!")
D.export("water_density_434.dx")
print("✓ Saved water_density.dx")
