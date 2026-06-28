# I wrote this script to calculate the number of water molecules around the important catalytic residues of my protein of interest. To see the how much water is accessing through this residues of active site pocket.
import MDAnalysis as mda
import numpy as np

u = mda.Universe("PR_Protein_bgc_434_center.gro", "PR_Protein_bgc_434_center.xtc")

water_counts = []
# Define the radius (e.g., 5.0 Angstroms) around the catalytic residues
dist_cutoff = 5.0

print(f"Analyzing glucose molecules within {dist_cutoff} Å of active site...")

for ts in u.trajectory[::10]: # Using a stride of 10 for speed
    # Select water oxygens (OW) within the cutoff of residues 153 and 320
    # SOL or WAT is common for water resname; check your .gro file
    sel_string = f"resname bgc and around {dist_cutoff} (resid 110 or resid 318 or resid 323 or resid 156 or resid 160 or resid 166 or resid 167 or resid 168 or resid 208 or resid 230 or resid 233)"
    waters = u.select_atoms(sel_string)

    water_counts.append(len(waters))


mean_water = np.mean(water_counts)
print(f"Average glucose molecules within {dist_cutoff} Å of active site Pocket`: {mean_water}")
