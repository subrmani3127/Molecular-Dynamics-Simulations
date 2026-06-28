# Activate environment in current directory
import Pkg;
Pkg.activate(".");

# Run this once, to install necessary packages:
# Pkg.add(["ComplexMixtures", "PDBTools", "Plots", "LaTeXStrings"])

# Load packages
using ComplexMixtures
using PDBTools
using Plots, Plots.Measures
using LaTeXStrings

# The complete trajectory file can be downloaded from (1Gb):
# https://www.dropbox.com/scl/fi/zfq4o21dkttobg2pqd41m/glyc50_traj.dcd?rlkey=el3k6t0fx6w5yiqktyx96gzg6&dl=0

# The example output file is available at:
# 
# Load PDB file of the system
atoms = readPDB("PR_protein_bgc_41_center.pdb")

# Select the protein and the GLYC molecules
protein = select(atoms, "protein")
glyc = select(atoms, "resname SOL")

# Setup solute and solvent structures
solute = AtomSelection(protein, nmols=1)
solvent = AtomSelection(glyc, natomspermol=3)

# Path to the trajectory file
trajectory_file = "PR_protein_bgc_41_center.xtc"

# Run mddf calculation, and save results
results = mddf(trajectory_file, solute, solvent, Options(bulk_range=(8.0, 12.0)))
save(results, "PR_PB_41_water_results.json")
println("Results saved to PR_41_results.json")
