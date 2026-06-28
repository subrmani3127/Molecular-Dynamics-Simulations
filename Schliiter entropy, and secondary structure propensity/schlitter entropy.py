import numpy as np
import MDAnalysis as mda
from numpy.linalg import eigvalsh

def calculate_schlitter_entropy(cov_matrix, temperature=300):
    k_B = 1.380649e-23  # Boltzmann constant (J/K)
    hbar = 1.0545718e-34  # Reduced Planck constant (J·s)
    beta = k_B * temperature / (hbar ** 2)
    I = np.identity(cov_matrix.shape[0])
    mat = I + beta * cov_matrix
    eigenvalues = eigvalsh(mat)
    eigenvalues = eigenvalues[eigenvalues > 0]
    entropy = 0.5 * np.sum(np.log(eigenvalues)) * k_B
    return entropy

# Load your structure and trajectory files
# Replace 'protein.pdb' and 'traj.xtc' with your actual file names
u = mda.Universe('O08_c.gro', 'O08_control-sys.xtc')
atom_selection = u.select_atoms('protein and name CA')  # Example: C-alpha atoms

entropies = []
frame_indices = []

# Collect coordinates for all frames skipping 50 frames interval to computes the coordinates faster
coords_all = []
for i, ts in enumerate(u.trajectory[:50]):
    coords_all.append(atom_selection.positions.flatten())
    frame_indices.append(i+1)
coords_all = np.array(coords_all)

# Calculate the covariance matrix over all frames (as in the second image)
mean_coords = np.mean(coords_all, axis=0)
centered = coords_all - mean_coords
cov = np.cov(centered, rowvar=False)

# Calculate Schlitter entropy for the whole trajectory
S = calculate_schlitter_entropy(cov)
NA = 6.02214076e23  # Avogadro's number
S_mol = S * NA  # J/(mol·K)

# For a cumulative curve, calculate entropy for increasing numbers of frames
with open('schlitter_entropy_output.txt', 'w') as f:
    f.write('Frame\tSchlitter_Entropy_J_per_mol_K\n')
    for i in range(5, len(coords_all)+1):
        cov_cum = np.cov(centered[:i], rowvar=False)
        S_cum = calculate_schlitter_entropy(cov_cum)
        S_cum_mol = S_cum * NA
        f.write(f'{i}\t{S_cum_mol:.4f}\n')
