#!/usr/bin/env python3
"""
MSM + TPT Analysis - FULLY CORRECTED
"""

import deeptime as dt
import mdtraj as md
import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
from deeptime.decomposition import TICA
from deeptime.clustering import KMeans
from deeptime.markov.msm import MaximumLikelihoodMSM
import pandas as pd

# ============================================================================
# CONFIGURATION
# ============================================================================
topology_file = 'PR_Protein_bgc_434_center.gro'
trajectory_file = 'PR_Protein_bgc_434_center.xtc'

CHUNK_SIZE = 1000
STRIDE = 10
TICA_LAG = 10
TICA_DIM = 2
N_CLUSTERS = 100
N_METASTABLE = 5
MSM_LAG = 10

u = mda.Universe(topology_file, trajectory_file)

print("="*60)
print("MSM + TPT ANALYSIS - CORRECTED")
print("="*60)

# ============================================================================
# 1. LOAD & EXTRACT FEATURES
# ============================================================================
print("\n1. Extracting distance features...")

cat_atoms = u.select_atoms("resname GLU and (resid 153 or resid 320) and name CA")
glucose_atoms = u.select_atoms("resname bgc")

print(f"   CAT atoms: {len(cat_atoms)}")
print(f"   Glucose atoms: {len(glucose_atoms)}")

features_list = []
for ts in u.trajectory[::STRIDE]:
    cat_coords = cat_atoms.positions
    glucose_coords = glucose_atoms.positions
    
    dists = mda.lib.distances.distance_array(
        cat_coords,
        glucose_coords,
        box=ts.dimensions
    )
    
    features_list.append(dists.ravel())

features = np.stack(features_list)
print(f"   Features shape: {features.shape}")
print(f"   Distance range: {features.min():.2f} - {features.max():.2f} Å")

# ============================================================================
# 2. TICA
# ============================================================================
print("\n2. Running TICA...")

tica_est = TICA(lagtime=TICA_LAG, dim=TICA_DIM)
tica_model = tica_est.fit(features).fetch_model()
tica_output = tica_model.transform(features)

print(f"   TICA output: {tica_output.shape}")

# ============================================================================
# 3. CLUSTERING
# ============================================================================
print("\n3. Clustering...")

cluster_est = KMeans(n_clusters=N_CLUSTERS, max_iter=100)
cluster_model = cluster_est.fit(tica_output).fetch_model()
dtraj = cluster_model.transform(tica_output)

print(f"   Discrete trajectory: {dtraj.shape}")

# ============================================================================
# 4. IMPLIED TIMESCALES
# ============================================================================
print("\n4. Computing implied timescales...")

lag_times = [10, 20, 50, 100, 200]
its_data = []

for lag in lag_times:
    try:
        msm_temp = MaximumLikelihoodMSM(lagtime=lag, reversible=True)
        msm_temp_model = msm_temp.fit(dtraj).fetch_model()
        
        # Get actual MSM from collection
        if hasattr(msm_temp_model, 'msms'):
            actual_msm = msm_temp_model[0]
        else:
            actual_msm = msm_temp_model
            
        its = actual_msm.timescales()[:3]
        its_data.append(its)
        print(f"   Lag {lag}: {its[:3]}")
    except Exception as e:
        print(f"   Lag {lag}: Failed - {e}")
        its_data.append(np.full(3, np.nan))

its_data = np.array(its_data)

# ============================================================================
# 5. BUILD FINAL MSM
# ============================================================================
print("\n5. Building MSM...")

msm_est = MaximumLikelihoodMSM(lagtime=MSM_LAG, reversible=True)
msm_collection = msm_est.fit(dtraj).fetch_model()

# CRITICAL: Extract actual MSM from collection
if hasattr(msm_collection, 'msms'):
    msm_model = msm_collection[0]  # Largest connected component
    print(f"   Using largest MSM from collection")
else:
    msm_model = msm_collection

print(f"   Active states: {msm_model.n_states}")
print(f"   Lag time: {msm_model.lagtime}")

print("\n   Top 5 Timescales:")
for i, ts in enumerate(msm_model.timescales()[:5]):
    print(f"     {i+1}. {ts:.2f} frames")

# ============================================================================
# 6. METASTABLE STATES (PCCA)
# ============================================================================
print("\n6. PCCA+ Metastable States...")

pcca = msm_model.pcca(N_METASTABLE)

active_set = msm_model.count_model.state_symbols
dtraj_active = dtraj[np.isin(dtraj, active_set)]
state_to_meta = {state: pcca.assignments[i] for i, state in enumerate(active_set)}
metastable_traj = np.array([state_to_meta[s] for s in dtraj_active])

print(f"\n   Metastable state populations:")
for i in range(N_METASTABLE):
    pop = np.sum(metastable_traj == i) / len(metastable_traj)
    print(f"     State {i}: {pop*100:.2f}%")

T_coarse = pcca.coarse_grained_transition_matrix
pi_coarse = pcca.coarse_grained_stationary_probability

# Save PCCA results
np.save('metastable_trajectory.npy', metastable_traj)
np.save('metastable_assignments.npy', pcca.assignments)
np.save('transition_matrix_coarse.npy', T_coarse)
np.save('metastable_memberships.npy', pcca.memberships)

print("   ✓ PCCA complete")

# ============================================================================
# 7. SAVE REPRESENTATIVE STRUCTURES
# ============================================================================
print("\n7. Extracting representative structures...")

for state in range(N_METASTABLE):
    state_frames = np.where(metastable_traj == state)[0]
    if len(state_frames) > 0:
        frame_idx = state_frames[0] * STRIDE
        try:
            frame = md.load_frame(trajectory_file, frame_idx, top=topology_file)
            frame.save_pdb(f'state_{state}.pdb')
            print(f"   State {state}: Saved")
        except:
            print(f"   State {state}: Failed to load frame")

# ============================================================================
# 8. TPT DIAGNOSTICS
# ============================================================================
print("\n" + "="*60)
print("8. TPT DIAGNOSTICS")
print("="*60)

def diagnose_msm(msm_model, pcca_assignments, source_meta=0, sink_meta=4):
    """Diagnose MSM connectivity."""
    
    source_states = np.where(pcca_assignments == source_meta)[0]
    sink_states = np.where(pcca_assignments == sink_meta)[0]
    
    print(f"\n   Diagnostics for {source_meta} → {sink_meta}:")
    print(f"   MSM states: {msm_model.n_states}")
    print(f"   Source microstates: {len(source_states)}")
    print(f"   Sink microstates: {len(sink_states)}")
    
    # Check connectivity
    source_to_sink = msm_model.transition_matrix[np.ix_(source_states, sink_states)].sum()
    
    # Populations
    source_pop = msm_model.stationary_distribution[source_states].sum()
    sink_pop = msm_model.stationary_distribution[sink_states].sum()
    
    print(f"   Source population: {source_pop:.3f}")
    print(f"   Sink population: {sink_pop:.3f}")
    print(f"   Direct flux: {source_to_sink:.6f}")
    
    if source_to_sink < 1e-6:
        print("   ⚠ WARNING: Very weak/no direct connection!")
        return False
    else:
        print("   ✓ States connected")
        return True

# Run diagnostics
diagnose_msm(msm_model, pcca.assignments, source_meta=0, sink_meta=4)

# ============================================================================
# 9. ALL-PAIRS TPT ANALYSIS
# ============================================================================
print("\n" + "="*60)
print("9. ALL-PAIRS TPT ANALYSIS")
print("="*60)
    
from deeptime.markov import ReactiveFlux

def all_pairs_tpt_analysis(msm_model, pcca_assignments, stride=10, frame_ps=10.0):
    """Compute TPT flux + MFPT for all pairs."""

    n_states = len(np.unique(pcca_assignments))
    print(f"\n   Analyzing {n_states} metastable states...")

    results = []

    for i in range(n_states):
        for j in range(n_states):
            if i == j:
                continue

            source_idx = np.where(pcca_assignments == i)[0]
            sink_idx = np.where(pcca_assignments == j)[0]

            if len(source_idx) < 10 or len(sink_idx) < 10:  # Min population
                continue

            try:
                # ✅ FIXED: Use np.array() not .tolist()
                flux = ReactiveFlux(
                    msm_model.transition_matrix,
                    source_idx,        # np.array ✅
                    sink_idx           # np.array ✅
                )

                total_flux = flux.total_flux

                # Simplified MFPT (using deeptime)
                mfpt_frames = msm_model.mfpt(source_idx, sink_idx)
                mfpt_ns = (mfpt_frames * stride * frame_ps) / 1000.0

                results.append({
                    'Source': i,
                    'Sink': j,
                    'Flux': total_flux,
                    'MFPT_ns': mfpt_ns,
                    'Source_pop': msm_model.stationary_distribution[source_idx].sum(),
                    'Sink_pop': msm_model.stationary_distribution[sink_idx].sum()
                })

            except Exception as e:
                print(f"   Skipped {i}→{j}: {str(e)[:40]}")
                continue

    df = pd.DataFrame(results)
    
    if len(df) == 0:
        print("   ⚠ No successful TPT calculations!")
        return df

    # Sort + Print
    df = df.sort_values('Flux', ascending=False).reset_index(drop=True)
    
    print(f"\n   📊 TOP PATHWAYS BY FLUX:")
    print(df[['Source', 'Sink', 'Flux', 'MFPT_ns']].round(4).to_string(index=False))
    
    df.to_csv('all_pairs_tpt.csv', index=False)
    print(f"\n   ✓ Saved: all_pairs_tpt.csv ({len(df)} pairs)")
    
    return df

# Run TPT analysis
flux_table = all_pairs_tpt_analysis(msm_model, pcca.assignments, stride=STRIDE)

if len(flux_table) > 0:
    print(f"\n   📈 NETWORK SUMMARY:")
    print(f"   Average flux: {flux_table['Flux'].mean():.6f}")
    print(f"   Strongest flux: {flux_table['Flux'].max():.6f}")
    print(f"   Connected pairs: {len(flux_table[flux_table['Flux'] > 1e-6])}")
    if not flux_table['MFPT_ns'].isna().all():
        finite_mfpt = flux_table['MFPT_ns'][np.isfinite(flux_table['MFPT_ns'])]
        if len(finite_mfpt) > 0:
            print(f"   Fastest transition: {finite_mfpt.min():.2f} ns")
            print(f"   Median transition: {finite_mfpt.median():.2f} ns")

# ============================================================================
# 10. VISUALIZATION
# ============================================================================
print("\n" + "="*60)
print("10. CREATING VISUALIZATIONS")
print("="*60)

fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Plot 1: Implied Timescales
ax = axes[0, 0]
for i in range(min(3, its_data.shape[1])):
    mask = ~np.isnan(its_data[:, i])
    if np.any(mask):
        ax.plot(np.array(lag_times)[mask],
                its_data[mask, i],
                'o-', label=f'ITS {i+1}')

ax.plot(lag_times, lag_times, 'k--', label='y=x')
ax.set_xlabel('Lag time (frames)')
ax.set_ylabel('Timescale (frames)')
ax.set_xscale('log')
ax.set_yscale('log')
ax.legend()
ax.grid(True, alpha=0.3)
ax.set_title('Implied Timescales')

# Plot 2: Free Energy Surface
ax = axes[0, 1]
h, xe, ye = np.histogram2d(tica_output[:, 0], tica_output[:, 1], bins=50)
fe = -np.log(h + 1e-10)
fe[fe == np.inf] = np.nan

im = ax.imshow(fe.T, origin='lower', aspect='auto',
               extent=[xe[0], xe[-1], ye[0], ye[-1]], cmap='viridis')
ax.set_xlabel("TICA 1")
ax.set_ylabel("TICA 2")
ax.set_title("Free Energy Surface")
plt.colorbar(im, ax=ax, label="Free Energy (kT)")

# Plot 3: Metastable States
ax = axes[1, 0]
full_meta_traj = np.full(len(dtraj), -1)
active_mask = np.isin(dtraj, active_set)
full_meta_traj[active_mask] = metastable_traj

scatter = ax.scatter(tica_output[:, 0], tica_output[:, 1],
                    c=full_meta_traj, cmap='tab10', s=1, alpha=0.5,
                    vmin=0, vmax=N_METASTABLE-1)
ax.set_xlabel("TICA 1")
ax.set_ylabel("TICA 2")
ax.set_title("Metastable States")
plt.colorbar(scatter, ax=ax, label="State")

# Plot 4: Transition Network
ax = axes[1, 1]
im = ax.imshow(T_coarse, cmap='Blues', vmin=0, vmax=1)
ax.set_xlabel('To State')
ax.set_ylabel('From State')
ax.set_title('Transition Matrix')
for i in range(N_METASTABLE):
    for j in range(N_METASTABLE):
        color = 'white' if T_coarse[i, j] < 0.5 else 'black'
        ax.text(j, i, f'{T_coarse[i, j]:.2f}',
               ha="center", va="center", fontsize=8, color=color)
plt.colorbar(im, ax=ax)

plt.tight_layout()
plt.savefig('msm_tpt_analysis.png', dpi=300, bbox_inches='tight')
print("   ✓ Saved: msm_tpt_analysis.png")

# ============================================================================
# 11. SAVE DATA
# ============================================================================
print("\n11. Saving all data...")

np.save('tica_output.npy', tica_output)
np.save('dtraj.npy', dtraj)
np.save('features.npy', features)

with open('msm_summary.txt', 'w') as f:
    f.write("MSM + TPT Analysis Summary\n")
    f.write("="*60 + "\n\n")
    f.write(f"Frames analyzed: {features.shape[0]}\n")
    f.write(f"Features: {features.shape[1]}\n")
    f.write(f"Clusters: {N_CLUSTERS}\n")
    f.write(f"Active states: {msm_model.n_states}\n")
    f.write(f"Metastable states: {N_METASTABLE}\n")
    f.write(f"MSM lag: {MSM_LAG}\n\n")
    f.write("Timescales (frames):\n")
    for i, ts in enumerate(msm_model.timescales()[:5]):
        f.write(f"  {i+1}. {ts:.2f}\n")
    f.write("\nMetastable populations:\n")
    for i in range(N_METASTABLE):
        pop = np.sum(metastable_traj == i) / len(metastable_traj)
        f.write(f"  State {i}: {pop*100:.2f}%\n")

print("   ✓ Saved: msm_summary.txt")

print("\n" + "="*60)
print("✅ ANALYSIS COMPLETE!")
print("="*60)
print("\n📁 Output files:")
print("   • msm_tpt_analysis.png")
print("   • all_pairs_tpt.csv")
print("   • msm_summary.txt")
print("   • state_*.pdb")
print("   • *.npy arrays")
