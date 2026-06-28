# This python script calculates the survival probability of the water moelcules around 4.0 angstrom from the protein center. I calculated this to measure against ligand suvival probability.

import MDAnalysis as mda
import matplotlib.pyplot as plt
from MDAnalysis.analysis.waterdynamics import SurvivalProbability as SP

# this variable consists of topology and trajectory file generated from simulations data.
pro = mda.Universe("PR_Protein_bgc_434_part1_center.gro","PR_Protein_bgc_434_part1_center.xtc")

# selecting this SOL mentioning water particularly
print(pro.select_atoms("resname SOL").atoms.names[10])

select = "byres name OW and around 4.0 protein"

sp = SP(pro,select,verbose=True)
sp.run(start=0,stop=100002,tau_max=20)

# computation of time series of tau and sp.
tau_timeseries = sp.tau_timeseries
sp_timeseries = sp.sp_timeseries

for tau,sp in zip(tau_timeseries,sp_timeseries):
    print("{time} {sp}".format(time=tau,sp=sp))


# Data Visualizations
plt.xlabel("Time")
plt.ylabel("SP")
plt.title("Survival Probability")
plt.plot(tau_timeseries,sp_timeseries)
plt.savefig("SP_PR_PB_434_part1_sol.png")
plt.show()
