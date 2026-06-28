import MDAnalysis as mda
import matplotlib.pyplot as plt
from MDAnalysis.analysis.waterdynamics import WaterOrientationalRelaxation as WOR

u = mda.Universe("PR_Protein_bgc_434_center.gro","PR_Protein_bgc_434_center.xtc")

select = "byres name OH2 and sphzone 6.0 protein"

WOR_analysis = WOR(u,select,0,100002,20)
WOR_analysis.run()
time =0


for WOR_OH, WOR_HH, WOR_dip in WOR_analysis.timeseries:
    print("{time} {WOR_OH} {time} {WOR_HH} {time} { WOR_dip}".format(time=time,WOR_OH = WOR_OH, WOR_HH = WOR_HH, WOR_dip = WOR_dip))
    time+=1

plt.figure(1,figsize=(18, 6))

#WOR OH
plt.subplot(131)
plt.xlabel('time')
plt.ylabel('WOR')
plt.title('WOR OH')
plt.plot(range(0,time),[column[0] for column in WOR_analysis.timeseries])

#WOR HH
plt.subplot(132)
plt.xlabel('time')
plt.ylabel('WOR')
plt.title('WOR HH')
plt.plot(range(0,time),[column[1] for column in WOR_analysis.timeseries])

#WOR dip
plt.subplot(133)
plt.xlabel('time')
plt.ylabel('WOR')
plt.title('WOR HH')
plt.plot(range(0,time),[column[1] for column in WOR_analysis.timeseries])

#WOR dip
plt.subplot(133)
plt.xlabel('time')
plt.ylabel('WOR')
plt.title('WOR dip')
plt.plot(range(0,time),[column[2] for column in WOR_analysis.timeseries])

plt.show()


