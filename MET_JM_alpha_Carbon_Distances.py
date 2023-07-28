# MET_JM_alpha_Carbon_Distances.py
# created by: Gabriella Estevam @ UCSF
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import colors
from matplotlib.ticker import PercentFormatter

#csv_array=pd.read_csv("DFGin_MET_aligned_strutures.csv", sep=',', header=0)
#csv_array=pd.read_csv("DFGin_MET_aligned_strutures.csv", sep=',', header=0)
csv_array=pd.read_csv("MET_ensemblle_JM_C_helix_resi_dist.csv", sep=',', header=0)

'''
dist_array0=csv_array["PDB"].to_numpy()
dist_array1=csv_array["1069 to 1121 distance"].to_numpy() #creats a numpy array - matplotlib readability
dist_array2=csv_array["1066 to 1129 distance"].to_numpy()
dist_array3=csv_array["1062 to 1126 distance"].to_numpy()
'''

dist_array0=csv_array["PDB"].to_numpy()
dist_array1=csv_array["1066 to 1129 distance"].to_numpy() #creats a numpy array - matplotlib readability
dist_array2=csv_array["1062 to 1125 distance"].to_numpy()
dist_array3=csv_array["1058 to 1121 distance"].to_numpy()

defined_bins1=np.arange(5,13,0.9)
defined_bins2=np.arange(5,13,0.9)
defined_bins3=np.arange(5,13,0.9)
#print defined_bins

# all aligned structures as histograms and scatters

fig,axs = plt.subplots(1,3)
fig.suptitle("JM and C-helix residue distances ")
axs[0].hist(dist_array1, bins=defined_bins1, color='steelblue',ec="black")
axs[0].set_title("V1066-I1129")
axs[0].set_ylim([0, 40])
axs[0].set_xlim([4, 12])
axs[0].set_ylabel("Count")
axs[0].set_xlabel("Alpha Carbon Distance ($\AA$)")

axs[1].hist(dist_array2, bins=defined_bins2,color='steelblue',ec="black")
axs[1].set_title("L1062-L1125")
axs[1].set_ylim([0, 40])
axs[1].set_xlim([4, 12])
axs[1].set_ylabel("Count")
axs[1].set_xlabel("Alpha Carbon Distance ($\AA$)")

axs[2].hist(dist_array3, bins=defined_bins3,color='steelblue',ec="black")
axs[2].set_title("L1058-V1121")
axs[2].set_ylim([0, 40])
axs[2].set_xlim([4, 12])
axs[2].set_ylabel("Count")
axs[2].set_xlabel("Alpha Carbon Distance ($\AA$)")
'''
axs[1,0].scatter(dist_array0, dist_array1)
axs[1,0].set_ylim([5, 15])

axs[1,1].scatter(dist_array0, dist_array2)
axs[1,1].set_ylim([5, 15])

axs[1,2].scatter(dist_array0, dist_array3)

axs[1,2].set_ylim([5, 15])
'''
plt.show()
