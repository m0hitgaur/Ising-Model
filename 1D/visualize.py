import matplotlib.pyplot as plt
import numpy as np
import glob
from itertools import product
import os
import matplotlib.animation as animation
import pandas as pd
import seaborn as sns

from matplotlib.ticker import (MultipleLocator, AutoMinorLocator)

current_directory = os.getcwd()

plt.rcParams['axes.linewidth'] = 2.0
plt.rcParams["legend.labelspacing"]=0.1
#plt.rc_context({"xtick.major.pad": 8})
#plt.rc_context({"ytick.major.pad": 5})
plt.rc('xtick',labelsize=15)
plt.rc('ytick',labelsize=15)
#plt.rc('legend',fontsize=30)
plt.rcParams["font.family"] = "serif"
plt.rcParams['mathtext.fontset'] ="cm"
plt.rcParams['xtick.direction'] = 'in'
plt.rcParams['ytick.direction'] = 'in'
plt.rcParams['xtick.major.size'] = 10
#plt.rcParams['xtick.major.width'] = 2
plt.rcParams['xtick.minor.size'] = 5
#plt.rcParams['xtick.minor.width'] = 1.5
plt.rcParams['ytick.major.size'] = 10
#plt.rcParams['ytick.major.width'] = 2
plt.rcParams['ytick.minor.size'] = 5
#plt.rcParams['ytick.minor.width'] = 1.5
plt.rcParams["legend.handlelength"] = 1.0
plt.rcParams["legend.handletextpad"] = 0.2

sns.set_theme(style="white") # 'whitegrid' provides a grid, 'deep' is a good default color palette
sns.set_context("paper")



run=input("Run : ")
temp="0.10"


# Function to load data from files
def load_data_int(path):
    data = []
    with open(path, 'r') as f:
            for line in f:
                data.append(int(line.strip()))
    return np.array(data)

def load_data_float(path):
    data = []
    with open(path, 'r') as f:
            for line in f:
                data.append(line.strip().split())
    return np.array(data)

detail=[]
with open(current_directory+f'/details.txt', 'r') as file:
    for line in file:
        detail.append(line)
N = int(detail[0])   # Size of the grid
J=float(detail[1])
maxiter=int(detail[2])
total_runs=int(detail[3])

times=[]
for i in range( maxiter):
    if(i<10):tf=1
    if(i>10):tf=10
    if(i>100):tf=50
    if(i>1000):tf=100
    if(i%tf==0):times.append(i)

image=[]
for time in times: 
    path_lattice = current_directory+f'/data/Temp_{temp}/run_{run}/lattice_{time}.dat'
    lattice = load_data_int(path_lattice)
    image.append(lattice)
image=np.array(image)
fig,axes=plt.subplots(2,2,figsize=(10,10))
data=pd.read_csv("ising_1d_averaged_results.csv")

T=data["T"]
Energy=data["Energy"]
Energy_err=data["E_err"]
Mag=data["Mag_abs"]
Mag_err=data["M_err"]
C=data["HeatCap"]
C_err=data["C_err"]
Chi=data["Suscept"]
Chi_err=data["Chi_err"]
axes[0,0].imshow(image.T,origin="lower",aspect="auto")
axes[0,1].plot(T,Energy,label="Energy")
axes[0,1].plot(T,Mag,label="$|Mag|$")
axes[1,0].plot(T,Chi,label=r"$\chi$")
axes[1,1].plot(T,C,label=r"$C$")
axes[1,1].set_xlabel(r"$T$",size=25)
axes[1,0].set_xlabel(r"$T$",size=25)
axes[1,1].set_ylabel(r"$C$",size=25)
axes[1,0].set_ylabel(r"$\chi$",size=25)
axes[0,1].legend()

# Show the animation
plt.show()
