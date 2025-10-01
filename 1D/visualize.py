import matplotlib.pyplot as plt
import numpy as np
import glob
from itertools import product
import os
import matplotlib.animation as animation
import pandas as pd
current_directory = os.getcwd()

run=input("Run : ")
temp=input("Temp : ")


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
fig,axes=plt.subplots(1,2,figsize=(10,10))
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
axes[0].imshow(image.T,origin="lower",aspect="auto")
axes[1].plot(T,Energy)

# Show the animation
plt.show()
