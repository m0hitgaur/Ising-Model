import matplotlib.pyplot as plt
import numpy as np
import glob
from itertools import product
import os
import matplotlib.animation as animation

current_directory = os.getcwd()

run=input("Run : ")



# Function to load data from files
def load_data(path):
    data = []
    with open(path, 'r') as f:
            for line in f:
                data.append(int(line.strip()))
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
    path_lattice = current_directory+f'/data/run_{run}/lattice_{time}.dat'
    lattice = load_data(path_lattice)
    image.append(lattice)

plt.imshow(image)

# Show the animation
plt.show()
