import numpy as np
import matplotlib.pyplot as plt
import imageio.v2 as imageio # Use v2 for clearer intent
import os
import glob
import re # To extract numbers from filenames

# --- Configuration ---

# The base path where your 'data' folder is located
# Example: If your files are in /home/user/sims/data/lattice_...txt, set this to '/home/user/sims'
base_folder = '/home/mohit-gaur/Desktop/Angle_90/Ising_model/2D/'
numberoftrials=1
# Folder to save the generated GIFs
gif_output_folder = 'ising_gifs'

# Duration of each frame in the GIF (in milliseconds)
# 100 ms = 0.1 seconds = 10 frames per second
frame_duration_ms = 100

# Colormap for visualizing the lattice (e.g., 'coolwarm' for -1/1, 'binary' for 0/1)
# Adjust vmin/vmax if your data is not strictly -1 and 1
colormap = 'coolwarm'
vmin_val = -1
vmax_val = 1

# --- Process each trial ---
for trial in range(numberoftrials):
    print(f"\n--- Processing Trial {trial} ---")

    temp_frame_folder = os.path.join(gif_output_folder, f'temp_frames_trial_{trial}')
    os.makedirs(temp_frame_folder, exist_ok=True) # Create temporary folder for frames
    
    frame_filepaths = [] # To store paths of generated frames

plus_time=[]
minus_time=[]
time_plot=[]
for time in range(50000):
    tf=1
    if(time<10):tf=1
    if(time>10):tf=10
    if(time>100):tf=50
    if(time>1000):tf=100
    if(time>10000):tf=500
    if(time%tf==0):


        # Read the lattice data from the file
        lattice_data = np.loadtxt(base_folder+f'/data/lattice_{trial}_{time}_.txt')
        plus=0
        minus=0
        for row in lattice_data:
            for element in row:
                if(element==1):plus+=1
                else:minus+=1
        plus_time.append(plus)
        minus_time.append(minus)
        time_plot.append(time)

plt.plot(time_plot,plus_time)
plt.plot(time_plot,minus_time)
# Ensure tight layout to prevent labels/title overlapping
plt.tight_layout()
plt.show()
            
print("\n--- Script Finished ---")