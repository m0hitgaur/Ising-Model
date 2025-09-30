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

            # Create a matplotlib figure and axes
            fig, ax = plt.subplots(figsize=(10,10)) # Adjust figure size as needed

            # Visualize the 2D lattice data
            # Use imshow with colormap and set vmin/vmax for consistent color scaling
            img = ax.imshow(lattice_data, cmap=colormap, vmin=vmin_val, vmax=vmax_val)

            # Set title and remove axis ticks/labels for a cleaner look
            ax.set_title(f'Trial: {trial}, Time: {time}')
            ax.set_xticks([])
            ax.set_yticks([])

            # Ensure tight layout to prevent labels/title overlapping
            plt.tight_layout()

            # Save the plot as an image file
            # dpi controls the resolution (dots per inch)
            plt.savefig(base_folder+gif_output_folder+f'/lattice_{trial}_{time}_.png', dpi=100)

            # Close the figure to free up memory
            plt.close(fig)

            
print("\n--- Script Finished ---")