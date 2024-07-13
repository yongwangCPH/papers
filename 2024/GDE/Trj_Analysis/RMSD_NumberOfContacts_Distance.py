#!/usr/bin/env python
# coding: utf-8

######For one trajectory######
import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from MDAnalysis.analysis import rms
from MDAnalysis.lib.distances import distance_array
import pandas as pd
from matplotlib.ticker import MultipleLocator

# Load the simulation files
u = mda.Universe('../md_center.gro', '../md_conti_center.xtc')

# Define the ligand and protein selections
ligand = u.select_atoms('resname AGLC')  # Adjust 'AGLC' to the actual ligand residue name
protein = u.select_atoms('protein')
selected_group = u.select_atoms('resid 173 406 461 557 487')  # Adjust 'resid 173 406 461 557 487' to your selected group

# Check if selections are empty
if len(ligand) == 0:
    raise ValueError("Ligand selection is empty. Please check the residue name.")
if len(selected_group) == 0:
    raise ValueError("Selected group is empty. Please check the residue IDs.")

# Initialize arrays to store the results
rmsd_values = []
com_distances = []
contact_counts = []

# RMSD calculation object
R = rms.RMSD(u, ligand, select='resname AGLC')  # Adjust selection if needed
R.run()

# Iterate through each frame in the trajectory
for ts in u.trajectory:
    # Calculate the RMSD for the current frame
    rmsd_values.append(R.rmsd[ts.frame, 2])  # RMSD values are stored in the 3rd column

    # Calculate the center of mass distance between the ligand and the selected group
    com_distance = np.linalg.norm(ligand.center_of_mass() - selected_group.center_of_mass())
    com_distances.append(com_distance)

    # Calculate the number of contacts within 5 Å between the ligand and the protein
    contact_matrix = distance_array(ligand.positions, protein.positions)
    contact_count = np.sum(contact_matrix < 5.0)
    contact_counts.append(contact_count)

# Create a dataframe for seaborn
data = pd.DataFrame({
    'COM Distance (Å)': com_distances,
    'RMSD of Ligand (Å)': rmsd_values,
    'Number of Contacts': contact_counts
})

# Customize plot settings
x_range = (2, 12)  # Change as needed
y_range = (0, 4)   # Change as needed
x_major_ticks = np.arange(4, 12, 2)  # Change as needed
x_minor_ticks = np.arange(4, 12, 2)  # Change as needed
y_major_ticks = np.arange(0.5, 4.5, 1)   # Change as needed
y_minor_ticks = np.arange(0.5, 4.5, 1) # Change as needed
v_min = 0  # Change as needed
v_max = 1000  # Change as needed

# Plot the results
plt.figure(figsize=(12, 8))
sns.set(style="white")

# Create a JointGrid
g = sns.JointGrid(data=data, x="COM Distance (Å)", y="RMSD of Ligand (Å)")

# Add a kde plot with color mapping to the density
kde_plot = g.plot_joint(sns.kdeplot, cmap="YlOrBr", fill=True, cbar=False)
g.plot_marginals(sns.histplot, kde=True, color="y")

# Set axis limits
g.ax_joint.set_xlim(x_range)
g.ax_joint.set_ylim(y_range)

# Set major and minor ticks
g.ax_joint.set_xticks(x_major_ticks)
g.ax_joint.set_xticks(x_minor_ticks, minor=True)
g.ax_joint.set_yticks(y_major_ticks)
g.ax_joint.set_yticks(y_minor_ticks, minor=True)

# Add grid
g.ax_joint.grid(which='both')

# Add a separate colorbar for the KDE plot
cbar_ax = kde_plot.figure.add_axes([.93, .25, .02, .4])
norm = plt.Normalize(v_min, v_max)
sm = plt.cm.ScalarMappable(cmap="YlOrBr", norm=norm)
sm.set_array(data['Number of Contacts'])
kde_plot.figure.colorbar(sm, cax=cbar_ax, label='Number of Contacts')
plt.subplots_adjust(right=0.9)  # Adjust the plot to make space for the color bar
plt.savefig('../R1_RMSD_COM_Contacts.png', dpi=300, bbox_inches='tight')
plt.show()

######For averaging three trajectories######
import MDAnalysis as mda
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from MDAnalysis.analysis import rms
from MDAnalysis.lib.distances import distance_array
import pandas as pd
from matplotlib.ticker import MultipleLocator

# Function to analyze a single trajectory
def analyze_trajectory(gro_file, xtc_file):
    u = mda.Universe(gro_file, xtc_file)
    ligand = u.select_atoms('resname AGLC')  # Adjust 'AGLC' to the actual ligand residue name
    protein = u.select_atoms('protein')
    selected_group = u.select_atoms('resid 173 406 461 557 487')  # Adjust 'resid 173 406 461 557 487' to your selected group

    if len(ligand) == 0:
        raise ValueError("Ligand selection is empty. Please check the residue name.")
    if len(selected_group) == 0:
        raise ValueError("Selected group is empty. Please check the residue IDs.")

    rmsd_values = []
    com_distances = []
    contact_counts = []

    R = rms.RMSD(u, ligand, select='resname AGLC')  # Adjust selection if needed
    R.run()

    for ts in u.trajectory:
        rmsd_values.append(R.rmsd[ts.frame, 2])  # RMSD values are stored in the 3rd column
        com_distance = np.linalg.norm(ligand.center_of_mass() - selected_group.center_of_mass())
        com_distances.append(com_distance)
        contact_matrix = distance_array(ligand.positions, protein.positions)
        contact_count = np.sum(contact_matrix < 5.0)
        contact_counts.append(contact_count)

    return rmsd_values, com_distances, contact_counts

# Load the simulation files for three trajectories
gro_files = [
    '../R1/md_center.gro',
    '../R2/md_center.gro',
    '../R3/md_center.gro'
]

xtc_files = [
    '../R1/md_conti_center.xtc',
    '../R2/md_conti_center.xtc',
    '../R3/md_conti_center.xtc'
]

# Analyze each trajectory
results = [analyze_trajectory(gro_file, xtc_file) for gro_file, xtc_file in zip(gro_files, xtc_files)]

# Find the minimum length of the trajectories
min_length = min(len(res[0]) for res in results)

# Truncate all results to the minimum length
truncated_results = [(rmsd[:min_length], com[:min_length], contact[:min_length]) for rmsd, com, contact in results]

# Stack and average the results
stacked_rmsd_values = np.vstack([res[0] for res in truncated_results])
stacked_com_distances = np.vstack([res[1] for res in truncated_results])
stacked_contact_counts = np.vstack([res[2] for res in truncated_results])

avg_rmsd_values = np.mean(stacked_rmsd_values, axis=0)
avg_com_distances = np.mean(stacked_com_distances, axis=0)
avg_contact_counts = np.mean(stacked_contact_counts, axis=0)

# Create a dataframe for seaborn
data = pd.DataFrame({
    'COM Distance (Å)': avg_com_distances,
    'RMSD of Ligand (Å)': avg_rmsd_values,
    'Number of Contacts': avg_contact_counts
})

# Customize plot settings
x_range = (4, 9)  # Change as needed
y_range = (0, 2.8)   # Change as needed
x_major_ticks = np.arange(4.5, 9, 2)  # Change as needed
x_minor_ticks = np.arange(4.5, 9, 2)  # Change as needed
y_major_ticks = np.arange(0.4, 2.5, 1)   # Change as needed
y_minor_ticks = np.arange(0.4, 2.5, 1) # Change as needed
v_min = 0  # Change as needed
v_max = 1000  # Change as needed

# Plot the results
plt.figure(figsize=(12, 8))
sns.set(style="white")

# Create a JointGrid
g = sns.JointGrid(data=data, x="COM Distance (Å)", y="RMSD of Ligand (Å)")

# Add a kde plot with color mapping to the density
kde_plot = g.plot_joint(sns.kdeplot, cmap="YlOrBr", fill=True, cbar=False)
g.plot_marginals(sns.histplot, kde=True, color="y")

# Set axis limits
g.ax_joint.set_xlim(x_range)
g.ax_joint.set_ylim(y_range)

# Set major and minor ticks
g.ax_joint.set_xticks(x_major_ticks)
g.ax_joint.set_xticks(x_minor_ticks, minor=True)
g.ax_joint.set_yticks(y_major_ticks)
g.ax_joint.set_yticks(y_minor_ticks, minor=True)

# Add grid
g.ax_joint.grid(which='both')

# Add a separate colorbar for the KDE plot
cbar_ax = kde_plot.figure.add_axes([.93, .25, .02, .4])
norm = plt.Normalize(v_min, v_max)
sm = plt.cm.ScalarMappable(cmap="YlOrBr", norm=norm)
sm.set_array(data['Number of Contacts'])
kde_plot.figure.colorbar(sm, cax=cbar_ax, label='Number of Contacts')
plt.subplots_adjust(right=0.9)  # Adjust the plot to make space for the color bar
plt.savefig('../average_RMSD_COM_Contacts.png', dpi=300, bbox_inches='tight')
plt.show()

