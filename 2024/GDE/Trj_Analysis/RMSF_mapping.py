#!/usr/bin/env python
# coding: utf-8

#######RMSF Mapping######
import pandas as pd

# Read the data files into DataFrames
data1 = pd.read_csv("../R1/rmsf_ProCa.xvg", comment="#", delim_whitespace=True, skiprows=17, names=["Residue", "RMSF_1"])
data2 = pd.read_csv("../R2/rmsf_ProCa.xvg", comment="#", delim_whitespace=True, skiprows=17, names=["Residue", "RMSF_2"])
data3 = pd.read_csv("../R3/rmsf_ProCa.xvg", comment="#", delim_whitespace=True, skiprows=17, names=["Residue", "RMSF_3"])
data4 = pd.read_csv("../R4/rmsf_ProCa.xvg", comment="#", delim_whitespace=True, skiprows=17, names=["Residue", "RMSF_4"])
data5 = pd.read_csv("../R5/rmsf_ProCa.xvg", comment="#", delim_whitespace=True, skiprows=17, names=["Residue", "RMSF_5"])
data6 = pd.read_csv("../R6/rmsf_ProCa.xvg", comment="#", delim_whitespace=True, skiprows=17, names=["Residue", "RMSF_6"])

# Merge the dataframes on 'Residue'
merged_data = pd.merge(data1, data2, on='Residue', suffixes=('_1', '_2'))
merged_data = pd.merge(merged_data, data3, on='Residue', suffixes=('', '_3'))
merged_data = pd.merge(merged_data, data4, on='Residue', suffixes=('', '_4'))
merged_data = pd.merge(merged_data, data5, on='Residue', suffixes=('', '_5'))
merged_data = pd.merge(merged_data, data6, on='Residue', suffixes=('', '_6'))

# Calculate the average RMSF
merged_data['Average RMSF'] = merged_data[['RMSF_1', 'RMSF_2', 'RMSF_3', 'RMSF_4', 'RMSF_5', 'RMSF_6']].mean(axis=1)

# Create a new dataframe with Residue and Average RMSF
new_dataframe = merged_data[['Residue', 'Average RMSF']]

# Print the new dataframe
print(new_dataframe)

# Uncomment to save to a CSV file
new_dataframe.to_csv("../average.csv", index=False)

from Bio.PDB import PDBParser, PDBIO
parser = PDBParser()
structure = parser.get_structure("GDE", "../GDE.pdb")

n = 1
column = 1
bfactor = []

# Read average RMSF values from the CSV file
with open('../average.csv', 'r') as f:
    next(f)  # Skip the header line
    for row in f:  # Iterate directly over the file object to avoid reading all lines into memory
        columns = row.split(',')
        if len(columns) > column:  # Ensure that there are enough columns
            bfactor.append(float(columns[column]))
          
bfactor

for residue, bf in zip(structure.get_residues(), bfactor):
    for atom in residue:
        atom.set_bfactor(bf)

io = PDBIO()
io.set_structure(structure)
io.save("../GDE_RMSF.pdb")  

