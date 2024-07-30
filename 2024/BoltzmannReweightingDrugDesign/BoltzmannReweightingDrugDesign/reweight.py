import os
import json
import pandas as pd
import numpy as np
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--file', help='The path of table including smiles and docking score',required=True)
parser.add_argument('--N', help='We will sample N*Num_of_mol times according to Boltzmann Distribution',default=30)

args = parser.parse_args()

def softmax(x):
    e_x = np.exp(x)
    return e_x / np.sum(e_x, axis=-1, keepdims=True)

word = ["3H","2H","13C","+","14","9","18","O-","-","P","N","."," ","14C","[2H]","Ca", "Mg", "K", "Na","Co",  "Cu",  "Fe", "Al", "Mo", "Rb", "Mn", "Se", "Zr", "Ga", "Te",  "Sn", "Ti", "W", "Si", "Bi", "As", "B", "Be", "Ni", "Ge", "V", "Ag", "Hg", "Cd", "Sr", "Sb", "Au", "Ba", "Ta", "Pb", "Li",  "Y", "Ru", "Cs", "Pd", "Pt", "Ce", "La", "Nd", "Re", "Tl", "Gd", "Zn", "Hf", "Th", "He", "Ar", "Nb", "Sm","Rh", "Rf", "Cr", "Eu", "In", "Ir", "Tm", "Pr", "Ho", "Lu", "Tb", "Er", "Dy", "Os"]

def get_smile(cids,df):
    kekuleSmiles = True
    smile = []
    vina_s = []
    vina_m = []
    CID = []
    for cid in cids:
        #sdf_file_path = f"sdf/{cid}.sdf"
        mol = Chem.SDMolSupplier(sdf_file_path)[0]
        smi = Chem.MolToSmiles(mol, kekuleSmiles=kekuleSmiles)
        has_isotope = any(atom.GetIsotope() != 0 for atom in mol.GetAtoms())
        if any(element in smi for element in word) or len(smi)>49 and not has_isotope:
            continue
        smile.append(smi)
        vina_m.append(df.loc[df["ligand"]==cid]["M"].item())
        vina_s.append(df.loc[df["ligand"] == cid]["S"].item())
        CID.append(cid)
    return smile,vina_m,vina_s,CID

N = int(args.N)  #copies
df = pd.read_excel(args.file)  #input


smile = list(df.smile)
cid = list(df.CID)
vina_s = list(df.score)







T=0.6
p = softmax([-x/T for x in vina_s])
re_cid = np.random.choice(cid, p=p,size = N*len(smile),replace=True)
index = [cid.index(i) for i in re_cid]
re_vina_s=[vina_s[i] for i in index]
re_smile = [smile[i] for i in index]
out_re = re_smile
out_re.insert(0,"smiles")
with open('reweighting.csv', 'w') as file:
    for item in out_re:
        file.write('%s\n' % item)




