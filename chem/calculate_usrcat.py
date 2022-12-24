#!/usr/bin/env python3
# title           :calculate_usrcat.py
# description     :Calculate USRCAT similarity matrix
# author          :Carlos Vigil Vásquez
# date            :20220717
# version         :20220717a
# notes           :Requires rdkit, numpy & pandas
# copyright       :Copyright (C) 2022 Carlos Vigil Vásquez (cvigil2@uc.cl).
# license         :Permission to copy and modify is granted under the MIT license

import sys
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdBase
from rdkit.Chem import RDConfig
from rdkit.Chem import AllChem
from rdkit.Chem.rdMolDescriptors import GetUSRScore, GetUSRCAT

print(rdBase.rdkitVersion)
SMILES = sys.argv[1]
SMILES2 = sys.argv[2]
USRCAT = SMILES2.replace("SMILES", "DS").replace("smi", "usrcat.txt")

print(f"Calculating RDKit's USRCAT similarity for {SMILES}")

# Load SMILES and prepare molecules
mols = [
    Chem.AddHs(m)
    for m in Chem.SmilesMolSupplier(SMILES, delimiter="\t", titleLine=False)
    if m != None
]
for i, mol in enumerate(mols):
    if i % 50 == 0:
        print(i)
    AllChem.EmbedMolecule(
        mol, maxAttempts=100, useExpTorsionAnglePrefs=True, useBasicKnowledge=True
    )
mols2 = [
    Chem.AddHs(m)
    for m in Chem.SmilesMolSupplier(SMILES2, delimiter="\t", titleLine=False)
    if m != None
]
for i, mol in enumerate(mols2):
    if i % 50 == 0:
        print(i)
    AllChem.EmbedMolecule(
        mol, maxAttempts=100, useExpTorsionAnglePrefs=True, useBasicKnowledge=True
    )

# Calculate USRCAT
features = [GetUSRCAT(molecule) for molecule in mols]
features2 = [GetUSRCAT(molecule) for molecule in mols2]
similarities = []
for i in range(len(features2)):
    for j in range(len(features)):
        similarities.append(
            dict(
                source=mols2[i].GetProp("_Name"),
                target=mols[j].GetProp("_Name"),
                usrscore=GetUSRScore(features2[i], features[j]),
            )
        )
edgelist = pd.DataFrame(similarities)

# Save as matrix
sim_mat = edgelist.pivot_table(index="source", columns="target", values="usrscore")
sim_mat.to_csv(USRCAT, sep=" ", header=True, index=True)
