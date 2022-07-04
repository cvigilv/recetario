#!/usr/bin/env python3
#title           :prepare_fcfp4.py
#description     :Calculate FCFP4 molecular descriptor
#author          :Carlos Vigil Vásquez
#date            :20220511
#version         :20220511a
#notes           :Requires rdkit, numpy & pandas
#copyright       :Copyright (C) 2022 Carlos Vigil Vásquez (cvigil2@uc.cl).
#license         :Permission to copy and modify is granted under the MIT license

import sys

from numpy import array
from pandas import DataFrame
from rdkit.Chem import SmilesMolSupplier
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect

SMILES = sys.argv[1]
FCFP4 = SMILES.replace("SMILES", "DS").replace("smi", "fcfp4.txt")

print(f"Calculating RDKit's FCFP4 fingerprint for {SMILES}")

# Load SMILES
all_smiles = SmilesMolSupplier(SMILES, delimiter="\t\t", titleLine=False)

# Calculate FCFP4 and clean-up data
args = dict(useChirality=True, radius=2, nBits=2048, useFeatures=True)
bits = [GetMorganFingerprintAsBitVect(molecule, **args) for molecule in all_smiles]
bits = [array(fingerprint) for fingerprint in bits]

# Save fingerprint
DataFrame(bits).to_csv(FCFP4, sep=" ", header=False, index=False)
