#!/usr/bin/env python3
# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "numpy",
#     "pandas",
#     "rdkit",
# ]
# ///

import sys

from numpy import array
from pandas import DataFrame
from rdkit.Chem import SmilesMolSupplier
from rdkit.Chem.AllChem import GetMorganFingerprintAsBitVect

SMILES = sys.argv[1]
FCFP4 = sys.argv[2]

print(f"Calculating RDKit's FCFP4 fingerprint for {SMILES}")

# Load SMILES
all_smiles = SmilesMolSupplier(SMILES, delimiter="\t\t", titleLine=False)

# Calculate FCFP4 and clean-up data
args = dict(useChirality=True, radius=2, nBits=2048, useFeatures=True)
bits = [GetMorganFingerprintAsBitVect(molecule, **args) for molecule in all_smiles]
bits = [array(fingerprint) for fingerprint in bits]


# Save fingerprint
DataFrame(bits).to_csv(FCFP4, sep=" ", header=False, index=False)
