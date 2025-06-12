#!/usr/bin/env -S uv run --script
# title           :getmurcko
# description     :Calculate Murcko scaffolds for molecules in SMILES file
# author          :Carlos Vigil Vásquez
# date            :20220912
# version         :20250430a
# notes           :Requires pandas, tqdm & rdkit. Run with `uv run getmurcko.py`
# copyright       :Copyright (C) 2025 Carlos Vigil Vásquez (carlos.vigil.v@gmail.com).
# license         :Permission to copy and modify is granted under the MIT license
# /// script
# requires-python = ">=3.13"
# dependencies = [
#     "rdkit",
#     "tqdm",
# ]
# ///

import argparse
from datetime import datetime
from rdkit import Chem
from rdkit.Chem.Scaffolds.MurckoScaffold import GetScaffoldForMol, MakeScaffoldGeneric

__version__ = "20250430a"
__date__ = datetime.now().strftime("%Y-%m-%d")
__title__ = "getmurcko.py"
__desc__ = "Calculate Murcko scaffolds from SMILES"

def get_murcko_scaffold(mol, generic=False):
    """
    Calculate Murcko molecular scaffold from molecule.
    """
    scaffold = GetScaffoldForMol(mol)

    if generic:
        scaffold = MakeScaffoldGeneric(scaffold)

    scaffold_smile = Chem.CanonSmiles(Chem.MolToSmiles(scaffold))

    return scaffold_smile, mol.GetProp("_Name")

def main():
    print(f"{__title__} v{__version__} - {__desc__}\n")
    # Argument parser
    parser = argparse.ArgumentParser(description="Calculate Murcko scaffolds")
    parser.add_argument("-i", "--input", type=str, nargs="+", help="SMILES files", required=True)
    parser.add_argument("-o", "--output", type=str, default="", help="Scaffolds SMILES file")
    parser.add_argument("--delim", type=str, default=" ", help="SMILES files delimiter")
    parser.add_argument("--generic", action="store_true", help="Generic Murcko scaffolds")
    parser.set_defaults(generic = False)
    args = parser.parse_args()

    all_mols = [
        m
        for m in Chem.SmilesMolSupplier(args.input[0], delimiter=args.delim, titleLine=False)
        if m != None
    ]

    all_scaffolds = [
        get_murcko_scaffold(m, generic = args.generic)
        for m in all_mols
    ]

    with open(args.output, "w+") as io:
        for scaffold, molid in all_scaffolds:
            io.write(f"{scaffold} {molid}\n")

    print(f"Number of compounds: {len(all_mols)}")
    print(f"Number of unique scaffolds: {len(set([scf for scf,_ in all_scaffolds]))}")

if __name__ == "__main__":
    main()
