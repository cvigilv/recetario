#!/usr/bin/env python3
# title           :calculate_shape.py
# description     :Calculate shape Tanimoto/Tversky similarity matrix
# author          :Carlos Vigil Vásquez
# date            :20220830
# version         :20220830a
# notes           :Requires rdkit & numpy
# copyright       :Copyright (C) 2022 Carlos Vigil Vásquez (cvigil2@uc.cl).
# license         :Permission to copy and modify is granted under the MIT license

import argparse
import itertools
import multiprocessing
import numpy as np
import pandas as pd
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import rdMolAlign, rdShapeHelpers, rdMolDescriptors, AllChem

__version__ = "20220830a"
__date__ = "2022.08.30"
__title__ = "calculate_shape.py"
__desc__ = "Calculate shape Tanimoto/Tversky similarity matrix"


"""
Based on https://iwatobipen.wordpress.com/2018/09/26/3d-alignment-function-of-rdkit-rdkit/
"""
def MMFP_ShapeTversky(ref, query, ref_mmfp, query_mmfp, conf_ids, alpha, beta):
    ref_confid, query_confid = conf_ids
    # Align molecules
    pyO3A = rdMolAlign.GetO3A(
        ref,
        query,
        ref_mmfp,
        query_mmfp,
        ref_confid,
        query_confid,
    )
    pyO3A.Align()

    # Calculate Shape Tanimoto
    similarity = rdShapeHelpers.ShapeTverskyIndex(
        ref, query, alpha, beta, ref_confid, query_confid
    )

    return similarity


"""
Based on https://iwatobipen.wordpress.com/2018/09/26/3d-alignment-function-of-rdkit-rdkit/
"""


def Crippen_ShapeTversky(ref, query, ref_params, query_params, conf_ids, alpha, beta):
    ref_confid, query_confid = conf_ids

    # Align molecules
    pyO3A = rdMolAlign.GetCrippenO3A(
        ref, query, ref_params, query_params, ref_confid, query_confid
    )
    pyO3A.Align()

    # Calculate shape similarity
    similarity = rdShapeHelpers.ShapeTverskyIndex(
        ref, query, alpha, beta, ref_confid, query_confid
    )

    return similarity


def main():
    print(f"{__title__} v{__version__} - {__desc__}\n")
    # Argument parser {{{
    parser = argparse.ArgumentParser(description="Calculate shape similarity matrix")
    parser.add_argument("-i", "--input", type=str, nargs="+", help="SMILES files")
    parser.add_argument(
        "-o", "--output", type=str, default="shapesimmat.txt", help="Similarity matrix"
    )
    parser.add_argument("--delim", type=str, default=" ", help="SMILES files delimiter")
    parser.add_argument("--alpha", type=float, default=1.0, help="Tversky α parameter")
    parser.add_argument("--beta", type=float, default=1.0, help="Tversky β parameter")
    parser.add_argument(
        "--method",
        type=str,
        default="crippen",
        help="Alignment method (options: ['crippen', 'mmfp'])",
    )
    parser.add_argument(
        "--nconfs", type=int, default=25, help="Number of conformers per compounds"
    )
    parser.add_argument(
        "--prunerms", type=float, default=0.5, help="Conformer prune RMS"
    )
    parser.add_argument(
        "--nthreads",
        type=int,
        default=1,
        help="Number of threads to use ('-1' for all threads)",
    )

    args = parser.parse_args()

    if 0 < len(args.input) < 3:
        pass
    else:
        print(f"[FATAL] Passed {len(args.i)} SMILES files, a maximum of 2 are allowed!")
        exit(1)

    if args.method not in ["crippen", "mmfp"]:
        print(
            f"[FATAL] Method {args.method} not available! Try with either 'crippen' or 'mmfp'"
        )
        exit(1)

    if args.nthreads == -1:
        args.nthreads = multiprocessing.cpu_count()

    print(":: Configuration ::")
    print("===================")
    for item, value in vars(args).items():
        if type(value) is list:
            value = "; ".join(value)
        print(f"{item.ljust(10)} = {value}")
    print()
    # }}}

    # Conformer generation
    print(":: Conformer generation ::")
    print("==========================")
    ref_mols = [
        Chem.AddHs(m)
        for m in Chem.SmilesMolSupplier(args.input[0], delimiter=args.delim, titleLine=False)
        if m != None
    ]
    for mol in tqdm(
        ref_mols,
        desc=f"Calculating {args.nconfs} conformers for reference molecules",
        smoothing=0.0,
    ):
        AllChem.EmbedMultipleConfs(
            mol,
            clearConfs=True,
            numConfs=args.nconfs,
            numThreads=args.nthreads,
            pruneRmsThresh=args.prunerms,
        )

    if len(args.input) == 1:
        query_mols = ref_mols
    else:
        query_mols = [
            Chem.AddHs(m)
            for m in Chem.SmilesMolSupplier(args.input[1], delimiter=args.delim, titleLine=False)
            if m != None
        ]
        for mol in tqdm(
            query_mols,
            desc=f"Calculating {args.nconfs} conformers for query molecules    ",
            smoothing=0.0,
        ):
            AllChem.EmbedMultipleConfs(
                mol,
                clearConfs=True,
                numConfs=args.nconfs,
                numThreads=args.nthreads,
                pruneRmsThresh=args.prunerms,
            )

    conformer_distribution = [mol.GetNumConformers() for mol in ref_mols] + [
        mol.GetNumConformers() for mol in query_mols
    ]
    print("Descriptive statistics:")
    print(f"{' - Number of reference molecules'.ljust(40)} = {len(ref_mols)}")
    if len(args.input) == 2:
        print(f"{' - Number of query molecules'.ljust(40)} = {len(query_mols)}")
    print("Number of conformers per molecule:")
    print(
        f"{' - Mean ± SD'.ljust(40)} = {np.nanmean(conformer_distribution):.2f} ± {np.nanstd(conformer_distribution):.2f}"
    )
    print(
        f"{' - Median (IQR)'.ljust(40)} = {np.nanmedian(conformer_distribution):.0f} ({np.nanpercentile(conformer_distribution, 25):.0f} - {np.nanpercentile(conformer_distribution, 75):.0f})"
    )

    # Alignment and similarity calculation
    print()
    print(":: Shape similarity calculation ::")
    print("==================================")

    # Initialization
    shape_func = None
    ref_params = []
    query_params = []

    # Parameter calculation
    if args.method == "mmfp":
        shape_func = MMFP_ShapeTversky
        ref_params = [AllChem.MMFFGetMoleculeProperties(mol) for mol in ref_mols]
        query_params = [AllChem.MMFFGetMoleculeProperties(mol) for mol in query_mols]
    elif args.method == "crippen":
        shape_func = Crippen_ShapeTversky
        ref_params = [rdMolDescriptors._CalcCrippenContribs(mol) for mol in ref_mols]
        query_params = [
            rdMolDescriptors._CalcCrippenContribs(mol) for mol in query_mols
        ]

    # Similarity matrix calculation
    simmat = np.zeros((len(ref_mols), len(query_mols)))
    pool = multiprocessing.Pool(args.nthreads)

    if len(args.input) == 1 and (args.alpha == args.beta):
        pbar = tqdm(
            total=int(len(ref_mols) * (len(ref_mols) - 1) / 2),
            desc="Calculating shape similarity",
            smoothing=0.0,
        )
        for i, ref in enumerate(ref_mols):
            for j, query in enumerate(ref_mols[0:i]):
                if ref == query:
                    all_tanimotos = [1]
                else:
                    # Initialize ShapeTanimoto
                    all_tanimotos = []

                    # Parallelized Shape Tanimoto for every possible pair of conformers
                    # between reference and query molecule
                    for similarity in pool.starmap(
                        shape_func,
                        zip(
                            itertools.repeat(ref),
                            itertools.repeat(query),
                            itertools.repeat(ref_params[i]),
                            itertools.repeat(query_params[j]),
                            itertools.product(
                                range(ref.GetNumConformers()),
                                range(query.GetNumConformers()),
                            ),
                            itertools.repeat(args.alpha),
                            itertools.repeat(args.beta),
                        ),
                    ):
                        all_tanimotos.append(similarity)

                # Assign max(ShapeTanimoto) to molecule pair
                simmat[i, j] = np.max(all_tanimotos)
                pbar.update(1)

        simmat = np.maximum(simmat, simmat.transpose())

    else:
        pbar = tqdm(
            total=len(ref_mols) * len(query_mols),
            desc="Calculating shape similarity",
            smoothing=0.0,
        )
        for i, ref in enumerate(ref_mols):
            for j, query in enumerate(query_mols):
                if ref == query:
                    all_tanimotos = [1]
                else:
                    # Initialize ShapeTanimoto
                    all_tanimotos = []

                    # Parallelized Shape Tanimoto for every possible pair of conformers
                    # between reference and query molecule
                    for similarity in pool.starmap(
                        shape_func,
                        zip(
                            itertools.repeat(ref),
                            itertools.repeat(query),
                            itertools.repeat(ref_params[i]),
                            itertools.repeat(query_params[j]),
                            itertools.product(
                                range(ref.GetNumConformers()),
                                range(query.GetNumConformers()),
                            ),
                            itertools.repeat(args.alpha),
                            itertools.repeat(args.beta),
                        ),
                    ):
                        all_tanimotos.append(similarity)

                # Assign max(ShapeTanimoto) to molecule pair
                simmat[i, j] = np.max(all_tanimotos)
                pbar.update(1)
    pbar.close()

    print(
        f"{'Mean ± SD'.ljust(15)} = {np.nanmean(simmat.flatten()):.3f} ± {np.nanstd(simmat.flatten()):.3f}"
    )
    print(
        f"{'Median (IQR)'.ljust(15)} = {np.nanmedian(simmat.flatten()):.3f} ({np.nanpercentile(simmat.flatten(), 25):.3f} - {np.nanpercentile(simmat.flatten(), 75):.3f})"
    )

    simmat_df = pd.DataFrame(simmat)
    simmat_df.index = [mol.GetProp("_Name") for mol in ref_mols]
    simmat_df.columns = [mol.GetProp("_Name") for mol in query_mols]
    print(simmat_df)

    simmat_df.to_csv(args.output, header=True, index=True, sep=args.delim)



if __name__ == "__main__":
    main()
