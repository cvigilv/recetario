#!/usr/bin/env python3
#title           :edgelist2matrix.py
#description     :Convert edge list to matrix
#author          :Carlos Vigil Vásquez
#date            :20220511
#version         :20220511a
#notes           :
#copyright       :Copyright (C) 2022 Carlos Vigil Vásquez (cvigil2@uc.cl).
#license         :Permission to copy and modify is granted under the MIT license

import sys
import pandas as pd

# Load edge list to dataframe
edgelist = pd.read_csv(
    sys.argv[1],
    sep="\t",
    header=None,
    index_col=None,
    names=["source", "target", "value"],
)
print(edgelist)
matrix = edgelist.pivot_table(index="source", columns="target", values="value")
matrix.fillna(0, inplace=True)
matrix = matrix.astype({c:int for c in matrix.columns})

matrix.to_csv(sys.argv[2], sep=" ")
