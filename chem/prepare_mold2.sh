#!/bin/sh
#title           :prepare_mold2.sh
#description     :Calculate Mold2 molecular descriptor
#author          :Carlos Vigil Vásquez
#date            :20220511
#version         :20220511a
#notes           :Requires Mold2 binary
#copyright       :Copyright (C) 2022 Carlos Vigil Vásquez (cvigil2@uc.cl).
#license         :Permission to copy and modify is granted under the MIT license

SMILES=$1
SDF="/tmp/$(basename $SMILES).sdf"
MOLD2="$(echo $SMILES | sed s/SMILES/DS/ | sed s/smi/mold2.txt/)"

# Convert SMILES to SDF
obabel -ismi "$SMILES" -O "$SDF" --gen2d

# Calculate and clean-up descriptor
echo $SDF
Mold2 -i "$SDF" -o "/tmp/$(basename $MOLD2)" -r "/tmp/$(basename $MOLD2).report.txt"
cat "/tmp/$(basename $MOLD2)" | tail -n+2 | cut -d"	" -f2- > "$MOLD2"
sed -i "s/\t/\ /g" "$MOLD2"
