#!/bin/sh
#title           :calculate_simcomp.sh
#description     :Calculate SIMCOMP similarity through API
#author          :Carlos Vigil Vásquez
#date            :20220511
#version         :20220511a
#notes           :
#copyright       :Copyright (C) 2022 Carlos Vigil Vásquez (cvigil2@uc.cl).
#license         :Permission to copy and modify is granted under the MIT license

SMILES=$1
SDF="/tmp/$(basename $SMILES).sdf"
SIMCOMP="$(echo $SMILES | sed s/SMILES/DD/ | sed s/smi/simcomp.txt/)"

echo "Calculating KEEG's SIMCOMP2 for $1"

# Convert SMILES to mol
obabel -ismi "$SMILES" -omol -O "$SDF" --gen2d

# Calculate SIMCOMP2 similarity measure
curl -F query_file=@"$SDF" \
	-F target_file=@"$SDF" \
	-F cutoff=0 \
	http://rest.genome.jp/simcomp2/ > "$SIMCOMP"
