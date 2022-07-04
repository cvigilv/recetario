#!/bin/sh
#title           :prepare_ecfp4.sh
#description     :Calculate ECFP4 molecular descriptor
#author          :Carlos Vigil Vásquez
#date            :20220511
#version         :20220511a
#notes           :Requires obabel
#copyright       :Copyright (C) 2022 Carlos Vigil Vásquez (cvigil2@uc.cl).
#license         :Permission to copy and modify is granted under the MIT license

SMILES=$1
ECFP4="$(echo $SMILES | sed s/SMILES/DS/ | sed s/smi/ecfp4.txt/)"

# Calculate and clean-up descriptor
obabel -ismi "$SMILES" -ofpt -O "/tmp/$(basename $ECFP4)" -xh -xfECFP4

echo "Calculating OpenBabel's ECFP4 fingerprint for $1"
cat "/tmp/$(basename $ECFP4)" |\
	grep -v 'Possible superstructure of' |\
	grep -v '>' | paste -d '\0'  - - - - - - - - - - - - - - - - - - - - - -  |\
	sed s/' '//g | tr "[a-z]" "[A-Z]" | \
	sed s/0/0000/g | \
	sed s/1/0001/g | \
	sed s/2/0010/g | \
	sed s/3/0011/g | \
	sed s/4/0100/g | \
	sed s/5/0101/g | \
	sed s/6/0110/g | \
	sed s/7/0111/g | \
	sed s/8/1000/g | \
	sed s/9/1001/g | \
	sed s/A/1010/g | \
	sed s/B/1011/g | \
	sed s/C/1100/g | \
	sed s/D/1101/g | \
	sed s/E/1110/g | \
	sed s/F/1111/g | \
	sed 's/./& /g' | \
	sed 's/[[:space:]]*$//' > "$ECFP4"
