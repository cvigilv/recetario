#!/usr/bin/env python3
#title           :prepare_kr.sh
#description     :Calculate Klethora-Roth molecular descriptor
#author          :Carlos Vigil Vásquez
#date            :20220511
#version         :20220511a
#notes           :Requires padeldescriptor
#copyright       :Copyright (C) 2022 Carlos Vigil Vásquez (cvigil2@uc.cl).
#license         :Permission to copy and modify is granted under the MIT license

import sys
import os
import pandas as pd
from padelpy import padeldescriptor

config = """
<Root>
    <Group name="2D">
        <Descriptor name="AcidicGroupCount" value="true"/>
        <Descriptor name="ALOGP" value="true"/>
        <Descriptor name="AminoAcidCount" value="false"/>
        <Descriptor name="APol" value="true"/>
        <Descriptor name="AromaticAtomsCount" value="true"/>
        <Descriptor name="AromaticBondsCount" value="true"/>
        <Descriptor name="AtomCount" value="true"/>
        <Descriptor name="Autocorrelation" value="true"/>
        <Descriptor name="BaryszMatrix" value="true"/>
        <Descriptor name="BasicGroupCount" value="true"/>
        <Descriptor name="BCUT" value="true"/>
        <Descriptor name="BondCount" value="true"/>
        <Descriptor name="BPol" value="true"/>
        <Descriptor name="BurdenModifiedEigenvalues" value="true"/>
        <Descriptor name="CarbonTypes" value="true"/>
        <Descriptor name="ChiChain" value="true"/>
        <Descriptor name="ChiCluster" value="true"/>
        <Descriptor name="ChiPathCluster" value="true"/>
        <Descriptor name="ChiPath" value="true"/>
        <Descriptor name="Constitutional" value="true"/>
        <Descriptor name="Crippen" value="true"/>
        <Descriptor name="DetourMatrix" value="true"/>
        <Descriptor name="EccentricConnectivityIndex" value="true"/>
        <Descriptor name="EStateAtomType" value="true"/>
        <Descriptor name="ExtendedTopochemicalAtom" value="true"/>
        <Descriptor name="FMF" value="true"/>
        <Descriptor name="FragmentComplexity" value="true"/>
        <Descriptor name="HBondAcceptorCount" value="true"/>
        <Descriptor name="HBondDonorCount" value="true"/>
        <Descriptor name="HybridizationRatio" value="true"/>
        <Descriptor name="InformationContent" value="true"/>
        <Descriptor name="IPMolecularLearning" value="false"/>
        <Descriptor name="KappaShapeIndices" value="true"/>
        <Descriptor name="KierHallSmarts" value="false"/>
        <Descriptor name="LargestChain" value="true"/>
        <Descriptor name="LargestPiSystem" value="true"/>
        <Descriptor name="LongestAliphaticChain" value="true"/>
        <Descriptor name="MannholdLogP" value="true"/>
        <Descriptor name="McGowanVolume" value="true"/>
        <Descriptor name="MDE" value="true"/>
        <Descriptor name="MLFER" value="true"/>
        <Descriptor name="PathCount" value="true"/>
        <Descriptor name="PetitjeanNumber" value="true"/>
        <Descriptor name="RingCount" value="true"/>
        <Descriptor name="RotatableBondsCount" value="true"/>
        <Descriptor name="RuleOfFive" value="true"/>
        <Descriptor name="Topological" value="true"/>
        <Descriptor name="TopologicalCharge" value="true"/>
        <Descriptor name="TopologicalDistanceMatrix" value="true"/>
        <Descriptor name="TPSA" value="true"/>
        <Descriptor name="VABC" value="true"/>
        <Descriptor name="VAdjMa" value="true"/>
        <Descriptor name="WalkCount" value="true"/>
        <Descriptor name="Weight" value="true"/>
        <Descriptor name="WeightedPath" value="true"/>
        <Descriptor name="WienerNumbers" value="true"/>
        <Descriptor name="XLogP" value="true"/>
        <Descriptor name="ZagrebIndex" value="true"/>
    </Group>
    <Group name="3D">
        <Descriptor name="Autocorrelation3D" value="true"/>
        <Descriptor name="CPSA" value="true"/>
        <Descriptor name="GravitationalIndex" value="true"/>
        <Descriptor name="LengthOverBreadth" value="true"/>
        <Descriptor name="MomentOfInertia" value="true"/>
        <Descriptor name="PetitjeanShapeIndex" value="true"/>
        <Descriptor name="RDF" value="true"/>
        <Descriptor name="WHIM" value="true"/>
    </Group>
    <Group name="Fingerprint">
        <Descriptor name="Fingerprinter" value="false"/>
        <Descriptor name="ExtendedFingerprinter" value="false"/>
        <Descriptor name="EStateFingerprinter" value="false"/>
        <Descriptor name="GraphOnlyFingerprinter" value="false"/>
        <Descriptor name="MACCSFingerprinter" value="false"/>
        <Descriptor name="PubchemFingerprinter" value="false"/>
        <Descriptor name="SubstructureFingerprinter" value="false"/>
        <Descriptor name="SubstructureFingerprintCount" value="false"/>
        <Descriptor name="KlekotaRothFingerprinter" value="true"/>
        <Descriptor name="KlekotaRothFingerprintCount" value="false"/>
        <Descriptor name="AtomPairs2DFingerprinter" value="false"/>
        <Descriptor name="AtomPairs2DFingerprintCount" value="false"/>
    </Group>
</Root>
"""
with open("/tmp/padeldescriptor.xml") as f:
    f.write(config)

SMILES = sys.argv[1]
KR = SMILES.replace("SMILES", "DS").replace("smi", "kr.txt")

print(f"Calculating PaDEL-Descriptor's KR fingerprint for {SMILES}")

# Calculate and clean-up KR descriptors
padeldescriptor(
    mol_dir=SMILES,
    d_file=f"/tmp/{os.path.basename(KR)}",
    descriptortypes="./descriptors.xml",
    d_2d=False,
    d_3d=False,
    fingerprints=True,
    headless=True,
    retainorder=True,
    threads=-1,
)

bits = pd.read_csv(f"/tmp/{os.path.basename(KR)}")
del bits["Name"]
bits.to_csv(KR, header=None, index=None, sep=' ')
