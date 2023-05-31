#General descriptors of molecules
#Read : test.xlsx
#Author: ylwu
#Date: 2023/5/30

import os
import xlrd
import numpy as np
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem import AllChem
from rdkit.Chem import MACCSkeys
from rdkit.Chem import DataStructs


def get_elements(input_file):
    data = xlrd.open_workbook(input_file)  # open excel
    table = data.sheet_by_name("Sheet1")  # read sheet
    nrows = table.nrows  # get the number of lines
    result = []

    for i in range(1, nrows):
        rows = table.row_values(i)  # put the data of line in array
        result.append(rows)

    line_number = np.shape(result)[0]
    out_file = open('RDKit_feature.txt', 'a')

    for line in range(line_number):
        SMILES = result[line][1]
        mol = Chem.MolFromSmiles(SMILES)
        MolWt = AllChem.CalcExactMolWt(mol)
        HeavyAtomMolWt = AllChem.CalcExactMolWt(mol, onlyHeavy = -1)
        HeavyAtomCount = mol.GetNumAtoms()
        NumHeteroatoms = AllChem.CalcNumHeteroatoms(mol)
        NumRotatableBonds = AllChem.CalcNumRotatableBonds(mol)
        FractionCSP3 = AllChem.CalcFractionCSP3(mol)
        Kappa1 = AllChem.CalcKappa1(mol)
        Kappa2 = AllChem.CalcKappa2(mol)
        Kappa3 = AllChem.CalcKappa3(mol)
        NumAromaticCarbocycles = AllChem.CalcNumAromaticCarbocycles(mol)
        NumAromaticRings = AllChem.CalcNumAromaticRings(mol)
        NumAmideBonds = AllChem.CalcNumAmideBonds(mol)
        NumAtomStereoCenters = AllChem.CalcNumAtomStereoCenters(mol)
        NumBridgeheadAtoms = AllChem.CalcNumBridgeheadAtoms(mol)
        NumSaturatedCarbocycles = AllChem.CalcNumSaturatedCarbocycles(mol)
        NumAliphaticCarbocycles = AllChem.CalcNumAliphaticCarbocycles(mol)
        NumAromaticHeterocycles = AllChem.CalcNumAromaticHeterocycles(mol)
        NumAliphaticRings = AllChem.CalcNumAliphaticRings(mol)
        NumLipinskiHBA = AllChem.CalcNumLipinskiHBA(mol)
        NumLipinskiHBD = AllChem.CalcNumLipinskiHBD(mol)
        NumRings = AllChem.CalcNumRings(mol)
        NumSaturatedHeterocycles = AllChem.CalcNumSaturatedHeterocycles(mol)
        NumSaturatedRings = AllChem.CalcNumSaturatedRings(mol)
        NumSpiroAtoms = AllChem.CalcNumSpiroAtoms(mol)
        NumUnspecifiedAtomStereoCenters = AllChem.CalcNumUnspecifiedAtomStereoCenters(mol)
        TPSA = AllChem.CalcTPSA(mol)
        NumHeterocycles = AllChem.CalcNumHeterocycles(mol)
        LabuteASA = AllChem.CalcLabuteASA(mol)



        print(MolWt,
              HeavyAtomMolWt,
              HeavyAtomCount,
              NumHeteroatoms,
              NumRotatableBonds,
              FractionCSP3,
              Kappa1,
              Kappa2,
              Kappa3,
              NumAromaticCarbocycles,
              NumAromaticRings,
              NumAmideBonds,
              NumAtomStereoCenters,
              NumBridgeheadAtoms,
              NumSaturatedCarbocycles,
              NumAliphaticCarbocycles,
              NumAromaticHeterocycles,
              NumAliphaticRings,
              NumLipinskiHBA,
              NumLipinskiHBD,
              NumRings,
              NumSaturatedHeterocycles,
              NumSaturatedRings,
              NumSpiroAtoms,
              NumUnspecifiedAtomStereoCenters,
              TPSA,
              NumHeterocycles,
              LabuteASA,
              file = out_file)

input_file = 'test.xlsx'
get_elements(input_file)



