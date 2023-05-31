import os
import re
import numpy as np
import warnings
import pandas as pd
import math
import rdkit
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw import rdMolDraw2D
from IPython.display import SVG

#Calculate the smallest distance between nitrogens
#When molecules do not contain nitrogens or only contain one nitrogen, DisNN=0
#When molecules contain more than three nitrogens, DisNN refers to the smallest value of distance
#Only suitable for molecules with up to 20 nerghboring non-hydrogen atoms
#python：3.7
#author：ylwu
#Date：2022/6/1



def get_disNN(smile):
    mol = Chem.MolFromSmiles(smile)
    #Record the ID and type of each atom in the molecule
    Symbol = []
    ID = []
    for atom in mol.GetAtoms():
        id = atom.GetIdx()
        symbol = atom.GetSymbol()
        Symbol.append(symbol)
        ID.append(id)

    N_num = 0
    for atom in mol.GetAtoms():
        if atom.GetAtomicNum() == 7:
            N_num += 1

    if N_num == 0:
        NN_is = 0

    elif N_num == 1:
        NN_is = 0

    elif N_num == 2:
        N_id = []
        iis_array = []
        #Find the ID of nitrogen

        for atom in mol.GetAtoms():
            id = atom.GetIdx()
            if atom.GetAtomicNum() == 7:
                N_id.append(id)

        for i in N_id:
            print(i)
            ID_it = ID.copy()
            ID_it.remove(i)
            print(ID, i, ID_it)
            #Find the first neighboring atoms of nitrogen
            distance = []
            Neighbor = []
            N = mol.GetAtomWithIdx(i)

            atom1 = []
            x1 = [i.GetIdx() for i in N.GetNeighbors()]
            print("i:", i)
            print("x1:", x1)
            for i in x1:
                d = 1
                if i in ID_it:
                    atom1.append(i)
                    Neighbor.append(i)
                    distance.append(d)
                    ID_it.remove(i)
                    if mol.GetAtomWithIdx(i).GetAtomicNum() == 7:
                        iis = d
                        iis_array.append(iis)
                    else:
                        continue
                else:
                    continue

            #Find the second neighboring atoms of N
            atom2 = []
            for i in atom1:
                d = 2
                a = mol.GetAtomWithIdx(i)
                x2 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x2:
                    if j in ID_it:
                        atom2.append(j)
                        Neighbor.append(j)
                        distance.append(d)
                        ID_it.remove(j)
                        if mol.GetAtomWithIdx(j).GetAtomicNum() == 7:
                            iis = d
                            iis_array.append(iis)
                        else:
                            continue
                    else:
                        continue
                else:
                    continue

            #Find the third neighboring atoms of N
            atom3 = []
            for i in atom2:
                d = 3
                a = mol.GetAtomWithIdx(i)
                x3 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x3:
                    if j in ID_it:
                        atom3.append(j)
                        Neighbor.append(j)
                        distance.append(d)
                        ID_it.remove(j)
                        if mol.GetAtomWithIdx(j).GetAtomicNum() == 7:
                            iis = d
                            iis_array.append(iis)
                        else:
                            continue
                    else:
                        continue
                else:
                    continue

            #Find the fourth neighboring atoms of N
            atom4 = []
            for i in atom3:
                d = 4
                a = mol.GetAtomWithIdx(i)
                x4 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x4:
                    if j in ID_it:
                        atom4.append(j)
                        Neighbor.append(j)
                        distance.append(d)
                        ID_it.remove(j)
                        if mol.GetAtomWithIdx(j).GetAtomicNum() == 7:
                            iis = d
                            iis_array.append(iis)
                        else:
                            continue
                    else:
                        continue
                else:
                    continue

            #Find the 5th neighboring atoms of N
            atom5 = []
            for i in atom4:
                d = 5
                a = mol.GetAtomWithIdx(i)
                x5 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x5:
                    if j in ID_it:
                        atom5.append(j)
                        Neighbor.append(j)
                        distance.append(d)
                        ID_it.remove(j)
                        if mol.GetAtomWithIdx(j).GetAtomicNum() == 7:
                            iis = d
                            iis_array.append(iis)
                        else:
                            continue
                    else:
                        continue
                else:
                    continue

            #Find the 6th neighboring atoms of N
            atom6 = []
            for i in atom5:
                d = 6
                a = mol.GetAtomWithIdx(i)
                x6 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x6:
                    if j in ID_it:
                        atom6.append(j)
                        Neighbor.append(j)
                        distance.append(d)
                        ID_it.remove(j)
                        if mol.GetAtomWithIdx(j).GetAtomicNum() == 7:
                            iis = d
                            iis_array.append(iis)
                        else:
                            continue
                    else:
                        continue
                else:
                    continue

            #Find the 7th neighboring atoms of N
            atom7 = []
            for i in atom6:
                d = 7
                a = mol.GetAtomWithIdx(i)
                x7 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x7:
                    if j in ID_it:
                        atom7.append(j)
                        Neighbor.append(j)
                        distance.append(d)
                        ID_it.remove(j)
                        if mol.GetAtomWithIdx(j).GetAtomicNum() == 7:
                            iis = d
                            iis_array.append(iis)
                        else:
                            continue
                    else:
                        continue
                else:
                    continue

            #Find the 8th neighboring atoms of N
            atom8 = []
            for i in atom7:
                d = 8
                a = mol.GetAtomWithIdx(i)
                x8 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x8:
                    if j in ID_it:
                        atom8.append(j)
                        Neighbor.append(j)
                        distance.append(d)
                        ID_it.remove(j)
                        if mol.GetAtomWithIdx(j).GetAtomicNum() == 7:
                            iis = d
                            iis_array.append(iis)
                        else:
                            continue
                    else:
                        continue
                else:
                    continue

            #Find the 9th neighboring atoms of N
            atom9 = []
            for i in atom8:
                d = 9
                a = mol.GetAtomWithIdx(i)
                x9 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x9:
                    if j in ID_it:
                        atom9.append(j)
                        Neighbor.append(j)
                        distance.append(d)
                        ID_it.remove(j)
                        if mol.GetAtomWithIdx(j).GetAtomicNum() == 7:
                            iis = d
                            iis_array.append(iis)
                        else:
                            continue
                    else:
                        continue
                else:
                    continue

            #Find the 10th neighboring atoms of N
            atom10 = []
            for i in atom9:
                d = 10
                a = mol.GetAtomWithIdx(i)
                x10 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x10:
                    if j in ID_it:
                        atom10.append(j)
                        Neighbor.append(j)
                        distance.append(d)
                        ID_it.remove(j)
                        if mol.GetAtomWithIdx(j).GetAtomicNum() == 7:
                            iis = d
                            iis_array.append(iis)
                        else:
                            continue
                    else:
                        continue
                else:
                    continue

            #Find the 11th neighboring atoms of N
            atom11 = []
            for i in atom10:
                d = 11
                a = mol.GetAtomWithIdx(i)
                x11 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x11:
                    if j in ID_it:
                        atom11.append(j)
                        Neighbor.append(j)
                        distance.append(d)
                        ID_it.remove(j)
                        if mol.GetAtomWithIdx(j).GetAtomicNum() == 7:
                            iis = d
                            iis_array.append(iis)
                        else:
                            continue
                    else:
                        continue
                else:
                    continue

            #Find the 12th neighboring atoms of N
            atom12 = []
            for i in atom11:
                d = 12
                a = mol.GetAtomWithIdx(i)
                x12 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x12:
                    if j in ID_it:
                        atom12.append(j)
                        Neighbor.append(j)
                        distance.append(d)
                        ID_it.remove(j)
                        if mol.GetAtomWithIdx(j).GetAtomicNum() == 7:
                            iis = d
                            iis_array.append(iis)
                        else:
                            continue
                    else:
                        continue
                else:
                    continue

            #Find the 13th neighboring atoms of N
            atom13 = []
            for i in atom12:
                d = 13
                a = mol.GetAtomWithIdx(i)
                x13 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x13:
                    if j in ID_it:
                        atom13.append(j)
                        Neighbor.append(j)
                        distance.append(d)
                        ID_it.remove(j)
                        if mol.GetAtomWithIdx(j).GetAtomicNum() == 7:
                            iis = d
                            iis_array.append(iis)
                        else:
                            continue
                    else:
                        continue
                else:
                    continue

            #寻Find the 14th neighboring atoms of N
            atom14 = []
            for i in atom13:
                d = 14
                a = mol.GetAtomWithIdx(i)
                x14 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x14:
                    if j in ID_it:
                        atom14.append(j)
                        Neighbor.append(j)
                        distance.append(d)
                        ID_it.remove(j)
                        if mol.GetAtomWithIdx(j).GetAtomicNum() == 7:
                            iis = d
                            iis_array.append(iis)
                        else:
                            continue
                    else:
                        continue
                else:
                    continue

            #Find the 15th neighboring atoms of N
            atom15 = []
            for i in atom14:
                d = 15
                a = mol.GetAtomWithIdx(i)
                x15 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x15:
                    if j in ID_it:
                        atom15.append(j)
                        Neighbor.append(j)
                        distance.append(d)
                        ID_it.remove(j)
                        if mol.GetAtomWithIdx(j).GetAtomicNum() == 7:
                            iis = d
                            iis_array.append(iis)
                        else:
                            continue
                    else:
                        continue
                else:
                    continue

            #Find the 16th neighboring atoms of N
            atom16 = []
            for i in atom15:
                d = 16
                a = mol.GetAtomWithIdx(i)
                x16 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x16:
                    if j in ID_it:
                        atom16.append(j)
                        Neighbor.append(j)
                        distance.append(d)
                        ID_it.remove(j)
                        if mol.GetAtomWithIdx(j).GetAtomicNum() == 7:
                            iis = d
                            iis_array.append(iis)
                        else:
                            continue
                    else:
                        continue
                else:
                    continue

            #Find the 17th neighboring atoms of N
            atom17 = []
            for i in atom16:
                d = 17
                a = mol.GetAtomWithIdx(i)
                x17 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x17:
                    if j in ID_it:
                        atom17.append(j)
                        Neighbor.append(j)
                        distance.append(d)
                        ID_it.remove(j)
                        if mol.GetAtomWithIdx(j).GetAtomicNum() == 7:
                            iis = d
                            iis_array.append(iis)
                        else:
                            continue
                    else:
                        continue
                else:
                    continue

            #Find the 18th neighboring atoms of N
            atom18 = []
            for i in atom17:
                d = 18
                a = mol.GetAtomWithIdx(i)
                x18 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x18:
                    if j in ID_it:
                        atom18.append(j)
                        Neighbor.append(j)
                        distance.append(d)
                        ID_it.remove(j)
                        if mol.GetAtomWithIdx(j).GetAtomicNum() == 7:
                            iis = d
                            iis_array.append(iis)
                        else:
                            continue
                    else:
                        continue
                else:
                    continue

            #Find the 19th neighboring atoms of N
            atom19 = []
            for i in atom18:
                d = 19
                a = mol.GetAtomWithIdx(i)
                x19 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x19:
                    if j in ID_it:
                        atom19.append(j)
                        Neighbor.append(j)
                        distance.append(d)
                        ID_it.remove(j)
                        if mol.GetAtomWithIdx(j).GetAtomicNum() == 7:
                            iis = d
                            iis_array.append(iis)
                        else:
                            continue
                    else:
                        continue
                else:
                    continue

            #Find the 20th neighboring atoms of N
            atom20 = []
            for i in atom19:
                d = 20
                a = mol.GetAtomWithIdx(i)
                x20 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x20:
                    if j in ID_it:
                        atom20.append(j)
                        Neighbor.append(j)
                        distance.append(d)
                        ID_it.remove(j)
                        if mol.GetAtomWithIdx(j).GetAtomicNum() == 7:
                            iis = d
                            iis_array.append(iis)
                        else:
                            continue
                    else:
                        continue
                else:
                    continue
        print("iis_array:", iis_array)
        NN_is = 1/(max(iis_array)*max(iis_array))

    else:
        NN_is = "NaN"

    return NN_is

smiles = open("smiles.txt",'r', encoding='utf-8')
disNN_file = open("disNN.txt", 'w', encoding='utf-8')
for smile in smiles.readlines():
    disNN = get_disNN(smile)
    print(disNN, file=disNN_file)
