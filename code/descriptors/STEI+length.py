import os
import re
import numpy as np
import warnings
import pandas as pd
import math
import rdkit
from rdkit import Chem

#STEI and length of molecules
#STEI：the steric effect index of nitrogen. When there are multiple nitrogen atoms, take the maximum value.
#Eccentricity：the maximum value in the row of nitrogen in the distance matrix. When there are multiple nitrogen atoms, take the maximum value.
#
#python：3.7
#author：ylwu
#Date：2022/6/7


def get_STEI(smile):
    mol = Chem.MolFromSmiles(smile)

    #将整个分子中每个原子对应的编号和原子类型进行记录
    Symbol = []
    ID = []
    for atom in mol.GetAtoms():
        id = atom.GetIdx()
        symbol = atom.GetSymbol()
        Symbol.append(symbol)
        ID.append(id)

    STEI = []
    N_id = []
    length = []
    for atom in mol.GetAtoms():
        id = atom.GetIdx()
        if atom.GetAtomicNum() == 7:
            N_id.append(id)

    for i in N_id:
        print(i)
        ID_it = ID.copy()
        ID_it.remove(i)
        print(ID, i, ID_it)

        #find the ID of nitrogen
        distance = []
        Neighbor = []
        N = mol.GetAtomWithIdx(i)

        #Find the first neighboring atoms of nitrogen
        atom1 = []
        x1 = [i.GetIdx() for i in N.GetNeighbors()]
        for i in x1:
            d = 1
            if i in ID_it:
                atom1.append(i)
                Neighbor.append(i)
                distance.append(d)
                ID_it.remove(i)
            else:
                continue
        print("distance:", distance)

        #Find the second neighboring atoms of nitrogen
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
                else:
                    continue
            print("x2:", x2)
        print("distance:", distance)

        #Find the third neighboring atoms of nitrogen
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
                else:
                    continue
            print("x3:", x3)
        print("distance:", distance)

        #Find the 4th neighboring atoms of nitrogen
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
                else:
                    continue
            print("x4:", x4)
        print("distance:", distance)

        #Find the 5th neighboring atoms of nitrogen
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
                else:
                    continue
            print("x5:", x5)
        print("distance:", distance)

        #Find the 6th neighboring atoms of nitrogen
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
                else:
                    continue
            print("x6:", x6)
        print("distance:", distance)

        #Find the 7th neighboring atoms of nitrogen
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
                else:
                    continue
            print("x7:", x7)
        print("distance:", distance)

        #Find the 8th neighboring atoms of nitrogen
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
                else:
                    continue
            print("x8:", x8)
        print("distance:", distance)

        #Find the 9th neighboring atoms of nitrogen
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
                else:
                    continue
            print("x9:", x9)
        print("distance:", distance)

        #Find the 10th neighboring atoms of nitrogen
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
                else:
                    continue
            print("x10:", x10)
        print("distance:", distance)

        #Find the 11th neighboring atoms of nitrogen
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
                else:
                    continue
            print("x11:", x11)
        print("distance:", distance)

        #Find the 12th neighboring atoms of nitrogen
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
                else:
                    continue
            print("x12:", x12)
        print("distance:", distance)

        #Find the 13th neighboring atoms of nitrogen
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
                else:
                    continue
            print("x13:", x13)
        print("distance:", distance)

        #Find the 14th neighboring atoms of nitrogen
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
                else:
                    continue
            print("x14:", x14)
        print("distance:", distance)

        #Find the 15th neighboring atoms of nitrogen
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
                else:
                    continue
            print("x15:", x15)
        print("distance:", distance)

        #Find the 16th neighboring atoms of nitrogen
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
                else:
                    continue
            print("x16:", x16)
        print("distance:", distance)

        #Find the 17th neighboring atoms of nitrogen
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
                else:
                    continue
            print("x17:", x17)
        print("distance:", distance)

        #Find the 18th neighboring atoms of nitrogen
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
                else:
                    continue
            print("x18:", x18)
        print("distance:", distance)

        #Find the 19th neighboring atoms of nitrogen
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
                else:
                    continue
            print("x19:", x19)
        print("distance:", distance)

        #Find the 20th neighboring atoms of nitrogen
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
                else:
                    continue
            print("x20:", x20)
        print("distance:", distance)

        stei = 0
        for i in distance:
            stei += 1/(i*i*i)
        STEI.append(stei)


        len = max(distance)
        length.append(len)
    STEI_max = max(STEI)
    length_max = max(length)


    return STEI_max, length_max


smiles = open("smiles.txt",'r', encoding='utf-8')
STEI_file = open("STEI+length.txt", 'w', encoding='utf-8')
for smile in smiles.readlines():
    STEI, length = get_STEI(smile)
    print(STEI, length, file=STEI_file)
