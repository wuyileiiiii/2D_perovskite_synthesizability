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

#Calculate the number of rotational bonds in the alkyl tail
#When molecules contain more than one nitrogen atoms, take the maximum value
#Only suitable for molecules with up to 20 nerghboring non-hydrogen atoms
#python：3.7
#author：ylwu
#Date：2022/6/10



def get_Rotation(smile):
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
    N_id = []
    trash = open("trash.txt",'w')
    for atom in mol.GetAtoms():
        id = atom.GetIdx()
        if atom.GetAtomicNum() == 7:
            N_num += 1
            N_id.append(id)

    Rotation = []
    for i in N_id:
        ID_it = ID.copy()
        ID_it.remove(i)
        #Find the ID of nitrogen
        distance = []
        Neighbor = []
        N = mol.GetAtomWithIdx(i)

        #Determine whether nitrogen atoms are on the ring
        if mol.GetAtomWithIdx(i).IsInRing():
            Num_rotation = 0
        else:
            Num_rotation = 1

        #Find the first neighboring atoms of N
        atom1 = []
        x1 = [i.GetIdx() for i in N.GetNeighbors()]
        ss = 0
        for i in x1:
            d = 1
            if i in ID_it:
                if mol.GetAtomWithIdx(i).IsInRing():
                    ss += 1
                else:
                    atom1.append(i)
                    Neighbor.append(i)
                    distance.append(d)
                    ID_it.remove(i)
            else:
                print("found NN", file=trash)

        if ss == 0:
            judge = 0
            Num_rotation = Num_rotation + 1
        else:
            judge = 1

        #Find the second neighboring atoms of N
        atom2 = []
        if judge == 1:
            print("found", file=trash)
        else:
            ss = 0
            for i in atom1:
                d = 2
                a = mol.GetAtomWithIdx(i)
                x2 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x2:
                    if j in ID_it:
                        if mol.GetAtomWithIdx(j).IsInRing():
                            ss += 1
                        else:
                            atom2.append(j)
                            Neighbor.append(j)
                            distance.append(d)
                            ID_it.remove(j)
                    else:
                        print("found NN", file=trash)
            if ss == 0:
                judge = 0
                Num_rotation = Num_rotation + 1
            else:
                judge = 1

        #Find the third neighboring atoms of N
        atom3 = []
        if judge == 1:
            print("found", file=trash)
        else:
            ss = 0
            for i in atom2:
                d = 3
                a = mol.GetAtomWithIdx(i)
                x3 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x3:
                    if j in ID_it:
                        if mol.GetAtomWithIdx(j).IsInRing():
                            ss += 1
                        else:
                            atom3.append(j)
                            Neighbor.append(j)
                            distance.append(d)
                            ID_it.remove(j)
                    else:
                        print("found NN", file=trash)
            if ss == 0:
                judge = 0
                Num_rotation = Num_rotation + 1
            else:
                judge = 1

        #Find the fourth neighboring atoms of N
        atom4 = []
        if judge == 1:
            print("found", file=trash)
        else:
            ss = 0
            for i in atom3:
                d = 4
                a = mol.GetAtomWithIdx(i)
                x4 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x4:
                    if j in ID_it:
                        if mol.GetAtomWithIdx(j).IsInRing():
                            ss += 1
                        else:
                            atom4.append(j)
                            Neighbor.append(j)
                            distance.append(d)
                            ID_it.remove(j)
                    else:
                        print("found NN", file=trash)
            if ss == 0:
                judge = 0
                Num_rotation = Num_rotation + 1
            else:
                judge = 1

        #Find the 5th neighboring atoms of N
        atom5 = []
        if judge == 1:
            print("found", file=trash)
        else:
            ss = 0
            for i in atom4:
                d = 5
                a = mol.GetAtomWithIdx(i)
                x5 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x5:
                    if j in ID_it:
                        if mol.GetAtomWithIdx(j).IsInRing():
                            ss += 1
                        else:
                            atom5.append(j)
                            Neighbor.append(j)
                            distance.append(d)
                            ID_it.remove(j)
                    else:
                        print("found NN", file=trash)
            if ss == 0:
                judge = 0
                Num_rotation = Num_rotation + 1
            else:
                judge = 1

        #Find the 6th neighboring atoms of N
        atom6 = []
        if judge == 1:
            print("found", file=trash)
        else:
            ss = 0
            for i in atom5:
                d = 6
                a = mol.GetAtomWithIdx(i)
                x6 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x6:
                    if j in ID_it:
                        if mol.GetAtomWithIdx(j).IsInRing():
                            ss += 1
                        else:
                            atom6.append(j)
                            Neighbor.append(j)
                            distance.append(d)
                            ID_it.remove(j)
                    else:
                        print("found NN", file=trash)
            if ss == 0:
                judge = 0
                Num_rotation = Num_rotation + 1
            else:
                judge = 1

        #Find the 7th neighboring atoms of N
        atom7 = []
        if judge == 1:
            print("found", file=trash)
        else:
            ss = 0
            for i in atom6:
                d = 7
                a = mol.GetAtomWithIdx(i)
                x7 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x7:
                    if j in ID_it:
                        if mol.GetAtomWithIdx(j).IsInRing():
                            ss += 1
                        else:
                            atom7.append(j)
                            Neighbor.append(j)
                            distance.append(d)
                            ID_it.remove(j)
                    else:
                        print("found NN", file=trash)
            if ss == 0:
                judge = 0
                Num_rotation = Num_rotation + 1
            else:
                judge = 1

        #Find the 8th neighboring atoms of N
        atom8 = []
        if judge == 1:
            print("found", file=trash)
        else:
            ss = 0
            for i in atom7:
                d = 8
                a = mol.GetAtomWithIdx(i)
                x8 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x8:
                    if j in ID_it:
                        if mol.GetAtomWithIdx(j).IsInRing():
                            ss += 1
                        else:
                            atom8.append(j)
                            Neighbor.append(j)
                            distance.append(d)
                            ID_it.remove(j)
                    else:
                        print("found NN", file=trash)
            if ss == 0:
                judge = 0
                Num_rotation = Num_rotation + 1
            else:
                judge = 1

        #Find the 9th neighboring atoms of N
        atom9 = []
        if judge == 1:
            print("found", file=trash)
        else:
            ss = 0
            for i in atom8:
                d = 9
                a = mol.GetAtomWithIdx(i)
                x9 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x9:
                    if j in ID_it:
                        if mol.GetAtomWithIdx(j).IsInRing():
                            ss += 1
                        else:
                            atom9.append(j)
                            Neighbor.append(j)
                            distance.append(d)
                            ID_it.remove(j)
                    else:
                        print("found NN", file=trash)
            if ss == 0:
                judge = 0
                Num_rotation = Num_rotation + 1
            else:
                judge = 1

        #Find the 10th neighboring atoms of N
        atom10 = []
        if judge == 1:
            print("found", file=trash)
        else:
            ss = 0
            for i in atom9:
                d = 10
                a = mol.GetAtomWithIdx(i)
                x10 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x10:
                    if j in ID_it:
                        if mol.GetAtomWithIdx(j).IsInRing():
                            ss += 1
                        else:
                            atom10.append(j)
                            Neighbor.append(j)
                            distance.append(d)
                            ID_it.remove(j)
                    else:
                        print("found NN", file=trash)
            if ss == 0:
                judge = 0
                Num_rotation = Num_rotation + 1
            else:
                judge = 1

        #Find the 11th neighboring atoms of N
        atom11 = []
        if judge == 1:
            print("found", file=trash)
        else:
            ss = 0
            for i in atom10:
                d = 11
                a = mol.GetAtomWithIdx(i)
                x11 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x11:
                    if j in ID_it:
                        if mol.GetAtomWithIdx(j).IsInRing():
                            ss += 1
                        else:
                            atom11.append(j)
                            Neighbor.append(j)
                            distance.append(d)
                            ID_it.remove(j)
                    else:
                        print("found NN", file=trash)
            if ss == 0:
                judge = 0
                Num_rotation = Num_rotation + 1
            else:
                judge = 1

        #Find the 12th neighboring atoms of N
        atom12 = []
        if judge == 1:
            print("found", file=trash)
        else:
            ss = 0
            for i in atom11:
                d = 12
                a = mol.GetAtomWithIdx(i)
                x12 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x12:
                    if j in ID_it:
                        if mol.GetAtomWithIdx(j).IsInRing():
                            ss += 1
                        else:
                            atom12.append(j)
                            Neighbor.append(j)
                            distance.append(d)
                            ID_it.remove(j)
                    else:
                        print("found NN", file=trash)
            if ss == 0:
                judge = 0
                Num_rotation = Num_rotation + 1
            else:
                judge = 1

        #Find the 13th neighboring atoms of N
        atom13 = []
        if judge == 1:
            print("found", file=trash)
        else:
            ss = 0
            for i in atom12:
                d = 13
                a = mol.GetAtomWithIdx(i)
                x13 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x13:
                    if j in ID_it:
                        if mol.GetAtomWithIdx(j).IsInRing():
                            ss += 1
                        else:
                            atom13.append(j)
                            Neighbor.append(j)
                            distance.append(d)
                            ID_it.remove(j)
                    else:
                        print("found NN", file=trash)
            if ss == 0:
                judge = 0
                Num_rotation = Num_rotation + 1
            else:
                judge = 1

        #Find the 14th neighboring atoms of N
        atom14 = []
        if judge == 1:
            print("found", file=trash)
        else:
            ss = 0
            for i in atom13:
                d = 14
                a = mol.GetAtomWithIdx(i)
                x14 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x14:
                    if j in ID_it:
                        if mol.GetAtomWithIdx(j).IsInRing():
                            ss += 1
                        else:
                            atom14.append(j)
                            Neighbor.append(j)
                            distance.append(d)
                            ID_it.remove(j)
                    else:
                        print("found NN", file=trash)
            if ss == 0:
                judge = 0
                Num_rotation = Num_rotation + 1
            else:
                judge = 1

        #Find the 15th neighboring atoms of N
        atom15 = []
        if judge == 1:
            print("found", file=trash)
        else:
            ss = 0
            for i in atom14:
                d = 15
                a = mol.GetAtomWithIdx(i)
                x15 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x15:
                    if j in ID_it:
                        if mol.GetAtomWithIdx(j).IsInRing():
                            ss += 1
                        else:
                            atom15.append(j)
                            Neighbor.append(j)
                            distance.append(d)
                            ID_it.remove(j)
                    else:
                        print("found NN", file=trash)
            if ss == 0:
                judge = 0
                Num_rotation = Num_rotation + 1
            else:
                judge = 1

        #Find the 16th neighboring atoms of N
        atom16 = []
        if judge == 1:
            print("found", file=trash)
        else:
            ss = 0
            for i in atom15:
                d = 16
                a = mol.GetAtomWithIdx(i)
                x16 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x16:
                    if j in ID_it:
                        if mol.GetAtomWithIdx(j).IsInRing():
                            ss += 1
                        else:
                            atom16.append(j)
                            Neighbor.append(j)
                            distance.append(d)
                            ID_it.remove(j)
                    else:
                        print("found NN", file=trash)
            if ss == 0:
                judge = 0
                Num_rotation = Num_rotation + 1
            else:
                judge = 1

        #Find the 17th neighboring atoms of N
        atom17 = []
        if judge == 1:
            print("found", file=trash)
        else:
            ss = 0
            for i in atom16:
                d = 17
                a = mol.GetAtomWithIdx(i)
                x17 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x17:
                    if j in ID_it:
                        if mol.GetAtomWithIdx(j).IsInRing():
                            ss += 1
                        else:
                            atom17.append(j)
                            Neighbor.append(j)
                            distance.append(d)
                            ID_it.remove(j)
                    else:
                        print("found NN", file=trash)
            if ss == 0:
                judge = 0
                Num_rotation = Num_rotation + 1
            else:
                judge = 1

        #Find the 18th neighboring atoms of N
        atom18 = []
        if judge == 1:
            print("found", file=trash)
        else:
            ss = 0
            for i in atom17:
                d = 18
                a = mol.GetAtomWithIdx(i)
                x18 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x18:
                    if j in ID_it:
                        if mol.GetAtomWithIdx(j).IsInRing():
                            ss += 1
                        else:
                            atom18.append(j)
                            Neighbor.append(j)
                            distance.append(d)
                            ID_it.remove(j)
                    else:
                        print("found NN", file=trash)
            if ss == 0:
                judge = 0
                Num_rotation = Num_rotation + 1
            else:
                judge = 1

        #Find the 19th neighboring atoms of N
        atom19 = []
        if judge == 1:
            print("found", file=trash)
        else:
            ss = 0
            for i in atom18:
                d = 19
                a = mol.GetAtomWithIdx(i)
                x19 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x19:
                    if j in ID_it:
                        if mol.GetAtomWithIdx(j).IsInRing():
                            ss += 1
                        else:
                            atom19.append(j)
                            Neighbor.append(j)
                            distance.append(d)
                            ID_it.remove(j)
                    else:
                        print("found NN", file=trash)
            if ss == 0:
                judge = 0
                Num_rotation = Num_rotation + 1
            else:
                judge = 1

        #Find the 20th neighboring atoms of N
        atom20 = []
        if judge == 1:
            print("found", file=trash)
        else:
            ss = 0
            for i in atom19:
                d = 20
                a = mol.GetAtomWithIdx(i)
                x20 = [i.GetIdx() for i in a.GetNeighbors()]
                for j in x20:
                    if j in ID_it:
                        if mol.GetAtomWithIdx(j).IsInRing():
                            ss += 1
                        else:
                            atom20.append(j)
                            Neighbor.append(j)
                            distance.append(d)
                            ID_it.remove(j)
                    else:
                        print("found NN", file=trash)
            if ss == 0:
                judge = 0
                Num_rotation = Num_rotation + 1
            else:
                judge = 1

        Rotation.append(Num_rotation)
    Num_Rotation = max(Rotation)
    trash.close()
    return N_num, Num_Rotation


smiles = open("smiles.txt",'r', encoding='utf-8')
Rotation_file = open("Rotation.txt", 'w', encoding='utf-8')
for smile in smiles.readlines():
    N_num, Num_Rotation = get_Rotation(smile)
    print(N_num, Num_Rotation, file=Rotation_file)

