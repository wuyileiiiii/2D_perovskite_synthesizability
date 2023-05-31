import os
import re
import numpy as np
import warnings
import pandas as pd
import math

CID_file = open("CID.txt",'r',encoding='utf-8')
cid_line = CID_file.readlines()


for i in range(len(cid_line)):
    sdf1 = "Conformer3D_CID_" + str(cid_line[i]) + ".sdf"
    sdf = sdf1.replace('\n','')
    out1 = str(cid_line[i])+ ".gjf"
    out = out1.replace('\n','')
    print(out)
    print(sdf)

    try:
        file_open = open(sdf, "r", encoding="utf-8")
        line_count = len(open(sdf).readlines())
        print(line_count)
        X_all = []
        Y_all = []
        Z_all = []
        E_all = []
        dataArr = []
        for line in file_open.readlines():
            lineArr = []
            curLine = line.strip().split()
            for i in range(len(curLine)):
                lineArr.append(curLine[i])
            dataArr.append(lineArr)

        num_atom = int(dataArr[3][0])
        for i in range(4, 4 + num_atom):
            X_all.append(dataArr[i][0])
            Y_all.append(dataArr[i][1])
            Z_all.append(dataArr[i][2])
            E_all.append(dataArr[i][3])

        tmp = open(out, 'w', encoding="utf-8")
        print('%chk=S0_opt.chk', file = tmp)
        print('%nprocshared=48',file = tmp)
        print('%mem=50GB',file = tmp)
        print('# opt b3lyp/6-31g(d)',file = tmp)
        print('    ',file = tmp)
        print('Title Card Required',file = tmp)
        print('    ',file = tmp)
        print('0 1',file = tmp)
        for i in range(num_atom):
            print(E_all[i], "   ", X_all[i], "   ", Y_all[i], "   ", Z_all[i],  file = tmp)
        print('    ',file = tmp)
        print('    ',file = tmp)
        print('    ',file = tmp)

        tmp.close()

    except:
        out2 = "error_" + str(cid_line[i]) + ".gjf"
        out = out2.replace('\n', '')
        tmp = open(out, 'w', encoding="utf-8")
        print('error', file=tmp)
        tmp.close()













