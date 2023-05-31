import os
import re
import numpy as np
import warnings
import pandas as pd
import math

CID_file = open("CID.txt",'r',encoding='utf-8')
cid_line = CID_file.readlines()
file_out = open('change.txt', 'w')

for i in range(len(cid_line)):
    sdf1 = "Conformer3D_CID_" + str(cid_line[i]) + ".sdf"
    sdf = sdf1.replace('\n','')
    out1 = str(cid_line[i])+ ".gjf"
    out = out1.replace('\n','')


    try:
        file_open = open(sdf, "r", encoding="utf-8")
        line_count = len(open(sdf).readlines())
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
            X_all.append(float(dataArr[i][0]))
            Y_all.append(float(dataArr[i][1]))
            Z_all.append(float(dataArr[i][2]))
            E_all.append(dataArr[i][3])

#X_max
        a = str(np.where(X_all==np.max(X_all)))
        n1 = a.replace('(array([','')
        n2 = n1.replace('], dtype=int64),)','')
        X_max = dataArr[int(n2)+4][0]
        E_X_max = dataArr[int(n2)+4][3]
        if E_X_max == 'H':
            C_X_max = -0.78
        elif E_X_max == 'C':
            C_X_max = -0.95
        elif E_X_max == 'N':
            C_X_max = -0.84
        elif E_X_max == 'O':
            C_X_max = -0.88
        elif E_X_max == 'F':
            C_X_max = -0.87
        else:
            C_X_max = 'NaN'
#       print("X_max: ", X_max,"E_X_max: ", E_X_max, "C_X_max: ", C_X_max)

#X_min
        a = str(np.where(X_all==np.min(X_all)))
        n1 = a.replace('(array([','')
        n2 = n1.replace('], dtype=int64),)','')
        X_min = dataArr[int(n2)+4][0]
        E_X_min = dataArr[int(n2)+4][3]
        if E_X_min == 'H':
            C_X_min = -0.78
        elif E_X_min == 'C':
            C_X_min = -0.95
        elif E_X_min == 'N':
            C_X_min = -0.84
        elif E_X_min == 'O':
            C_X_min = -0.88
        elif E_X_min == 'F':
            C_X_min = -0.87
        else:
            C_X_min = 'NaN'
#       print("X_min: ", X_min,"E_X_min: ", E_X_min, "C_X_min: ", C_X_min)

#Y_max
        a = str(np.where(Y_all==np.max(Y_all)))
        n1 = a.replace('(array([','')
        n2 = n1.replace('], dtype=int64),)','')
        Y_max = dataArr[int(n2)+4][0]
        E_Y_max = dataArr[int(n2)+4][3]
        if E_Y_max == 'H':
            C_Y_max = -0.78
        elif E_Y_max == 'C':
            C_Y_max = -0.95
        elif E_Y_max == 'N':
            C_Y_max = -0.84
        elif E_Y_max == 'O':
            C_Y_max = -0.88
        elif E_Y_max == 'F':
            C_Y_max = -0.87
        else:
            C_Y_max = 'NaN'
#       print("Y_max: ", Y_max,"E_Y_max: ", E_Y_max, "C_Y_max: ", C_Y_max)

#Y_min
        a = str(np.where(Y_all==np.min(Y_all)))
        n1 = a.replace('(array([','')
        n2 = n1.replace('], dtype=int64),)','')
        Y_min = dataArr[int(n2)+4][0]
        E_Y_min = dataArr[int(n2)+4][3]
        if E_Y_min == 'H':
            C_Y_min = -0.78
        elif E_Y_min == 'C':
            C_Y_min = -0.95
        elif E_Y_min == 'N':
            C_Y_min = -0.84
        elif E_Y_min == 'O':
            C_Y_min = -0.88
        elif E_Y_min == 'F':
            C_Y_min = -0.87
        else:
            C_Y_min = 'NaN'
#       print("Y_min: ", Y_min,"E_Y_min: ", E_Y_min, "C_Y_min: ", C_Y_min)

#Z_max
        a = str(np.where(Z_all==np.max(Z_all)))
        n1 = a.replace('(array([','')
        n2 = n1.replace('], dtype=int64),)','')
        Z_max = dataArr[int(n2)+4][0]
        E_Z_max = dataArr[int(n2)+4][3]
        if E_Z_max == 'H':
            C_Z_max = -0.78
        elif E_Z_max == 'C':
            C_Z_max = -0.95
        elif E_Z_max == 'N':
            C_Z_max = -0.84
        elif E_Z_max == 'O':
            C_Z_max = -0.88
        elif E_Z_max == 'F':
            C_Z_max = -0.87
        else:
            C_Z_max = 'NaN'
#       print("Z_max: ", Z_max,"E_Z_max: ", E_Z_max, "C_Z_max: ", C_Z_max)

#Z_min
        a = str(np.where(Z_all==np.min(Z_all)))
        n1 = a.replace('(array([','')
        n2 = n1.replace('], dtype=int64),)','')
        Z_min = dataArr[int(n2)+4][0]
        E_Z_min = dataArr[int(n2)+4][3]
        if E_Z_min == 'H':
            C_Z_min = -0.78
        elif E_Z_min == 'C':
            C_Z_min = -0.95
        elif E_Z_min == 'N':
            C_Z_min = -0.84
        elif E_Z_min == 'O':
            C_Z_min = -0.88
        elif E_Z_min == 'F':
            C_Z_min = -0.87
        else:
            C_Z_min = 'NaN'
#       print("Z_min: ", Z_min,"E_Z_min: ", E_Z_min, "C_Z_min: ", C_Z_min)

        print("X_max: ", X_max, "E_X_max: ", E_X_max, "C_X_max: ", C_X_max,
              "X_min: ", X_min,"E_X_min: ", E_X_min, "C_X_min: ", C_X_min,
              "Y_max: ", Y_max, "E_Y_max: ", E_Y_max, "C_Y_max: ", C_Y_max,
              "Y_min: ", Y_min, "E_Y_min: ", E_Y_min, "C_Y_min: ", C_Y_min,
              "Z_max: ", Z_max, "E_Z_max: ", E_Z_max, "C_Z_max: ", C_Z_max,
              "Z_min: ", Z_min, "E_Z_min: ", E_Z_min, "C_Z_min: ", C_Z_min,
               file = file_out)

    except:
            print("error", file = file_out)

file_out.close()











