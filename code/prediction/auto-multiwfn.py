import os
import re
import numpy as np
import warnings
import pandas as pd
import math

CID_file = open("CID.txt",'r',encoding='utf-8')
cid_line = CID_file.readlines()


for i in range(len(cid_line)):
    gjf1 = str(cid_line[i])+ ".gjf"
    gjf = gjf1.replace('\n','')
    out1 = "../result/" + str(cid_line[i])+ ".txt"
    out = out1.replace('\n','')

    print(gjf)
    print(out)

    try:
        path = "../" + gjf
        os.system("./Multiwfn " + path + "<script.txt >" + out)
    except:
        out3 = "../result/" + "error_" + str(cid_line[i]) + ".txt"
        out2 = out3.replace('\n', '')
        tmp = open(out2, 'w', encoding="utf-8")
        print('error', file=tmp)
        tmp.close()













