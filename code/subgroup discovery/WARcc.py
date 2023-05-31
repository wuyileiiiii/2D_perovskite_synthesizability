#Calculate the WARcc of each subgroup
#Author : ylwu
#2023/5/30


import numpy as np
import time

start=time.time()

inp=open('train.txt','r').readlines()
X0=[]
y=np.array([])
for i in range(1,len(inp)):   # from the second line to the last line
	line=str(inp[i]).split()
	feat=[]
	for j in range(2,len(line)):  # line[0]: name; line[1]: property; line[2] first feature
		feat.append(float(line[j]))
	X0.append(feat)
	y=np.append(y,float(line[1]))
X=np.array(X0)

Num_total_positive = 0
Num_total_negative = 0
Num_total = len(y)
print(Num_total)
sum = 0
for i in range(Num_total):
	sum += y[i]
Num_total_positive = sum
Num_total_negative = Num_total - Num_total_positive
print(Num_total_positive, Num_total_negative)

def subgroup(subx0,subx1,subx2,subx3):
	print(subx0,subx1,subx2,subx3)
	sub_x = []
	sub_y = []
	print(Num_total)
	for i in range(Num_total):
		line = []
		if round(X[i][0],2) >= subx0 and round(X[i][0],2) <= subx1:
			if round(X[i][1],2) >= subx2 and round(X[i][1],2) <= subx3:
				for j in range(2):
					line.append(X[i][j])
				sub_x.append(line)
				sub_y.append(y[i])
	print(sub_x,sub_y)
	tmp = 0
	Num_sub = len(sub_y)
	print(Num_sub)
	for i in range(Num_sub):
		tmp += sub_y[i]
	Num_subgroup_positive = tmp
	Num_subgroup_negative = len(sub_y) - Num_subgroup_positive

	complement_x = []
	complement_y = []
	for i in range(Num_total):
		line_com = []
		if round(X[i][0],2) < subx0 or round(X[i][0],2) > subx1:
			if round(X[i][1],2) < subx2 and round(X[i][1],2) > subx3:
				for j in range(2):
					line_com.append(X[i][j])
				complement_x.append(line)
				complement_y.append(y[i])
	print(complement_x, complement_y)
	tmp1 = 0
	Num_complement = len(complement_y)
	for i in range(Num_complement):
		tmp1 += complement_y[i]
	Num_complement_positive = tmp
	Num_complement_negative = len(sub_y) - Num_subgroup_positive

	p_subgroup = Num_subgroup_positive/Num_sub
	p_data = Num_total_positive/Num_total
	WARcc = (Num_sub/Num_total)*(p_subgroup-p_data)
	return WARcc

WARcc = subgroup(485.6, 540.1, 1.01, 1.72)
print(WARcc)







