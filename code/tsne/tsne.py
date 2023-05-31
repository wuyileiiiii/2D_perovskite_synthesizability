import numpy as np
from sklearn.manifold import TSNE

np.set_printoptions(threshold= 1000000)

inp=open('train.txt','r').readlines()
X=[]
y_train=np.array([])
for i in range(1,len(inp)):   # from the second line to the last line
	line=str(inp[i]).split()
	feat=[]
	for j in range(2,len(line)):  # line[0]: name; line[1]: property; line[2] first feature
		feat.append(float(line[j]))
	X.append(feat)
	y_train=np.append(y_train,float(line[1]))
X_train=np.array(X)

file_out = open("result1.txt",'w')

# 4个3维的数据

ts = TSNE(n_components=2)
# 训练模型
ts.fit_transform(X_train)
# 打印结果
print(ts.embedding_, file = file_out)




