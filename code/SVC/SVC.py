# SVM classfier with linear kernel
# read train.dat

import numpy as np
import time
from sklearn.model_selection import GridSearchCV
from sklearn.svm import LinearSVC
from sklearn.svm import SVC
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import  roc_curve, roc_auc_score
start=time.time()

inp=open('train.txt','r').readlines()
X=[]
y_train=np.array([])
for i in range(1,len(inp)):
	line=str(inp[i]).split()
	feat=[]	
	for j in range(2,len(line)):  # line[0]: name; line[1]: property; line[2] first feature
		feat.append(float(line[j]))
	X.append(feat)
	y_train=np.append(y_train,float(line[1]))
X_train=np.array(X)

# Linear kernel
clf = GridSearchCV(SVC(kernel='linear',tol=1e-4,max_iter=-1),
         cv=10, n_jobs=-1, verbose=1,
         param_grid={"C": np.logspace(1, 5, 10, base=10)},
         return_train_score=True,refit=True)

file_out = open("out.txt", "w")
print("training:", file=file_out)
clf.fit(X_train, y_train)
print("params: ", file=file_out)
print(clf.cv_results_['params'], file=file_out)
print("mean_test_score: ", file=file_out)
print(clf.cv_results_['mean_test_score'], file=file_out)
print("mean_train_score: ", file=file_out)
print(clf.cv_results_['mean_train_score'], file=file_out)
print("best CV score: ", file=file_out)
print(clf.best_score_, file=file_out)
print("best_params: ", file=file_out)
print(clf.best_params_, file=file_out)
print("best_estimator: ", file=file_out)
print(clf.best_estimator_, file=file_out)
print("coef: (sum(c_i*x_i)+c_0=0) ", file=file_out)
print(clf.best_estimator_.coef_, file=file_out)
print("intercept: ", file=file_out)
print(clf.best_estimator_.intercept_, file=file_out)
print("support: ", file=file_out)
print(clf.best_estimator_.support_, file=file_out)
print(clf.best_estimator_.support_vectors_, file=file_out)

def train_output(X_train,y_train):
    y_pred=clf.predict(X_train)
    er=0
    fail_index=[]
    for i in range(len(y_train)):
          print('Y_ture, Y_pred, data_point: ',np.around(y_train[i],decimals=3),np.around(y_pred[i],decimals=3),i, file=file_out)
          if np.sign(y_train[i]) != np.sign(y_pred[i]):
              er+=1
              fail_index.append(i)
    print('Wrong prediction: Data point ',fail_index, file=file_out)
    print("Total number of misclassified data: ",er, file=file_out)
    fpr, tpr, _ = roc_curve(y_train, y_pred)
    AUC = roc_auc_score(y_train, y_pred)
    ROC_file = open('tpr&fpr.txt', 'a')
    print('fpr: ', [i for i in fpr], file=ROC_file)
    print('tpr: ', [j for j in tpr], file=ROC_file)
    print('AUC = ', AUC, file=ROC_file)
    ROC_file.close()

train_output(X_train,y_train)
print("training_score (accuracy): ",clf.best_estimator_.score(X_train,y_train), file=file_out)

end=time.time()
print("wall-clock time (seconds):",end-start, file=file_out)
file_out.close()