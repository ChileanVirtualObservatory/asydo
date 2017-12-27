from sklearn import svm
import numpy as np

trainPerc=0.75

X=np.nan_to_num(np.load('../dataset-ML/exp2.npy'))
siz=len(X)
half=int(siz/2)
X0=X[0:half]
X1=X[half:siz]
print "len(X0)=" + str(len(X0))
print "len(X1)=" + str(len(X1))
X0train=X0[0:int(trainPerc*len(X0))]
Y0train=np.zeros(len(X0train))
X1train=X1[0:int(trainPerc*len(X1))]
Y1train=np.ones(len(X1train))
print "len(X0train)=" + str(len(X0train))
print "len(X1train)=" + str(len(X1train))
X0test=X0[int(trainPerc*len(X0)):len(X0)]
Y0test=np.zeros(len(X0test))
X1test=X1[int(trainPerc*len(X1)):len(X1)]
Y1test=np.ones(len(X1test))
print "len(X0test)=" + str(len(X0test))
print "len(X1test)=" + str(len(X1test))
Xtrain=np.vstack((X0train,X1train))
Ytrain=np.hstack((Y0train,Y1train))
print "len(Xtrain)=" + str(len(Xtrain))
Xtest=np.vstack((X0test,X1test))
Ytest=np.hstack((Y0test,Y1test))
print "len(Xtest)=" + str(len(Xtest))
clf = svm.SVC()
clf.fit(Xtrain, Ytrain)  
res=clf.predict(Xtest)
vect=np.abs(res - Ytest)
print str((len(Ytest) - np.sum(vect))/len(Ytest)*100) + "%"
