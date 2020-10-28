import pathlib
from SPaRTAN import SPaRTAN
from scipy.io import loadmat
import os
from random import sample
import numpy as np
import matplotlib.pyplot as plt


#load the data
# path = pathlib.Path(__file__).parent.absolute()
dataset = loadmat('../data/pbmc5kdc.mat')

D = dataset['D']
P = dataset['Ppbmc5kdc']
Y = dataset['Ypbmc5kdc']

#split into training and testing
nsamples = P.shape[0]    
indices = sample(range(nsamples),int(nsamples*0.8))
P_train = P[indices,:]
Y_train = Y[:,indices]
P_test = np.delete(P, indices, 0)
Y_test = np.delete(Y, indices, 1)

#create the object of AffReg
reg = SPaRTAN()

#train the model
reg.fit( D, P_train, Y_train, lamda = 0.001, rsL2 = 0, spectrumA = 1, spectrumB = 0.6 )

#predict test data
Y_pred = reg.predict(P_test)

#calculate the correlation between Y_pred with ground truth, and plot
corr = reg.get_corr(Y_pred, Y_test)
print("correlation = %.3f"%corr)

#retrive interaction between TF and protein (W)
W = reg.get_W()

# get projected protein expression on Y_train
projP = reg.get_projP()

# get projected TF activity on P_train
projD = reg.get_projD()
