# Run SPaRTAN in Python

### Introduction
This is the Python/Cython implementation of the SPaRTAN. In order to improve the running time performance, we convert some computationally intensive python modules into Cython modules.  All functionalities are integrated into the class SPaRTAN . Please check SPaRTAN_demo.py file for the basic use of this class. This tutorial focuses on how to apply SPaRTAN in a real project situation.

### Prerequisites
The code runs on Python 3, and the following packages are used:

pandas, pickle, numpy, random, os, sklearn, scipy, Cython, 

### Cython extension built
There are two Cython extension modules needed for running SPARTAN. We have built the extensions under Windows(.pyx files) and Linux/Mac (.so files) system. You can download ones based on your operating system. If they are not compatible with your platform, then you need to build Cython extension on site. The followings are the instruction on how to build the Cython extension

- build cythLeastR extension module 
    
	1. Go to folder "cythLeastR built", then execute the command:
	
	```sh
	python setup.py build_ext --inplace
	```
	
	 It will generate a .so (Mac or Linux), or a .pyd (Windows) file. 
	 
	2. Copy the .so(.pyd) file into "SPaRTAN_python" folder
    

- build cythKrnPlus extension module 
    
	1. Go to folder "cythKrnPlus built", then execute the command:
	
	```sh
	python setup.py build_ext --inplace
	```
	 It will generate a.so (Mac or Linux), or a.pyd (Windows) file. 
     
	2. Copy the .so(.pyd) file into "SPaRTAN_python" folder

### Cross-validation
SPaRTAN model has 4 parameters spectrumA, spectrumB, rsL2 and lambda. Their values are determined by the user input data D, P, and Y. We use cross-validation to determine the best values of those parameters. Here we explain step by step

**Load the data**

In the example, we load a Matlab format dataset
```sh
dataset = loadmat("../data/pbmc5kdc.mat")
D = dataset['D']
P = dataset['Ppbmc5kdc']
Y = dataset['Ypbmc5kdc']
```

**Split the samples of P and Y matrix into training and testing set:**

```sh
from sklearn.model_selection import KFold

fold = 5
kf = KFold(n_splits = fold)
for train_index, test_index in kf.split(PP):
    P_train, P_test = P[train_index,:], P[test_index,:]
    Y_train, Y_test = Y[:,train_index], Y[:,test_index]
     ....

```
**Apply normalization to the dataset**

```sh
from sklearn.preprocessing import normalize
D = normalize(D,axis=0)  #normalize by column
P_train = normalize(P_train, axis=1) #normalize by row
Y_train = normalize(Y_train, axis=0) #normalize by column
```
**Train the model with cross-validation**

SPaRTAN model is integrated into a python class. We first need to initialize e the object of the class

```sh
reg = SPaRTAN()
```
For each fold, we train the model with D, P_train, Y_train, and test values of 4 parameters by calling the function "fit" 

```sh
 reg.fit( D, P_train, Y_train, lamda = lamda_choice, rsL2 = rsL2_test, spectrumA = spectrumA_test, pectrumB = spectrumB_test)
```

After train the model with fit function, we then predict Y_test_predict with a test set of P(P_test) and extract the correlation between predicted Y and observed Y 
```sh
Y_test_pred = reg.predict(P_test)
```
For each fold and each combination of parameters, calculate the correlation between predicted Y_test_pred and observed Y_test. The best parameters will be determined d based on the average correlation  value of all folds on the test set.
```sh
corr = reg.corr(Y_test_pred, Y_test)
```
**Extract best parameters**

After the cross-validation, we got the correlation values of all possible combination of the parameters. Then we pick the best parameter combinations with which yield the biggest correlation

```sh
max_a,max_b,max_l,max_r = np.unravel_index(corr_all_spearman.argmax(), corr_all_spearman.shape)
lamda_best = lamdas[max_l]
rsL2_best = rsL2s[max_r]
spectrumA_best = spectrumAs[max_a]
spectrumB_best = spectrumBs[max_b]
```
**Complete implementation  of cross-validation**
```sh
D = normalize(D, axis=0)

#create the object of SPaRTAN
reg = SPaRTAN()

lamdas = [0.001, 0.01, 0.1, 0.2, 0.3 ]
rsL2s = [0, 0.001, 0.01]
spectrumAs = [1]
spectrumBs = [0.5, 0.6, 0.7 ]

lenlamdas = len(lamdas)
lenrsL2s = len(rsL2s)
lenspAs = len(spectrumAs)
lenspBs = len(spectrumBs)

# cross-validate to determine hyperparameters
fold = 5
corr_all_spearman = np.zeros((lenspAs, lenspBs, lenlamdas, lenrsL2s)   ) 
for a in range(0, lenspAs):
    for b in range(0, lenspBs):
        for l in range(0, lenlamdas):
            for r in range(0, lenrsL2s):
                print("cross validating spectrumA={}, spectrumB={}, lambda={}, rsL2={}".format(spectrumAs[a], spectrumBs[b], lamdas[l], rsL2s[r]))
                sum_corr_spearman = 0

                kf = KFold(n_splits = fold)
                for train_index, test_index in kf.split(P):
                    
		    # split dataset into train and test set
                    P_train, P_test = P[train_index,:], P[test_index,:]
                    Y_train, Y_test = Y[:,train_index], Y[:,test_index]
                    
		    # Apply normalization on train and test set
                    Y_train = normalize(Y_train, axis=0)
                    Y_test = normalize(Y_test, axis=0)
						
                    P_train = normalize(P_train, axis=1)
                    P_test = normalize(P_test, axis=1)
                    
		    # train the model
                    reg.fit( D, P_train, Y_train, lamda = lamdas[l], rsL2 = rsL2s[r], spectrumA = spectrumAs[a], spectrumB = spectrumBs[b]  )

                    Y_pred = reg.predict(P_test)
	
                    corr_spearman = reg.get_corr(Y_pred, Y_test)
                    
                    sum_corr_spearman = sum_corr_spearman + corr_spearman

                corr_all_spearman[a, b, l, r] = sum_corr_spearman/fold

#extract best parameters
max_a,max_b,max_l,max_r = np.unravel_index(corr_all_spearman.argmax(), corr_all_spearman.shape)
lamda_best = lamdas[max_l]
rsL2_best = rsL2s[max_r]
spectrumA_best = spectrumAs[max_a]
spectrumB_best = spectrumBs[max_b]
```

### Train the model and get the projected data matrices
**Train the model again with the whole dataset and best parameters**

Now we can use the best parameters to train the model with the whole dataset to predict the features of the interest
```sh
#normalize P and Y, D has been normalized previously
Y = normalize(Y, axis=0)
P = normalize(P, axis=1)
reg.fit( D, P, Y, lamda = lamda_best, rsL2 = rsL2_best, spectrumA = spectrumA_best, spectrumB = spectrumB_best  )
```

**Extract TF and protein interaction**

Get the interaction W matrix based on D, P, Y
```sh
W = reg.get_W()
```
**Extract projected protein expression**

Get the projected P based on D, W and Y
```sh
proj_P = reg.get_projP()
```
**Extract projected TF activities**

Get the projected TF by sample matrix based on W and P
```sh
proj_D = reg.get_projD()
```
