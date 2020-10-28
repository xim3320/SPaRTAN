# Run SPaRTAN in Python

### Introduction



### Prerequesities
The code runs on Python 3, and following packages are used:
pandas, pickle, numpy, random, os, sklearn, scipy, Cython, 

### Cython extension built
In order to improve the running time, we convert two computationally intensive python modules into Cython modules, which calls some of C codes as well. We have uploaded Cython extensions built on a Windows and a Linux system. If they are not compatible with your platform, then you need to build them on site.


- build cythLeastR extention modue 
    
	1. Goto folder "cythKrnPlus built", then execute command:
	
	```sh
	python setup.py build_ext --inplace
	```
	
	It will generate a .so (Mac or Linux), or a .pyd (Windows) file. 
	 
	2. Copy the .so(.pyd) file into SPaRTAN folder
    

- build cythKrnPlus extention modue 
    
	1. Goto folder "cythKrnPlus built", then execute command:
	
	```sh
	python setup.py build_ext --inplace
	```
	It will generate a.so (Mac or Linux), or a.pyd (Windows) file. 
     
	2. Copy the .so(.pyd) file into SPaRTAN folder

### Load the data
We load the data into dataset object
```sh
dataset = loadmat("../data/pbmc5kdc.mat")
D = dataset['D']
P = dataset['Ppbmc5kdc']
Y = dataset['Ypbmc5kdc']

```
### Cross validation
SPaRTan model has 4 parameters pectrumA, spectrumB, rsL2 and lambda. Their values are determined by the user input data D, P, and Y. We use cross-validation to tune the parameters.
#### Split samples into training and testing set
First we need to split the samples of P and Y matrix into training and testing set:

```sh
from sklearn.model_selection import KFold

fold = 5
kf = KFold(n_splits = fold)
for train_index, test_index in kf.split(PP):
    P_train, P_test = PP[train_index,:], PP[test_index,:]
    Y_train, Y_test = YY[:,train_index], YY[:,test_index]
     ....

```
#### nomralize data
To improve the performance, we apply normalization to the D matrix and traning, testing set of P and Y matrices

```sh
from sklearn.preprocessing import normalize
D = normalize(D,axis=0)  #normalize by column
P_train = normalize(D, axis=1) #normalize by row
Y_train = normalize(Y, axis=0) #normalize by column
```
#### Create object of the model class
SPaRTAN model is integrated into a python class. We first need to initilize the object of class with

```sh
reg = SPaRTAN()
```
#### train the model with possible combination of parameters of each spliting fold
For each fold, we train the model with D, P_train, Y_train and test value of 4 parameters by calling function fit of the class

```sh
 reg.fit( D, P_train, Y_train, lamda = lamda_choice, rsL2 = rsL2_test, spectrumA = spectrumA_test, pectrumB = spectrumB_test)
```

#### predict the gene expression
After train the model with fit function, we then predict Y_test_predict with test set of P(P_test) and extract the correlation between predicted Y and observed Y 
```sh
Y_test_pred = reg.predict(P_test)
```
#### get the correlation of predicted gene expression and observed expression
For each fold, and each combination of parameters, calculate the correlation between predicted Y_test_pred and observed Y_test. The best parameters will be dtermined based on the average correlatin value of all folds on test set.
```sh
corr = reg.corr(Y_test_pred, Y_test)
```
#### complete implementatioin of cross-validation
The complete implementation of cross-validation is
```sh
D = normalize(D, axis=0)

#create the object of SPaRTAN
reg = SPaRTAN()

lamdas = [0.001, 0.01]#, 0.1, 0.2, 0.3 ]
rsL2s = [0, 0.001]#, 0.01]
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
                         
                    P_train, P_test = P[train_index,:], P[test_index,:]
                    Y_train, Y_test = Y[:,train_index], Y[:,test_index]

                    Y_train = normalize(Y_train, axis=0)
                    Y_test = normalize(Y_test, axis=0)
						
                    P_train = normalize(P_train, axis=1)
                    P_test = normalize(P_test, axis=1)
	
                    reg.fit( D, P_train, Y_train, lamda = lamdas[l], rsL2 = rsL2s[r], spectrumA = spectrumAs[a], spectrumB = spectrumBs[b]  )

                    Y_pred = reg.predict(P_test)
	
                    corr_spearman = reg.get_corr(Y_pred, Y_test)
                    
                    sum_corr_spearman = sum_corr_spearman + corr_spearman

                corr_all_spearman[a, b, l, r] = sum_corr_spearman/fold
```
### Extract best parameters
After the cross-validation, we got the correlation values of all possible combination of the parameters. Then we pick the best parameter combination with which yield biggest correlation

```sh
max_a,max_b,max_l,max_r = np.unravel_index(corr_all_spearman.argmax(), corr_all_spearman.shape)
lamda_best = lamdas[max_l]
rsL2_best = rsL2s[max_r]
spectrumA_best = spectrumAs[max_a]
spectrumB_best = spectrumBs[max_b]
```
### Train and Predict Features with best parameters
#### train the model with whole dataset
Now we can use the best parameters to train the model with the whole datasetto predict the features in interest
```sh
#normalize P and Y, D has been normalized previously
Y = normalize(Y, axis=0)
P = normalize(P, axis=1)
reg.fit( D, P, Y, lamda = lamda_best, rsL2 = rsL2_best, spectrumA = spectrumA_best, spectrumB = spectrumB_best  )
```

#### Extract TF and protein interaction
Get the interaction W matrix based on D, P, Y
```sh
W = reg.get_W()
```
#### Extract projected protein expression
Get the projected P based on D, W and Y
```sh
proj_P = reg.get_projP()
```
#### Extract projected TF activities
Get the projected TF by sample matrix based on W and P
```sh
proj_D = reg.get_projD()
```
