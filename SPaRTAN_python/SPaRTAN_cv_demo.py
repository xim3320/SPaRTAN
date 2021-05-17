import sys
import os
import random
import numpy as np
import pandas as pd
from sklearn.model_selection import KFold
from sklearn.preprocessing import normalize
from scipy.io import loadmat
from scipy import spatial
from scipy import stats
from SPaRTAN import SPaRTAN


def output_corr(file_corr, corr):
    with open(file_corr, "w") as outfile:
        # Any line starting with "#" will be ignored by numpy.loadtxt
        outfile.write("# Array shape: {0}\n".format(corr.shape))

        for a in range(0, len(corr)):
            spA = corr[a]
            outfile.write("# spectrumA: {}\n".format(spectrumAs[a]))
            for b in range(0, len(spA)):
                spB = spA[b]
                outfile.write("# spectrumB: {}\n".format(spectrumBs[b]))
                # save the labmda by rsL2 cv performace array
                np.savetxt(outfile, spB, fmt="%-7.4f")
                outfile.write("# \n")


print("Loading dataset pbmc5kdc.mat")
dataset = loadmat("../data/pbmc5kdc.mat")

D_ori = dataset["D"]
P_ori = dataset["Ppbmc5kdc"]
Y_ori = dataset["Ypbmc5kdc"]

# normalize the dataset
D = normalize(D_ori, norm="l2", axis=0)
Y = normalize(Y_ori, norm="l2", axis=0)
P = normalize(P_ori, norm="l2", axis=1)

# create the object of SPaRTAN
reg = SPaRTAN()

lamdas = [0.001, 0.01, 0.1, 0.2, 0.3]
rsL2s = [0, 0.001, 0.01]
spectrumAs = [1]
spectrumBs = [0.5, 0.6, 0.7]

lenlamdas = len(lamdas)
lenrsL2s = len(rsL2s)
lenspAs = len(spectrumAs)
lenspBs = len(spectrumBs)

# cross-validate to determine hyperparameters
fold = 5
corr_all_spearman = np.zeros((lenspAs, lenspBs, lenlamdas, lenrsL2s))
for a in range(0, lenspAs):
    for b in range(0, lenspBs):
        for l in range(0, lenlamdas):
            for r in range(0, lenrsL2s):
                print(
                    "cross validating spectrumA={}, spectrumB={}, lambda={}, rsL2={}".format(
                        spectrumAs[a], spectrumBs[b], lamdas[l], rsL2s[r]
                    )
                )
                sum_corr_spearman = 0

                kf = KFold(n_splits=fold)
                for train_index, test_index in kf.split(P_ori):

                    # split the data into train and test set
                    P_train, P_test = P_ori[train_index, :], P_ori[test_index, :]
                    Y_train, Y_test = Y_ori[:, train_index], Y_ori[:, test_index]

                    # normalize the train and test set
                    Y_train = normalize(Y_train, axis=0)
                    Y_test = normalize(Y_test, axis=0)

                    P_train = normalize(P_train, axis=1)
                    P_test = normalize(P_test, axis=1)

                    # train the model
                    reg.fit(
                        D,
                        P_train,
                        Y_train,
                        lamda=lamdas[l],
                        rsL2=rsL2s[r],
                        spectrumA=spectrumAs[a],
                        spectrumB=spectrumBs[b],
                    )

                    # get predicted value Y_pred  on P_test
                    Y_pred = reg.predict(P_test)

                    # get the correlation bewteen Y_pred and Y_test
                    corr_spearman = reg.get_corr(Y_pred, Y_test)

                    sum_corr_spearman = sum_corr_spearman + corr_spearman

                corr_all_spearman[a, b, l, r] = sum_corr_spearman / fold

if not os.path.exists("./output"):
    os.makedirs("./output")

outfile_corr = "./output/cvPerform.txt"
output_corr(outfile_corr, corr_all_spearman)

# retrive the best parameters
max_a, max_b, max_l, max_r = np.unravel_index(
    corr_all_spearman.argmax(), corr_all_spearman.shape
)

lamda_best = lamdas[max_l]
rsL2_best = rsL2s[max_r]
spectrumA_best = spectrumAs[max_a]
spectrumB_best = spectrumBs[max_b]

# re-train the model
reg.fit(D, P, Y, lamda_best, rsL2_best, spectrumA_best, spectrumB_best)

# retrieve W, projD, projP
W = reg.get_W()
projD = reg.get_projD()
projP = reg.get_projP()

outfile_W = "./output/W.csv"
outfile_projP = "./output/projP.csv"
outfile_projD = "./output/projD.csv"

with open(outfile_W, "w") as outfile:
    np.savetxt(outfile, W, fmt="%-7.4f")

with open(outfile_projP, "w") as outfile:
    np.savetxt(outfile, projP, fmt="%-7.4f")

with open(outfile_projD, "w") as outfile:
    np.savetxt(outfile, projD, fmt="%-7.4f")

print("Process finished successfully!")
