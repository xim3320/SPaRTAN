# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 13:31:46 2020

@author: XIM33
"""

import numpy as np
import cythLeastR
from scipy.stats import norm
#X = np.random.rand(10,4)
#y = np.random.rand(10,1)


np.random.seed(1)
X = norm.ppf(np.random.rand(4,10))
X = X.T
y = norm.ppf(np.random.rand(1,10))
y = y.T


opts = dict()   
opts['rsL2']=0
#opts['ind']=np.asarray([0,10])


pylambda = 0.1

beta, b = cythLeastR.LeastR(X,y,pylambda, opts)
print(beta)
print(b)
