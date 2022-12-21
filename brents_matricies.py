# -*- coding: utf-8 -*-
"""
Created on Mon Dec 19 20:42:38 2022

@author: jogib
"""

import numpy as np
from itertools import product
#%%
N=3
rows_g = product([0,1],repeat=N)
rows= [i for i in rows_g]
#%%
matricies_g = product(rows,repeat=N)
matricies = [i for i in matricies_g if np.linalg.det(i)%2==1]
#%%
def check_rows(arr_1,arr_2):
    return 1
def check_cols(arr1,arr2):
    return 1
def check_rc(arr1,arrr2):
    return check_rows(arr1,arr2) | check_cols(arr1,arr2)

