# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 00:38:46 2021

@author: bucky
"""
import numpy as np
from qualitative_analysis import make_waves

de = int(input('Which PDE: L1, L2, or L3? Input 1, 2, or 3 \n'))

if de == 1:
    k = np.array([1])
    L = 1
elif de == 2:
    k = np.array([1])
    L = 2
elif de == 3:
    k = np.array([1, 2])
    L = 3

make_waves(k, L)