# -*- coding: utf-8 -*-
"""
Created on Wed Oct  6 13:01:08 2021

@author: bucky
"""

from collections import deque
import numpy as np

# calculate the l2 norm
def error(IC, BC, S):
    # for the one way wave equation the solution should just shift the IC at
    # each time step
    max_err = 0
    h_x = 1/len(IC)
    # compare true solution with calculated solution at each time step
    for i in range(1, len(BC)):
        # calculate the true solution by shifting the IC
        sol_i = deque(IC)
        sol_i.rotate(i)
        for j in range(0, i):
            sol_i[j] = 0
        
        # now compare sol_i with the corresponding row in S in l2 norm
        curr_err = np.sqrt(h_x*sum((np.array(sol_i)-S[i,:])**2))
        
        # assign max_err
        if max_err <= curr_err:
            max_err = curr_err
            
    return max_err
    