# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 11:41:07 2021

@author: bucky
"""
import numpy as np

# create a function which takes two grid spacings, and a vector of initial and boundary conditions
def two_level_schemes(x, t, IC, BC, method):
    # check to make sure the size of IC and BC are consistent with hx and ht
    m = len(IC)
    n = len(BC)
    if n != t or m != x:
        raise Exception('Length of arrays must match grid spacings')
    
    if IC[0] != BC[0]:
        raise Exception('Arrays must match in their first entry')
    
    ht = 1/t
    hx = 1/x
    # create the matrix Q which will solve the transport equation at each time
    if method == 'ETCS':
        diag1m = [ht/(2*hx) for i in range(m-1)]
        diag = [1 for i in range(m)]
        diag1p = [-ht/(2*hx) for i in range(m-1)]
        Q = np.diag(diag1m, -1) + np.diag(diag, 0) + np.diag(diag1p, 1)
        P = np.identity(m)
    elif method == 'ETBS':
        diag1m = [ht/hx for i in range(m-1)]
        diag = [1-(ht/hx) for i in range(m)]
        diag1p = [0 for i in range(m-1)]
        Q = np.diag(diag1m, -1) + np.diag(diag, 0) + np.diag(diag1p, 1) 
        P = np.identity(m)
    elif method == 'CN':
        q_diag1m = [ht/(4*hx) for i in range(m-2)]
        q_diag = [1 for i in range(m-1)]
        q_diag1p = [-ht/(4*hx) for i in range(m-2)]
        Q = np.diag(q_diag1m, -1) + np.diag(q_diag, 0) + np.diag(q_diag1p, 1) 
        p_diag1m = [-ht/(4*hx) for i in range(m-2)]
        p_diag = [1 for i in range(m-1)]
        p_diag1p = [ht/(4*hx) for i in range(m-2)]
        P = np.diag(p_diag1m, -1) + np.diag(p_diag, 0) + np.diag(p_diag1p, 1)
    elif method == 'ITCS':
        diag1m = [-ht/(2*hx) for i in range(m-1)]
        diag = [1 for i in range(m)]
        diag1p = [ht/(2*hx) for i in range(m-1)]
        Q = np.identity(m)
        P = np.diag(diag1m, -1) + np.diag(diag, 0) + np.diag(diag1p, 1)

    S = np.zeros((t,x))
    S[:, 0] = BC
    S[0, :] = IC
    
    # fill in S that will store the solution on the grid for each time
    for i in range(t-1):
        u_i = S[i,:]
        if method == 'ETCS' or method == 'ETBS':
            u_i1 = np.matmul(Q, u_i)
            # Impose BC on right side of grid
            S[i+1, :] = u_i1
        else:
            v_i = u_i[1:]
            Qv_i = np.matmul(Q, v_i)
            Qv_i[0] = Qv_i[0] + (ht/(4*hx))*BC[i+1] + (ht/(4*hx))*BC[i]
            Qv_i[-1] = Qv_i[-1] - (ht/(4*hx))*S[i-1, -1] - (ht/(4*hx))*S[i, -1]
            u_i1 = np.linalg.solve(P, Qv_i)
            S[i+1, 1:] = u_i1
    return S