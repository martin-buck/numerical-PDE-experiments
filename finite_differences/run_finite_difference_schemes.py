# -*- coding: utf-8 -*-
"""
Created on Thu Sep 23 12:56:58 2021

@author: bucky

Numerical PDE schemes and analysis
"""

# import
from finite_difference_schemes import two_level_schemes
import matplotlib.pyplot as plt
import matplotlib.animation as ani
from error import error
import numpy as np
import math
plt.style.use('seaborn')


# ask for input on IC, spacing, and scheme
ic_in = int(input('Linear, step, sin, tent, or wave packet for IC? (1,2,3,4,5) \n'))
x = int(input('Enter integer for resolution of space (e.g. 100): \n'))
t = int(input('Enter integer for resolution of time (e.g. 100): \n'))
scheme = input('Type ETCS, ETBS, ITCS, or CN: \n')

if ic_in == 1:
    IC = [i/x for i in range(x)]
    ic_str = 'Linear'
elif ic_in == 2:
    if x % 2 == 0:
        IC1 = np.zeros(math.floor(x/2))
    else:
        IC1 = np.zeros(math.floor(x/2)+1)
    IC2 = np.array([1/2 for i in range(math.floor(x/2))])
    IC = np.concatenate((IC1, IC2), axis=0)
    ic_str = 'Step'
elif ic_in == 3:
    IC = [np.sin(2*np.pi*i/x) for i in range(x)]
    ic_str = 'Sin'
elif ic_in == 4:
    if x % 2 == 0:
        IC1 = np.array([2*i/x for i in range(math.floor(x/2))])
    else:
        IC1 = np.array([2*i/x for i in range(math.floor(x/2))+1])
    IC2 = np.array([(-2*(i+x/2)/x)+2 for i in range(math.floor(x/2))])
    IC = np.concatenate((IC1, IC2), axis=0)
    ic_str = 'Tent'
elif ic_in == 5:
    IC = [np.sin(2*np.pi*i/x) + np.sin(100*np.pi*i/x) for i in range(x)]
    ic_str = 'Wave Packet'

# only tested with homogeneous BC at x=0
BC = [0 for i in range(t)]

# call solver
S = two_level_schemes(x, t, IC, BC, scheme)

# make gif for eyeball metric
fig = plt.figure()
def build_plot(i):
    plt.cla()
    if i < t:
        plt.plot([i/x for i in range(x)], S[i,:])
        plt.ylim(min(IC), max(IC))
        plt.xlabel('Interval [0,1]')
        plt.ylabel('Solution u(x)')
        plt.grid(True)
        plt.title('One-way wave equation, IC=' + ic_str + ', ' + scheme + ', ht=hx=1/100')

filename = ic_str + '_' + scheme
animator = ani.FuncAnimation(fig, build_plot)
animator.save(filename + '.gif')

# print out error for different time and grid spacings for ETBS
r = input('Make a plot for ETBS for different grid spacings? (0, 1) \n')

if r == '1':
    # try different ratios of h_t to h_x
    spacings_t = np.array([90,91,92,93,94,95,96,97,98,99,100])
    spacing_x = 100
    ratios = spacings_t/spacing_x
    errors = np.zeros(len(spacings_t))
    
    # call solver for these grid spacings
    for i in range(len(spacings_t)):
        t = spacings_t[i]
        x = 100
        BC = [0 for i in range(t)]
                
        S = two_level_schemes(x, t, IC, BC, 'ETBS')
        errors[i] = error(IC, BC, S)
    
    # plot error vs. spacing ratio
    fig = plt.figure()
    plt.plot(ratios, errors, marker='o')
    plt.xlabel('Grid spacing ratio (ht/hx)')
    plt.ylabel('Maximum l2 error over time')
    plt.title('Maximum l2 Error For CFL Condition, IC=' + ic_str)
    plt.yscale('log')
    plt.ylim(10**-5, 10**5)
    plt.grid(True)
    plt.savefig('l2_error_ETBS_CFL_' + ic_str)