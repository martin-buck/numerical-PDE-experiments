# -*- coding: utf-8 -*-
"""
Created on Mon Oct 11 00:12:16 2021

@author: bucky
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as ani
plt.style.use('seaborn')


def make_waves(k, L):
    # L is an integer that determines which PDE we want to visualize solutions
    # to
    if L == 1:
        # just make the wave speed 1
        alpha = k
        beta = np.array([0 for i in range(len(k))])
        filename = 'traveling_wave'
    elif L == 2:
        # ask for dissapation constant
        # D = float(input('Dissapation constant? \n'))
        D = .01
        alpha = k
        beta = -D*k**2
        filename = 'dissapation_wave'
    elif L == 3:
        # ask for dispersion constant
        # mu = float(input('Dispersion constant? \n'))
        mu = 1
        alpha = k+mu*k**3
        beta = np.array([0 for i in range(len(k))])
        filename = 'dispersion_wave'
        
    # make a traveling wave vector at various points in time
    window = 100
    T = 100
    X = 1000
    x = np.array([i/window for i in range(X)])
    U = np.zeros((T, X))
    for t in range(T):
        u = np.array([0 for i in range(window)])
        for j in range(len(k)):
            utemp = np.exp(beta[j]*t)*np.sin(k[j]*2*np.pi*x[t:100+t]-2*np.pi*alpha[j]*t/window)
            u = u + utemp
        U[t, t:100+t] = u
    
    # make gif for eyeball metric
    fig = plt.figure()
    def build_plot(i):
        plt.cla()
        plt.ylim(-2, 2)
        plt.plot(x, U[i,:])
        plt.xlabel('Interval [0,10]')
        plt.ylabel('Solution u(x)')
        plt.grid(True)
    
    animator = ani.FuncAnimation(fig, build_plot)
    animator.save(filename + '.gif')
    plt.show()



        