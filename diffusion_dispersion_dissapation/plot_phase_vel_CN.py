# -*- coding: utf-8 -*-
"""
Created on Wed Oct 13 18:12:51 2021

@author: bucky
"""

import numpy as np
import matplotlib.pyplot as plt
plt.style.use('seaborn')

# plot the phase velocity as a function of k
x = np.linspace(0, 10, 100)
y = 2*np.arctan2(2*np.sin(x/100), 1)/(x/100)

plt.plot(x,y)
plt.ylabel('Phase Velocity')
plt.xlabel('Wavenumber k')
plt.title('Phase Velocity vs. Wavenumber for CN, hx=ht=1/100')