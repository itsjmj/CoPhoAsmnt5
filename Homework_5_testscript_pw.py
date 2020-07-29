'''Test script for Homework 5, Computational Photonics, SS 2020:  Scattering method.
'''


import numpy as np
import scipy.special
from Homework_5_function_headers import pw_expansion, mie_pw
from matplotlib import pyplot as plt


plt.rcParams.update({
        'figure.figsize': (12/2.54, 9/2.54),
        'figure.subplot.bottom': 0.15,
        'figure.subplot.left': 0.165,
        'figure.subplot.right': 0.90,
        'figure.subplot.top': 0.9,
        'axes.grid': False,
})

plt.close('all')

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# geometry parameters
max_x = ??  # width of computatinal domain [µm]
n0 = 1     # refractive index of the surrounding
nc = 1.5   # refractive index of the cylinder 
Rad = 1    # radius of the cylinder [µm]

# simulation parameters
lam = 1   # wavelength [µm]
Nx = ??  # number of grid points
m = ??    # number of the expansion orders retained


# %% run simulations and create representative figures of the results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# please add your code here
