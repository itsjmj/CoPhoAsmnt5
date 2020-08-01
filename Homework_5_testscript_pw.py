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
max_x = 5 # width of computatinal domain [µm]
n0 = 1     # refractive index of the surrounding
nc = 1.5   # refractive index of the cylinder 
Rad = 1    # radius of the cylinder [µm]

# simulation parameters
lam = 1   # wavelength [µm]
Nx = 500  # number of grid points
m = 20   # number of the expansion orders retained


# %% run simulations and create representative figures of the results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

pol='TM' #polarization
u=mie_pw(lam, m, n0, nc, Rad, pol, max_x, Nx)


if pol=='TE':
    tit='$E_z$'
elif pol=='TM':
    tit='$H_z$'
else:
    raise ValueError('Please choose TE or TM mode.')

# plot field distribution
x = np.linspace(-max_x,max_x,Nx)
plt.figure(figsize=[10,7])
plt.pcolormesh(x,x,np.real(u),cmap='rainbow')
plt.colorbar()
plt.title('Scattered field disribution {}'.format(tit))
plt.xlabel('x/$\mu$m')
plt.ylabel('y/$\mu$m')
plt.show()

#plot field intensity
plt.figure(figsize=[10,7])
plt.pcolormesh(x,x,np.abs(u)**2,cmap='rainbow')
plt.colorbar()
plt.title('Scattered field intensity {}'.format(tit))
plt.xlabel('x/$\mu$m')
plt.ylabel('y/$\mu$m')
plt.show()
