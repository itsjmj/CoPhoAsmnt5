'''Test script for Homework 5, Computational Photonics, SS 2020:  Scattering method.
'''


import numpy as np
import scipy.special
from Homework_5_function_headers import pw_expansion, mie_pw, mie_gauss
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
max_x = 20  # width of computatinal domain [µm]
n0 = 1     # refractive index of the surrounding
nc = 1.5#np.sqrt(13)#1.5   # refractive index of the cylinder 
Rad = 1/4    # radius of the cylinder [µm]

# simulation parameters
lam = 1 #5/0.5   # wavelength [µm]
Nx = 200  # number of grid points
m = 10    # number of the expansion orders retained


# %% run simulations and create representative figures of the results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
x = np.linspace(-max_x,max_x,Nx)
y = np.linspace(-max_x,max_x,Nx)
X,Y = np.meshgrid(x,y)
R = np.sqrt(X**2 + Y**2)
# ui, ui_an = pw_expansion(lam, m, n0, max_x, Nx)
waist = 3
u = mie_gauss(lam, waist, m, n0, nc, Rad, 'TM', max_x, Nx)

#%%
f = plt.figure()
plt.imshow(np.abs(u))
plt.show()

# f = plt.figure()
# plt.imshow(u.real)
# plt.show()

# u = mie_pw(lam, m, n0, nc, Rad, 'TE', max_x, Nx)
# f = plt.figure()
# plt.imshow(u.real)
# plt.show()

# f = plt.figure()
# plt.imshow(np.abs(u)**2)
# plt.show()
# #%%
# ut = np.zeros((Nx,Nx),dtype=complex)
# #ut[R<=Rad] = u[R<=Rad]
# ut[R>Rad] = ui_an[R>Rad] + u[R>Rad]
# f = plt.figure()
# plt.imshow(np.angle(ut))
# plt.show()

# f = plt.figure()
# plt.imshow(np.angle(ui))
# plt.show()
# please add your code here
