'''Test script for Homework 5, Computational Photonics, SS 2020:  Scattering method.
'''


import numpy as np
import scipy.special
import bluered_dark
from Homework_5_function_headers import mie_pw, mie_gauss
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
max_x = 5  # width of computatinal domain [µm]
waist = 3  # waist of the Gaussian beam [µm]
n0 = 1     # refractive index of the surrounding
nc = 1.5   # refractive index of the cylinder 
Rad = 3    # radius of the cylinder [µm]

# simulation parameters
lam = 1   # wavelength [µm]
Nx = 200  # number of grid points
m = 10    # number of the expansion orders retained

x = np.linspace(-max_x,max_x,Nx)
y = np.linspace(-max_x,max_x,Nx)
X,Y = np.meshgrid(x,y)

# %% run simulations and create representative figures of the results %%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#%%
#Figures to test dky and order of gaussian incident beam
pol = 'TE'
u = mie_gauss(lam, waist, m, n0, nc, Rad, pol, max_x, Nx)

#%%
pol = 'TE'
if pol == 'TE':
    u_field = r'$E_z$'
elif pol == 'TM':
    u_field = r'$E_z$'


from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.ticker as tk
def colorbar(mappable):
    ax = mappable.axes
    fig = ax.figure
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.07)
    fmt = tk.ScalarFormatter(useOffset=True)
    fmt.set_scientific(1)
    return fig.colorbar(mappable, format=fmt,cax=cax)

#%%
plt.rcParams.update({'font.size': 10})
f1, ax = plt.subplots(dpi=130)
ax.axis('equal')
cm = ax.pcolormesh(X,Y,u.real,cmap = 'rainbow')
cb = colorbar(cm)
ax.set_xlabel(r'x [$\mu$m]')
ax.set_ylabel(r'y [$\mu$m]')
ax.set_title(r'm = 10, $\Delta k_y = \rho/2$')
cb.set_label('Re['+u_field+']'+' [A.u.]')
f1.tight_layout()
f1.savefig('ptest_m10_r2.png',bbox_inches='tight')
plt.show()
#%%
rdiv_range = np.rint(np.linspace(3,90,3))
print(rdiv_range)
#%%
#Figures showing Gaussian scattering over radius range
#rdiv_range = np.rint(np.linspace(3,95,4))
rdiv_range = np.array([2.5,50,99])
u_TE = []
u_TM = []
for rdiv in rdiv_range:
    ute = mie_gauss(lam, waist, m, n0, nc, lam/rdiv, 'TE', max_x, Nx)
    utm = mie_gauss(lam, waist, m, n0, nc, lam/rdiv, 'TM', max_x, Nx)
    u_TE.append(ute)
    u_TM.append(utm)
    
#%%

plt.rcParams.update({'font.size': 8})
f2, axs = plt.subplots(3,2,dpi=130)
axsf = axs.flatten()
for ax in axsf:
    ax.axis('equal')

for i in range(len(u_TE)):
    uTE = u_TE[i]
    uTM = u_TM[i]
    cm1 = axs[i,0].pcolormesh(X,Y,uTE.real,cmap = 'rainbow')
    axs[i,0].set_xlabel(r'x [$\mu$m]')
    axs[i,0].set_ylabel(r'y [$\mu$m]')
    axs[i,0].set_title(r'Rad = $\lambda$'+'/{}, TE'.format(rdiv_range[i]))
    axs[i,0].add_artist(plt.Circle((0, 0), lam/rdiv_range[i], facecolor='none', edgecolor='k',lw=0.3))
    cb1 = colorbar(cm1)
    cb1.set_label(r'Re[$E_z$]'+' [A.u.]')
    
    cm2 = axs[i,1].pcolormesh(X,Y,uTM.real,cmap = 'rainbow')
    axs[i,1].set_xlabel(r'x [$\mu$m]')
    axs[i,1].set_ylabel(r'y [$\mu$m]')
    axs[i,1].set_title(r'Rad = $\lambda$'+'/{}, TM'.format(rdiv_range[i]))
    axs[i,1].add_artist(plt.Circle((0, 0), lam/rdiv_range[i], facecolor='none', edgecolor='k',lw=0.3))
    cb2 = colorbar(cm2)
    cb2.set_label(r'Re[$H_z$]'+' [A.u.]')
    
    
    
f2.tight_layout()
#f2.savefig('re_gauss_vary_rad.png',bbox_inches='tight')
plt.show()

#%%
f3, axs = plt.subplots(3,2,dpi=130)
axsf = axs.flatten()
for ax in axsf:
    ax.axis('equal')

for i in range(len(u_TE)):
    uTE = u_TE[i]
    uTM = u_TM[i]
    cm1 = axs[i,0].pcolormesh(X,Y,np.abs(uTE)**2,cmap = 'magma')
    axs[i,0].set_xlabel(r'x [$\mu$m]')
    axs[i,0].set_ylabel(r'y [$\mu$m]')
    axs[i,0].set_title(r'Rad = $\lambda$'+'/{}, TE'.format(rdiv_range[i]))
    axs[i,0].add_artist(plt.Circle((0, 0), lam/rdiv_range[i], facecolor='none', edgecolor='k',lw=0.3))
    cb1 = colorbar(cm1)
    cb1.set_label(r'$|E_z|^2$'+' [A.u.]')
    
    cm2 = axs[i,1].pcolormesh(X,Y,np.abs(uTM)**2,cmap = 'magma')
    axs[i,1].set_xlabel(r'x [$\mu$m]')
    axs[i,1].set_ylabel(r'y [$\mu$m]')
    axs[i,1].set_title(r'Rad = $\lambda$'+'/{}, TM'.format(rdiv_range[i]))
    axs[i,1].add_artist(plt.Circle((0, 0), lam/rdiv_range[i], facecolor='none', edgecolor='k',lw=0.3))
    cb2 = colorbar(cm2)
    cb2.set_label(r'$|H_z|^2$'+' [A.u.]')
    
    
    
f3.tight_layout()
#f3.savefig('abs_gauss_vary_radpng',bbox_inches='tight')
plt.show()
    

