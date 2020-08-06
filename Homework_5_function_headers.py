'''Homework 5, Computational Photonics, SS 2020:  Scattering method.
'''

import numpy as np
import scipy.special
from scipy.special import jv
from scipy.special import jvp
from scipy.special import hankel1 as hv
from scipy.special import h1vp as hvp


def pw_expansion(lam, m, n0, max_x, Nx):
    '''Calculates the amplitude distribution of a plane wave expanded in cylindrical functions.
       Computational domain has dimensions x = [-max_x, max_x], y = [-max_x, max_x].
       All lengths have to be specified in µm.


    Arguments
    ---------
        lam: float
            wavelength
        m:  int
            number of the expansion orders retained
        n0: float
            refractive index of the surrounding
        max_x: float
            absolute extension of the spatial domain considered (equal extension into x and y direction);
        Nx:  int          
            Number of points used to discretize the space
            

    Returns
    -------
        u: 2d-array
            Matrix containing the amplitude of the plane wave as calculated from the expansion
            [matrix with dimension Nx x Nx]
        u_an: 2d-array
            Matrix containing the analytical solution of the plane wave
            [matrix with dimensions Nx x Nx]            
    '''
    #Defining initial parameters
    k = n0*2*np.pi/lam
    
    #Setting up arrays of coordinates
    x = np.linspace(-max_x,max_x,Nx)
    y = np.linspace(-max_x,max_x,Nx)
    X,Y = np.meshgrid(x,y)
    R = np.sqrt(X**2 + Y**2)
    
    #theta is calculated using np.arctan2 to ensure correct sign and angle within [-pi,pi]
    theta = np.arctan2(Y,X)
    
    #angular orders as vector
    m_vec = np.arange(-m,m+1,1)
    
    #defining analytical field
    u_an = np.exp(1j*k*X)
    u = np.zeros((u_an.shape),dtype=complex)
    
    #Summing terms over for loop
    for m in m_vec:
       u += (1j**m)*jv(m,k*R)*np.exp(-1j*m*theta)
    
    return u,u_an

#%%

#%%
def mie_pw(lam, m, n0, nc, Rad, pol, max_x, Nx):
    '''Calculates the amplitude distribution of a plane wave scattered at a cylinder.
       Computational domain has dimensions x = [-max_x, max_x], y = [-max_x, max_x].
       Center of the cylinder is located at the origin of the coordinate system.
       All lengths have to be specified in µm.


    Arguments
    ---------
        lam: float
            wavelength
        m:  int
            number of the expansion orders retained
        n0: float
            refractive index of the surrounding
        nc: float
            refractive index of the cylinder            
        Rad: float
            radius of the cylinder
        pol: string
            polarization indicator (either 'TE' or 'TM')
        max_x: float
            absolute extension of the spatial domain considerd (equal extension into x and y direction)
        Nx:  int          
            Number of points used to discretize the space
            

    Returns
    -------
        u: 2d-array
            Matrix containing the amplitude of the plane wave scattered at a cylinder
            [matrix with dimensions Nx x Nx]
    '''
    
    #Defining initial parameters
    n1 = n0
    n2 = nc
    k1 = n1*2*np.pi/lam
    k2 = n2*2*np.pi/lam
    z1 = k1*Rad
    z2 = k2*Rad   
    if pol == 'TE':
        p1 = p2 = 1
    if pol == 'TM':
        p1 = n1**2
        p2 = n2**2
        
    #Setting up coordinate arrays
    x = np.linspace(-max_x,max_x,Nx)
    y = np.linspace(-max_x,max_x,Nx)
    X,Y = np.meshgrid(x,y)
    R = np.sqrt(X**2 + Y**2)
    
    #theta is calculated using np.arctan2 to ensure correct sign and angle within [-pi,pi]
    theta = np.arctan2(Y,X)
    
    #angular orders as vector
    m_vec = np.arange(-m,m+1,1)
    
    #empty field array
    u = np.zeros((Nx,Nx),dtype=complex)
    
    #performing sum over orders and calcuation am and bm coefficients at each step
    for m in m_vec:
       am = (jvp(m,z2)*jv(m,z1)*p1*n2 - jvp(m,z1)*jv(m,z2)*p2*n1)/(hvp(m,z1)*jv(m,z2)*p2*n1 - hv(m,z1)*jvp(m,z2)*p1*n2)
       bm = (hvp(m,z1)*jv(m,z1)*p2*n1 - hv(m,z1)*jvp(m,z1)*p2*n1)/(hvp(m,z1)*jv(m,z2)*p2*n1 - hv(m,z1)*jvp(m,z2)*p1*n2)
       
       #only field in correct part of domain are summed
       u[R>Rad] += (1j**m)*am*hv(m,k1*R[R>Rad])*np.exp(-1j*m*theta[R>Rad])
       u[R<=Rad] += (1j**m)*bm*jv(m,k2*R[R<=Rad])*np.exp(-1j*m*theta[R<=Rad])
    
    #analytical field used as incident field for better accuracy
    u[R>Rad] += np.exp(1j*k1*X[R>Rad])
    return u

            


def mie_gauss(lam, waist, m, n0, nc, Rad, pol, max_x, Nx):
    '''Calculates the amplitude distribution of a Gaussian beam scattered at a cylinder.
       Beam waist is located in the plane x = 0.
       Computational domain has dimensions x = [-max_x, max_x], y = [-max_x, max_x].
       Center of the cylinder is located at the origin of the coordinate system.
       All lengths have to be specified in µm.


    Arguments
    ---------
        lam: float
            wavelength
        waist: float
            waist of the Gaussian beam         
        m:  int
            number of the expansion orders retained
        n0: float
            refractive index of the surrounding
        nc: float
            refractive index of the cylinder            
        Rad: float
            radius of the cylinder
        pol: string
            polarization indicator (either 'TE' or 'TM')
        max_x: float
            absolute extension of the spatial domain considerd (equal extension into x and y direction)
        Nx:  int          
            Number of points used to discretize the space
            

    Returns
    -------
        u: 2d-array
            Matrix containing the amplitude of the Gaussian beam scattered at a cylinder
            [matrix with dimensions Nx x Nx]
    '''
    #Defining initial parameters 
    n1 = n0
    n2 = nc
    k1 = n1*2*np.pi/lam
    k2 = n2*2*np.pi/lam
    z1 = k1*Rad
    z2 = k2*Rad   
    if pol == 'TE':
        p1 = p2 = 1
    if pol == 'TM':
        p1 = n1**2
        p2 = n2**2
        
    #Setting up coordinate arrays
    x = np.linspace(-max_x,max_x,Nx)
    y = np.linspace(-max_x,max_x,Nx)
    X,Y = np.meshgrid(x,y)
    R = np.sqrt(X**2 + Y**2)
    
    #theta is calculated using np.arctan2 to ensure correct sign and angle within [-pi,pi]
    theta = np.arctan2(Y,X)
    m_vec = np.arange(-m,m+1,1)
    
    #Empty initial field array
    u = np.zeros((Nx,Nx),dtype=complex)
    w = waist
    
    #Calculating frequency spacing and array of ky values to sum over
    rho = 2/w
    dky = rho/5
    Nk = np.rint(k1/dky)*2 + 1
    Nk = int(Nk)
    
    #ky in [-k1,k1] to avoid imaginary kx
    ky_vec = np.linspace(-k1,k1,Nk) 
    dky = np.abs(ky_vec[1] - ky_vec[0])
    
    #Calculating double sum to solve gaussian scattering problem 
    for ky in ky_vec:
        print(ky)
        kx = np.sqrt(k1**2 - ky**2)
        thetak = np.arctan2(ky,kx)
        ck = (w/np.sqrt(4*np.pi))*np.exp(-(w*ky/2)**2)*dky
        
        #The incident field is calculated from cartesian plane waves to avoid error due to partial wave expansion
        u[R>Rad] += ck*np.exp(1j*(ky*Y[R>Rad] + kx*X[R>Rad]))
        for m in m_vec:
            am = (jvp(m,z2)*jv(m,z1)*p1*n2 - jvp(m,z1)*jv(m,z2)*p2*n1)/(hvp(m,z1)*jv(m,z2)*p2*n1 - hv(m,z1)*jvp(m,z2)*p1*n2)
            bm = (hvp(m,z1)*jv(m,z1)*p2*n1 - hv(m,z1)*jvp(m,z1)*p2*n1)/(hvp(m,z1)*jv(m,z2)*p2*n1 - hv(m,z1)*jvp(m,z2)*p1*n2)
            
            u[R>Rad] += ck*(1j**m)*am*hv(m,k1*R[R>Rad])*np.exp(-1j*m*(theta[R>Rad]-thetak)) 
            u[R<=Rad] += ck*bm*(1j**m)*jv(m,k2*R[R<=Rad])*np.exp(-1j*m*(theta[R<=Rad]-thetak))

    return u