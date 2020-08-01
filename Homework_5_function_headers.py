'''Homework 5, Computational Photonics, SS 2020:  Scattering method.
'''

import numpy as np
import scipy.special as sp



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
    pass
            
def dr_jv(m,x):
    return x*sp.jv(m-1,x)-m*sp.jv(m,x)
def dr_hv(m,x):
    return m*sp.hankel1(m,x)/x-sp.hankel1(m+1,x)

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
    #create coordinates
    x = np.linspace(-max_x,max_x,Nx)
    xv,yv = np.meshgrid(x,x)
    r=(xv**2+yv**2)**0.5
    theta=np.arctan(yv/xv)
    theta[xv<=0]=theta[xv<=0]+np.pi
    
    #new matrices
    u = np.zeros([Nx,Nx],complex)
    
    #identify polarization
    if pol=='TE':
        p0=1
        pc=1
    elif pol=='TM':
        p0=1/n0**2
        pc=1/nc**2
    else:
        raise ValueError('Please choose TE or TM mode.')
    
    #create new marices for mie coefficients
    am=np.zeros(2*m+1,complex)
    bm=np.zeros(2*m+1,complex)
    
    #calculate mie coefficients
    k=2*np.pi/lam
    for i in range(-m,m+1): 
        am[i+m]=(sp.jv(i,k*n0*Rad)*dr_jv(i,k*nc*Rad)*p0-
            sp.jv(i,k*nc*Rad)*dr_jv(i,k*n0*Rad)*pc)/(sp.jv(i,k*nc*Rad)*dr_hv(i,k*n0*Rad)*pc-
            sp.hankel1(i,k*n0*Rad)*dr_jv(i,k*nc*Rad)*p0)
        bm[i+m]=(sp.jv(i,k*n0*Rad)*dr_hv(i,k*n0*Rad)*pc-
            sp.hankel1(i,k*n0*Rad)*dr_jv(i,k*nc*Rad)*pc)/(sp.jv(i,k*nc*Rad)*dr_hv(i,k*n0*Rad)*pc-
            sp.hankel1(i,k*n0*Rad)*dr_jv(i,k*nc*Rad)*p0)
   
        u[r<=Rad]=u[r<=Rad]+bm[i+m]*1j**i*sp.jv(i,k*nc*r[r<=Rad])*np.exp(-1j*i*theta[r<=Rad])  #internal field
        u[r>Rad]=u[r>Rad]+am[i+m]*1j**i*sp.hankel1(i,k*n0*r[r>Rad])*np.exp(-1j*i*theta[r>Rad]) #scattered field
    
    u=u+np.exp(1j*k*n0*xv)  #total field
                                                          
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
    pass
