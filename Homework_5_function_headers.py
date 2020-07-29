'''Homework 5, Computational Photonics, SS 2020:  Scattering method.
'''

import numpy as np
import scipy.special



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
    pass

            


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