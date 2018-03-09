#!/usr/bin/env python

'''
NEED TO ENFORCE PERIODIC BOUNDARY CONDITIONS ON THE POSITION
This only works for 2d
'''

from functions import *
import argparse
import h5py
import math
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
from skimage import io
from scipy import fftpack
import sys
from numpy import linalg as LA

def full_ordered_frequencies_1d(array):
    return np.concatenate( ( -array[:0:-1], array[:-1]) )

def full_ordered_modes_2d(array):
    # this takes an array of positive frequencies and mirrors it about both
    # (2d) axes to have something that fftshift can work with
    array = np.concatenate( ( np.conj( array[:,:0:-1] ), array[:,:-1]), axis = 1)
    array = np.concatenate( ( np.conj( array[:0:-1] ), array[:-1]), axis = 0)
    return array

def gaussian_spacing_dx( sigma ):
    return 2 * sigma * np.sqrt( 2 * np.log( 2 ) )

def gaussian_spacing_sigma( dx ):
    return dx / ( 2 * np.sqrt( 2 * np.log( 2 ) ) )


def gaussian_convolution(positions_snapshot, r_0, sigma):
    # this returns the convolution of a gaussian with a distribution of
    # particles represented by delta peaks

    return np.exp( 
           np.sum( 
           - ( positions_snapshot  - r_0 ) ** 2
           , axis = 1 )  # this sums the components of the position vector for a single particle
           / ( 2 * sigma ** 2 ) ).sum() # this final sum sums over all particles
    


def smart_analytical_option(positions_snapshot, box_length, cell_length):

    nbins = box_length / cell_length
    sigma = gaussian_spacing_sigma( cell_length )

    if nbins != int( nbins ):
        print 'Error: box length not an exact multiple of cell length.\n\
                Box length: \t{0}\n\
                Cell length: \t{1}'.format(box_length, cell_length)
        exit()


    # these are shifted so that the coordinates of the bins are in the middle of the bins
    full_x = np.linspace( - box_length / 2, box_length / 2, nbins + 1 )[:-1] + 0.5 * cell_length
    full_y = np.linspace( - box_length / 2, box_length / 2, nbins + 1 )[:-1] + 0.5 * cell_length
    X, Y = np.meshgrid( full_x, full_y )
    Z = np.zeros_like( X )

    for iy, y in enumerate( full_y ):
        for ix, x in enumerate( full_x ):
            r_0 = np.array( [x, y] )
            Z[ iy, ix ] = gaussian_convolution( positions_snapshot, r_0, sigma)



    #plotting stuff
    plt.scatter(positions_snapshot[:,0], positions_snapshot[:,1], zorder = 10)
    ax = plt.contourf(X, Y, Z)
    plt.title("Analytical Postprocessing Initial Distribution")
    plt.axes().set_aspect('equal')
    cax = plt.axes([0.90, 0.1, 0.025, 0.8])
    plt.colorbar(cax=cax)
    #plot_periodic_density_contour_2d(box, Z, cell_length)
    return X, Y, Z




def main(filename):
    box = get_box_dims(filename)
    p_step, p_time, position     =  get_position_data(filename)
    o_step, o_time, orientation  =  get_orientation_data(filename)
    v_step, v_time, velocity     =  get_velocity_data(filename)
    #bounded_position = enforce_periodic_boundary_conditions(position, box)
    dims = get_dimension(filename)
    dm_step, dm_time, modes, wavevector = get_density_mode_data(filename)

    sigma = 0.1
    sigma = np.pi / ( max(wavevector[:,0]) * np.sqrt( 2 * np.log(2) ) )
    cell_length = 0.2 

    #ANALYTICAL OPTION
    smart_analytical_option(position[0], box[0], cell_length)
    plt.show()
    exit()
    #END ANALYTICAL OPTION



if __name__ == '__main__':
    for filename in sys.argv[1:]:
        main(filename)

