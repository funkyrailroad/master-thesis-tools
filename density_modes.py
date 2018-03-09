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
        
def new_wavevector_module(full_kx, full_ky, full_Z, axarray, positions, box_length):


    full_x = np.linspace(-box_length / 2, box_length / 2, full_Z.shape[0])
    full_y = np.linspace(-box_length / 2, box_length / 2, full_Z.shape[0])

    print
    print 'here'
    print
    k1max = max(abs(full_kx))
    print k1max
    print 2 * np.pi / box_length * 39
    print k1max == 2 * np.pi / 0.1
    print full_x.shape
    print "supposed real space spacing from k1max"
    print 2 * np.pi / k1max
    print "real space spacing"
    print full_x[1] - full_x[0]
    

    # quadrant plot
    ax = axarray[0]

    ax[0].scatter( positions[:,0], positions[:,1] )
    ax[0].set_title("Actual Particle Distribution")
    ax[0].set_aspect('equal')
    ax[0].set_xlim([ - box_length / 2, box_length / 2 ] )
    ax[0].set_ylim([ - box_length / 2, box_length / 2 ] )

    cont = ax[1].contourf( full_x, full_y, np.fft.fftshift(np.fft.ifft2( np.fft.ifftshift( full_Z ) ))) 
    #cont = ax[1].contourf( np.fft.ifft2( np.fft.ifftshift( full_Z ) ))
    ax[1].set_aspect('equal')
    ax[1].set_title("IFT (obtained particle distribution)")
    cax = plt.axes([0.90, 0.1, 0.025, 0.35])
    plt.colorbar(cont, cax=cax)


    ax = axarray[1]

    ax[0].set_aspect('equal')
    ax[0].contourf(full_kx, full_ky, np.abs( full_Z ))
    ax[0].set_title('Positive and negative frequencies')

    ax[1].scatter( positions[:,0], positions[:,1], zorder = 10)
    cont = ax[1].contourf( full_x, full_y, np.fft.fftshift(np.fft.ifft2( np.fft.ifftshift( full_Z ) ))) 
    #cont = ax[1].contourf( np.fft.ifft2( np.fft.ifftshift( full_Z ) ))
    ax[1].set_aspect('equal')
    ax[1].set_title("IFT and Actual Positions")
    cax = plt.axes([0.90, 0.1, 0.025, 0.35])
    plt.colorbar(cont, cax=cax)


def quadrant_to_full(kx, ky, Z, axarray, positions, box_length):
    '''
    x and y are one dimensional positive frequency arrays
    Z is the complex valued matrix that corresponds to the frequencies of the x
    and y vectors
    '''


    full_kx = full_ordered_frequencies_1d( kx )
    full_ky = full_ordered_frequencies_1d( ky )
    full_Z = full_ordered_modes_2d( Z )

    print 2 * np.pi / kx[1]
    #box_length = 2 * np.pi / kx[1] 
    full_x = np.linspace(-box_length / 2, box_length / 2, full_Z.shape[0])
    full_y = np.linspace(-box_length / 2, box_length / 2, full_Z.shape[0])

    # quadrant plot
    ax = axarray[0]
    ax[0].set_aspect('equal')
    ax[0].contourf(kx, ky, np.abs( Z ))
    ax[0].set_title('Positive Frequencies')

    # this ifft doesn't make any sense, Z is not in the right format
    #cont = ax[1].contourf( np.fft.ifft2( Z ) )
    #cax = plt.axes([0.90, 0.55, 0.025, 0.35])
    #plt.colorbar(cont, cax=cax)
    print box_length
    ax[1].scatter( positions[:,0], positions[:,1] )
    ax[1].set_title("Example Particle Distribution")
    ax[1].set_aspect('equal')
    ax[1].set_xlim([ - box_length / 2, box_length / 2])
    ax[1].set_ylim([ - box_length / 2, box_length / 2])

    # full plot
    ax = axarray[1]
    ax[0].set_aspect('equal')
    ax[0].contourf(full_kx, full_ky, np.abs( full_Z ))
    ax[0].set_title('Positive and negative frequencies')

    cont = ax[1].contourf( full_x, full_y, np.fft.fftshift(np.fft.ifft2( np.fft.ifftshift( full_Z ) ))) 
    #cont = ax[1].contourf( np.fft.ifft2( np.fft.ifftshift( full_Z ) ))
    ax[1].set_aspect('equal')
    ax[1].set_title("IFT (obtained particle distribution)")
    cax = plt.axes([0.90, 0.1, 0.025, 0.35])
    plt.colorbar(cont, cax=cax)

def populate_modes_matrix(wavevector, mode_snapshot, sigma):

    # make ordered lists of the values of k that were used
    # make it some sort of array so it's more easily generalized into 3d
    # k$ind are the indicies of wavevector that contain the unique
    # values of wavevector
    kx = np.unique( wavevector[:,0] )
    ky = np.unique( wavevector[:,1] )

    # make an n-dimensional matrix for holding the mode data
    modes_matrix = np.zeros((len(ky), len(kx)), dtype=complex)

    for i in range(mode_snapshot.shape[0]):
        # the indicies of the area in which the components of wavevector should
        # go are found
        k1_ind = np.where(kx <= wavevector[i,0])[0][-1]
        k2_ind = np.where(ky <= wavevector[i,1])[0][-1]
        print wavevector[0]
        print np.where(kx <= wavevector[i,0])[0][-1]
        print np.where(ky <= wavevector[i,1])[0][-1]
        print kx[k1_ind]
        print ky[k2_ind]
        print kx[k1_ind] == wavevector[0,0]
        print ky[k2_ind] == wavevector[0,1]
        #exit()
        k_vec = np.array( [ kx[ k1_ind ], ky[ k2_ind ] ] )
        #r_0 = np.array( [1.0, 0.0], dtype = float )
        if modes_matrix[k2_ind, k1_ind] == 0:

            # fully fourier transformed gaussian with r_0
            #gaussian = np.exp( -1j * ( np.dot(k_vec, r_0) ) - 0.5 * sigma ** 2 * ( kx[k1_ind]**2 + ky[k2_ind]**2 ) )

            # only real part of exponent ( without r_0 )
            gaussian = np.exp( - 0.5 * sigma ** 2 * ( kx[k1_ind]**2 + ky[k2_ind]**2 ) )

            # with prefactors and only real part of exponent
            #gaussian = sigma * np.sqrt( 2 * np.pi) * np.exp( - 0.5 * sigma ** 2 * ( kx[k1_ind]**2 + ky[k2_ind]**2 ) )

            #modes_matrix[k2_ind, k1_ind] = gaussian * mode_snapshot[i]
            modes_matrix[k2_ind, k1_ind] = mode_snapshot[i]
        else:
            print 'Warning: Uneccesary mode recalculation.'
            if modes_matrix[k2_ind, k1_ind] != mode_snapshot[i]:
                print "Error: Mode caclulated doubly and differently for a\
                        wavevector."
                #exit()
        #exit()
    return kx, ky, modes_matrix

def gaussian_spacing_dx( sigma ):
    return 2 * sigma * np.sqrt( 2 * np.log( 2 ) )

def gaussian_spacing_sigma( dx ):
    return dx / ( 2 * np.sqrt( 2 * np.log( 2 ) ) )

def print_modes_sorted_by_wavevector(wavevector, modes):
    # PURELY FOR OUTPUT, NOT RELEVANT FOR THE CALCULATION
    # sorting the data first by k1, then by k2
    ind = np.lexsort((wavevector[:,1], wavevector[:,0]))
    array = wavevector[ind]
    modes_ordered = modes[ind]
    # outputting the ordered k vectors and their corresponding absolute values
    # amplitudes
    for i, val in enumerate(array):
        print val, np.abs(modes_ordered[i])

def main(filename):
    box = get_box_dims(filename)
    p_step, p_time, position     =  get_position_data(filename)
    o_step, o_time, orientation  =  get_orientation_data(filename)
    v_step, v_time, velocity     =  get_velocity_data(filename)
    #bounded_position = enforce_periodic_boundary_conditions(position, box)
    dims = get_dimension(filename)
    dm_step, dm_time, modes, wavevector = get_density_mode_data(filename)

    modes = np.copy(modes)
    # getting the data of one timestep to play around with
    modes_snapshot = np.copy(modes[0,:])

    sigma = 0.05

    #print_modes_sorted_by_wavevector(wavevector, modes_snapshot)

    k1_u, k2_u, modes_matrix = populate_modes_matrix(wavevector, modes_snapshot, sigma)
    #exit()

    sigma = np.pi / ( max(wavevector[:,0]) * np.sqrt( 2 * np.log(2) ) )
    print sigma

    k1_u, k2_u, modes_matrix = populate_modes_matrix(wavevector, snapshot, sigma)
    print gaussian_spacing_dx(sigma)


    f, axarr = plt.subplots(2, 2)
    new_wavevector_module(k1_u, k2_u, modes_matrix, axarr, position[0], box[0]) 
    #plt.show()
    plt.savefig('123.eps')




if __name__ == '__main__':
    for filename in sys.argv[1:]:
        main(filename)

