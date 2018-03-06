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
        


def quadrant_to_full(kx, ky, Z, axarray, positions):
    '''
    x and y are one dimensional positive frequency arrays
    Z is the complex valued matrix that corresponds to the frequencies of the x
    and y vectors
    '''


    full_kx = full_ordered_frequencies_1d( kx )
    full_ky = full_ordered_frequencies_1d( ky )
    full_Z = full_ordered_modes_2d( Z )

    print 2 * np.pi / kx[1]
    box_length = 2 * np.pi / kx[1] 
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
        k_vec = np.array( [ kx[ k1_ind ], ky[ k2_ind ] ] )
        r_0 = np.array( [1.0, 0.0], dtype = float )
        if modes_matrix[k2_ind, k1_ind] == 0:
            gaussian = np.exp( -1j * np.pi * ( np.dot(k_vec, r_0) ) - 0.5 * np.pi
                ** 2 * sigma**2 * ( kx[k1_ind]**2 + ky[k2_ind]**2 ) )
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

def main(filename):
    box = get_box_dims(filename)
    p_step, p_time, position     =  get_position_data(filename)
    o_step, o_time, orientation  =  get_orientation_data(filename)
    v_step, v_time, velocity     =  get_velocity_data(filename)
    #bounded_position = enforce_periodic_boundary_conditions(position, box)
    dims = get_dimension(filename)
    dm_step, dm_time, modes, wavevector = get_density_mode_data(filename)

    modes = np.copy(modes)
    modes_shifted = np.fft.fftshift(modes)

    sigma = 0.1
    cell_length = 0.2 

    # getting the data of one timestep to play around with
    snapshot = np.copy(modes[0,:])
    k1_u, k2_u, modes_matrix = populate_modes_matrix(wavevector, snapshot, sigma)

    # PURELY FOR OUTPUT, NOT RELEVANT FOR THE CALCULATION
    # sorting the data first by k1, then by k2
    ind = np.lexsort((wavevector[:,1], wavevector[:,0]))
    array = wavevector[ind]
    snapshot_ordered = snapshot[ind]
    # outputting the ordered k vectors and their corresponding (complex)
    # amplitudes
    for i, val in enumerate(array):
        print val, np.abs(snapshot_ordered[i])

    # the associated positions for the wavevectors
    # how to calculate these?
    






    # testing with manually calculating the fourier modes on positions
    test_positions = np.array([ 
                            [0, -0]
                           ,[0, 1.5]
                           #,[0, -1.5]
                           ,[-1., -1.]
                           #,[1., 1.7]
                           #,[-1, -1] 
                           ], dtype=float)
    box = np.array( [4, 4] )

    test_modes = np.zeros(wavevector.shape[0], dtype=complex)
    for ik, k in enumerate( wavevector ):
        for ir, r in enumerate( test_positions ):
            test_modes[ik] += np.exp(1j * np.dot(k, r))


    f, axarr = plt.subplots(2, 2)
    kx, ky, test_modes_matrix = populate_modes_matrix(wavevector, test_modes, sigma)



    for i in range(kx.shape[0]):
        wavelength =  2 * np.pi / np.array([kx[i], ky[i]]) 
    #exit()
    quadrant_to_full(kx, ky, test_modes_matrix, axarr, test_positions) 
    plt.show()
    exit()

    










    f, axarr = plt.subplots(2, 2)
    quadrant_to_full(k1_u, k2_u, modes_matrix, axarr) 
    exit()

    # testing with an image of a grid
    image = io.imread('/home/mi/jatwell/playground/images/horizontal-vertical-lines.png')[11:-16, 13:-17]
    M, N = image.shape
    F2_shifted = fftpack.fft2(image)
    F2_ordered = fftpack.ifftshift(F2_shifted)
    print F2_shifted.shape
    print M, N
    #axarr[0,0].imshow(image, cmap ='gray')
    #axarr[0,1].imshow( np.log( 1 + np.abs( F2_ordered ) ), cmap ='viridis' )
    #axarr[0,1].contourf( np.log( 1 + np.abs( F2_ordered ) ), cmap ='viridis' )
    #axarr[1,0].contourf( np.fft.ifft2( F2_shifted ) , cmap ='viridis' )
    plt.show()


    # only the positive frequency modes are picked out
    plus =  F2_shifted[:M // 2 + 1, : N // 2 + 1]

    print F2_ordered[0] == full_ordered_modes_2d(plus)[0]
    print max(np.stack( ( full_ordered_modes_2d(plus)[0], F2_ordered[0], np.abs(full_ordered_modes_2d(plus)[0] - F2_ordered[0]) ) ).T[:,2])
        


    kx_full = np.fft.fftfreq(N)
    ky_full = np.fft.fftfreq(M)
    kx = kx_full[: N // 2 + 1]
    kx[-1] *= -1
    ky = ky_full[: M // 2 + 1]
    ky[-1] *= -1

    f, axarr = plt.subplots(2, 2)
    quadrant_to_full(kx, ky, plus, axarr) 
    #plt.show()




    exit()
    




    # modes is rho_k!
    for snapshot in modes:
        print snapshot.shape
        print wavevector.shape
        print snapshot
        back_transformed = np.fft.ifft2(snapshot)
        exit()
        exit()
        plt.contourf(wavevector[:,0], wavevector[:,1], snapshot)
        plt.show()
    

        exit()


    nsteps = position.shape[0]
    nparticles = position.shape[1]

    ncells = int(math.floor( box[0] / cell_length ))







if __name__ == '__main__':
    for filename in sys.argv[1:]:
        main(filename)

