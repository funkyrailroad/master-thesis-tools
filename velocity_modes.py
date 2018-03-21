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

        


def main(filename):
    box = get_box_dims(filename)
    p_step, p_time, position     =  get_position_data(filename)
    o_step, o_time, orientation  =  get_orientation_data(filename)
    v_step, v_time, velocity     =  get_velocity_data(filename)
    bounded_position = enforce_periodic_boundary_conditions(position, box)
    dims = get_dimension(filename)
    v_0 = get_v_0(filename)[-1]
    dm_step, dm_time, d_modes, wavevector = get_density_mode_data(filename)
    #vm_step, vm_time, v_modes, v_wavevector = get_velocity_mode_data(filename)
    print np.sum(orientation[1]**2, axis = 1)
    exit()


    nsteps = orientation.shape[0]
    sigma = 0.1
    max_steps = 5

    orientation = -orientation
    v_modes = get_velocity_modes( position, orientation, v_0, wavevector, max_steps )
    print d_modes.shape
    print v_modes.shape

    i = 0
    for position_snapshot, modes_snapshot in zip(bounded_position, v_modes):
        print orientation[i]
        k1_u, k2_u, d_modes_matrix = populate_modes_matrix(wavevector,
                d_modes[i,:], sigma)
        k1_u, k2_u, v_modes_matrix_x = populate_modes_matrix(wavevector,
                modes_snapshot[:,0], sigma)
        k1_u, k2_u, v_modes_matrix_y = populate_modes_matrix(wavevector,
                modes_snapshot[:,1], sigma)
        #k1_u, k2_u, v_modes_matrix_z = populate_modes_matrix(wavevector,
                #modes_snapshot[:,2], sigma)

        # what is going on with this hack? 
        # double check it when switching to halmd velocity mode module
        v_modes_matrix_x = np.flip(v_modes_matrix_x, axis = 0)
        v_modes_matrix_x = np.flip(v_modes_matrix_x, axis = 1)
        v_modes_matrix_y = np.flip(v_modes_matrix_y, axis = 0)
        v_modes_matrix_y = np.flip(v_modes_matrix_y, axis = 1)
        #v_modes_matrix_z = np.flip(v_modes_matrix_z, axis = 0)
        #v_modes_matrix_z = np.flip(v_modes_matrix_z, axis = 1)


        ft_d_modes = np.real( np.fft.fftshift(np.fft.ifft2( np.fft.ifftshift(
            d_modes_matrix ) )) )
        ft_v_modes_x = np.real( np.fft.fftshift(np.fft.ifft2( np.fft.ifftshift(
            v_modes_matrix_x ) )) )
        ft_v_modes_y = np.real( np.fft.fftshift(np.fft.ifft2( np.fft.ifftshift(
            v_modes_matrix_y ) )) )
        #ft_v_modes_z = np.real( np.fft.fftshift(np.fft.ifft2( np.fft.ifftshift(
            #v_modes_matrix_z ) )) )

        ft_d_modes = normalize_density_mode_matrix( ft_d_modes )
        ft_v_modes_x = normalize_other_mode_matrix( ft_v_modes_x, ft_d_modes)
        ft_v_modes_y = normalize_other_mode_matrix( ft_v_modes_y, ft_d_modes)
        #ft_v_modes_z = normalize_other_mode_matrix( ft_v_modes_z, ft_d_modes)

        f, axarr = plt.subplots(2, 2)
        #new_wavevector_module( k1_u, k2_u, d_modes_matrix, ft_d_modes, axarr,
                #position_snapshot, box[0] )
        #new_wavevector_module( k1_u, k2_u, v_modes_matrix_x, ft_v_modes_x, axarr,
                #position_snapshot, box[0] )
        new_wavevector_module( k1_u, k2_u, v_modes_matrix_y, ft_v_modes_y, axarr,
                position_snapshot, box[0] )
        #new_wavevector_module( k1_u, k2_u, v_modes_matrix_z, ft_v_modes_z, axarr,
                #position_snapshot, box[0] )
        plt.savefig("antialignment/movie-test/velocity-modes-%06d.png" %  i)
        plt.close()
        i += 1
        if i % 5 == 0:
            print "Completion: {:.2f}%".format(float(i) /  nsteps * 100)
        if i == max_steps:
            exit()


if __name__ == '__main__':
    for filename in sys.argv[1:]:
        main(filename)

