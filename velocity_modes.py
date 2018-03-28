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
    vm_step, vm_time, v_modes, v_wavevector = get_velocity_mode_data(filename)

    nsteps = orientation.shape[0]
    nparticles = orientation.shape[1]
    sigma = 0.1
    max_steps = 1

    for i in v_modes[0]:
        print i
    exit()

    '''
    bounded_position = np.array([[   [ 0,    0  ]
                                    ,[ 0.5,  0.0]
                                    ,[ 1.5,  0.0]
                                                 ]])
    orientation = np.array([[  [ 0,    0  ]
                              ,[-0.5,  1.0]
                              ,[-1.5,  0.2]
                                            ]])
    '''

    #v_modes = get_velocity_modes( bounded_position, orientation, v_0, wavevector, max_steps )

    

    i = 0
    for position_snapshot, orientation_snapshot, modes_snapshot in zip(bounded_position, orientation, v_modes):

        k1_u, k2_u, d_modes_matrix = populate_modes_matrix(wavevector,
                d_modes[i,:], sigma)
        k1_u, k2_u, v_modes_matrix_x = populate_modes_matrix(wavevector,
                modes_snapshot[:,0], sigma)
        k1_u, k2_u, v_modes_matrix_y = populate_modes_matrix(wavevector,
                modes_snapshot[:,1], sigma)
        #k1_u, k2_u, v_modes_matrix_z = populate_modes_matrix(wavevector,
                #modes_snapshot[:,2], sigma)

        ft_d_modes   = np.real( np.fft.fftshift(np.fft.ifft2( np.fft.ifftshift(
            d_modes_matrix ) )) )

        # NOTE the arrays here are spatial but pertain to the orientation along
        # a given axis. The arrays are indexed like [y,x] which is opposite the
        # traditional convention
        ft_v_modes_x = np.real( np.fft.fftshift(np.fft.ifft2( np.fft.ifftshift(
            v_modes_matrix_x ) )) )
        ft_v_modes_y = np.real( np.fft.fftshift(np.fft.ifft2( np.fft.ifftshift(
            v_modes_matrix_y ) )) )
        #ft_v_modes_z = np.real( np.fft.fftshift(np.fft.ifft2( np.fft.ifftshift(
            #v_modes_matrix_z ) )) )

        #ft_d_modes = normalize_density_mode_matrix( ft_d_modes, nparticles)
        #ft_v_modes_x = normalize_other_mode_matrix( ft_v_modes_x, ft_d_modes)
        #ft_v_modes_y = normalize_other_mode_matrix( ft_v_modes_y, ft_d_modes)
        #ft_v_modes_z = normalize_other_mode_matrix( ft_v_modes_z, ft_d_modes)

        maxind = np.unravel_index(np.argmax(ft_v_modes_x) , ft_v_modes_x.shape)
        a = slice(38, 41)
        b = slice(58, 61)
        print np.sum(ft_v_modes_x[a, b])


        r = np.array([1, 2.])
        ind = coord_to_ind(ft_v_modes_x, box, r)
        b = 2
        #ft_v_modes_x[ind[0]-b:ind[0]+b, ind[1] - b: ind[1] + b] = 1

        #ft_v_modes_x[a, b] = 1

        f, axarr = plt.subplots(2, 2)
        #full_velocity_orientation_plots( k1_u, k2_u, d_modes_matrix,
                #ft_d_modes, axarr, position_snapshot, box )
        full_velocity_orientation_plots( k1_u, k2_u, v_modes_matrix_x,
                ft_v_modes_x, axarr, position_snapshot, box )
        #full_velocity_orientation_plots( k1_u, k2_u, v_modes_matrix_y,
                #ft_v_modes_y, axarr, position_snapshot, box )
        #full_velocity_orientation_plots( k1_u, k2_u, v_modes_matrix_z,
                #ft_v_modes_z, axarr, position_snapshot, box )
        plt.savefig("antialignment/movie-test/velocity-modes-%06d.png" %  i)
        plt.show()
        plt.close()


        '''
        plt.quiver(position_snapshot[:,0], position_snapshot[:,1],
                orientation_snapshot[:,0], orientation_snapshot[:,1], 
                units = 'height')
        '''

        plt.quiver( ft_v_modes_x, ft_v_modes_y, units = 'height')

        plt.show()
        exit()

        i += 1
        if i % 5 == 0:
            print "Completion: {:.2f}%".format(float(i) /  nsteps * 100)
        if i == max_steps:
            exit()


if __name__ == '__main__':
    for filename in sys.argv[1:]:
        main(filename)

