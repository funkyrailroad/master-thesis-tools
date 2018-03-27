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
    dm_step, dm_time, modes, wavevector = get_density_mode_data(filename)
    nsteps = orientation.shape[0]
    nparticles = orientation.shape[1]
    sigma = 0.1
    max_steps = 1

    i = 0
    for position_snapshot, modes_snapshot in zip(bounded_position, modes):

        position_snapshot = np.array([ [0, 0]
                                      ,[1, 0]
                                      ,[2.5, 0.0]
                                        ])
        snap2 = calculate_density_mode_snapshot(wavevector, position_snapshot)
        modes_snapshot = snap2

        k1_u, k2_u, modes_matrix = populate_modes_matrix(wavevector,
                modes_snapshot, sigma)

        # taking the real part because it started as a real image, although the
        # imaginary part is very small to begin with
        ft_modes = np.real( np.fft.fftshift(np.fft.ifft2( np.fft.ifftshift(
            modes_matrix ) )) )

        # plotting square on contour plot ( unneeded )
        maxind = np.unravel_index(np.argmax(ft_modes) , ft_modes.shape)

        #ft_modes = normalize_density_mode_matrix(ft_modes, nparticles)

        f, axarr = plt.subplots(2, 2)
        full_density_position_plots( k1_u, k2_u, modes_matrix, ft_modes, axarr,
                position_snapshot, box[0] ) 
        plt.savefig("antialignment/movie-test/density-modes-%06d.png" %  i)
        plt.show()
        #plt.close()
        i += 1
        if i % 5 == 0:
            print "Completion: {:.2f}%".format(float(i) /  nsteps * 100)
        if i == max_steps:
            exit()







if __name__ == '__main__':
    for filename in sys.argv[1:]:
        main(filename)

