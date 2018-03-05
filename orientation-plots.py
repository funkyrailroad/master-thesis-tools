#!/usr/bin/env python

'''
This generates pictures of orientation plots and saves them as .png files
This only works for 2d
'''

from functions import *
from numpy import linalg as LA
import argparse
import h5py
import matplotlib.patches as patches
import matplotlib.pyplot as plt
import numpy as np
import sys


def main(filenames):
    for filename in filenames:
        box = get_box_dims(filename)
        p_step, p_time, position     =  get_position_data(filename)
        o_step, o_time, orientation  =  get_orientation_data(filename)
        v_step, v_time, velocity     =  get_velocity_data(filename)
        dims = get_dimension(filename)
        bounded_position = enforce_periodic_boundary_conditions(position, box)
        nsteps = orientation.shape[0]
        nparticles = orientation.shape[1]

    
        # how much bigger whole figure is that the simulation box
        figure_to_box_scale = 2.
        max_steps = 200

        if dims != 2:
            print "Houston, we have a problem. I can only make plots for 2d simulations."
            exit()

        i = 0
        for p, o in zip(bounded_position, orientation):
            fig = plt.figure()
            ax = fig.add_subplot(111)
            plt.quiver(p[:,0], p[:,1], o[:,0], o[:,1], units = 'height')
            ax.add_patch(patches.Rectangle(
                    (-box[0] / 2., -box[1] / 2.),
                    box[0],
                    box[1],
                    fill = False
                    )
                )
            plt.axis(figure_to_box_scale * np.array([- box[0] / 2., box[0] / 2., - box[1] / 2.,
                box[1] / 2.]))
            patches.Rectangle((- box[0] / 2., - box[1] / 2.), box[0], box[1] )
            plt.axes().set_aspect('equal')
            plt.savefig("antialignment/movie-test/%06d" %  i)
            plt.close()
            i += 1
            if i % 50 == 0:
                print "Completion: {:.2f}%".format(float(i) /  nsteps * 100)
            if i == max_steps:
                exit()



if __name__ == '__main__':
    main(sys.argv[1:])
