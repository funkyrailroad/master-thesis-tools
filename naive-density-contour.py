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
import sys
from numpy import linalg as LA

        
def main(filenames):
    for filename in filenames:

        box = get_box_dims(filename)
        p_step, p_time, position     =  get_position_data(filename)
        o_step, o_time, orientation  =  get_orientation_data(filename)
        v_step, v_time, velocity     =  get_velocity_data(filename)
        #bounded_position = enforce_periodic_boundary_conditions(position, box)
        dims = get_dimension(filename)

        nsteps = position.shape[0]
        nparticles = position.shape[1]

        cell_length = 0.2 
        ncells = int(math.floor( box[0] / cell_length ))


        for (t, timestep) in enumerate(position):

            # making the cell lists
            cells, particles = create_cell_linked_list(timestep, box,
                    cell_length, dims)

            # calculating the average velocity and density of each cell
            avg_rho, avg_v = particle_and_velocity_binning(cells,
                    particles, dims, orientation[t], cell_length, ncells)

            if dims == 2:

                fig = plt.figure()
                ax = fig.add_subplot(111)
                plot_periodic_density_contour_2d(box, avg_rho, cell_length)
                plt.axes().set_aspect('equal')
                plot_cells(box, cell_length)
                ax.add_patch(patches.Rectangle( -box / 2, *box, fill = False,
                    color = 'r'))
                for part in timestep:
                    plt.scatter(*part, color = 'k')
                plt.show() 

                exit()





if __name__ == '__main__':
    main(sys.argv[1:])

