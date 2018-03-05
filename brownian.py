#!/usr/bin/python
'''
The current condition of this code requires specific input. An equal number of
gpu and host files are to be supplied as arguments with the gpu files first
and the host files immediately after.
'''


import sys
import numpy as np
import argparse
import h5py
import matplotlib.pyplot as plt
import copy
from scipy.optimize import curve_fit
from math import log10, floor

def round_to_n(x, n): 
    return round(x, -int(floor(log10(x))) + (n - 1))

def lin_fit(x, a, b):
    return a*x + b

def lin_fit_through_null(x, a):
    return a*x

def logy_func(x, a, b):
    return a*np.exp(- b*x)

def calc_data(data):
    return new_data


def logy_fit(x, y):
    '''
    pcov is useless right now because of the transform of the first optimal
    parameter
    '''

    logy = np.log(y)
    logy = np.nan_to_num(logy)
    popt, pcov = curve_fit(lin_fit, x, logy)
    dummy = popt[1]
    popt[1] = popt[0]
    popt[0] = np.exp(dummy)
    return popt, pcov


def plot_data(x, y, title=None, xlabel=None, ylabel=None, legend=None,
        xscale='linear', yscale='linear', legend_loc=None):

    lineObjects = plt.plot(x, y)
    plt.title(title)
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.xscale(xscale)
    plt.yscale(yscale)
    plt.legend(lineObjects, legend, loc = legend_loc)

def main():

    ##############################
    # Clean data
    ##############################

    parser = argparse.ArgumentParser(prog='2d-simulation.py')
    # define and parse command line arguments
    parser.add_argument('input', metavar='INPUT', help='H5MD input file with data for state variables')
    args = parser.parse_args()

    # open and read data file
    H5 = h5py.File(args.input, 'r')
    
    p_full = H5['particles/all/position']
    p_step = np.array(p_full['step'])
    p_time = np.array(p_full['time'])
    position  = np.array(p_full['value'])

    o_full = H5['particles/all/orientation']
    o_step = np.array(o_full['step'])
    o_time = np.array(o_full['time'])
    orientation  = np.array(o_full['value'])

    v_full = H5['particles/all/velocity']
    v_step = np.array(v_full['step'])
    v_time = np.array(v_full['time'])
    velocity  = np.array(v_full['value'])


    # dynamics/* are in different format that particles/*
    msd_full = H5['dynamics/all/mean_square_displacement']
    msd_count = np.array(msd_full['count']).flatten() # what is count?
    msd_time = np.array(msd_full['time']).flatten()
    msd  = np.array(msd_full['value']).flatten()
    d_msd  = np.array(msd_full['error']).flatten()


    ocf_full = H5['dynamics/all/orientational_autocorrelation']
    ocf_count = np.array(ocf_full['count']).flatten() # what is count?
    ocf_time = np.array(ocf_full['time']).flatten()
    ocf  = np.array(ocf_full['value']).flatten()
    d_ocf  = np.array(ocf_full['error']).flatten()


    # organize data into a single array so it can be sorted
    sort_column_index = 0 # index of colum by which data should be sorted (0-sorting by time, 1-sorting by value)

    data = np.array((ocf_time, ocf, d_ocf)).T 
    ocf_time, ocf, d_ocf = data[np.argsort(data[:,sort_column_index])].T
    
    data = np.array((msd_time, msd, d_msd)).T 
    msd_time, msd, d_msd = data[np.argsort(data[:,sort_column_index])].T

    step = copy.deepcopy(p_step)
    time = copy.deepcopy(p_time)

    




    nsteps = orientation.shape[0]
    nparticles = orientation.shape[1]
    dims = orientation.shape[2]


    num_of_vars = 4
    num_of_files = len(sys.argv) - 1


    ###################
    ## Plotting
    ###################

    plt.figure()
    D_t = 1         # translational diffusion constant
    D_r = 1         # rotational diffusion constant
    D_par = D_t
    D_perp = D_t
    D_eff = ( D_par + 2 * D_perp ) / 3 # effective diffusion constant 
    v_0 = 30        # active motion
    plt.plot(time, np.exp(- 1. * D_t * time), label = "Theory exp(-Dt)") # prefactor is -1 for 2d, -2 for 3d

    plot_data(ocf_time, ocf, "GPU Side OCF",
        "Time", "OCF", ['Simulation'], 'log', 'linear', legend_loc = 'lower left')
    plt.savefig("ocf_gpu.png")
    plt.close('all')


    '''
    plt.figure()
    plot_data(msd_time, msd, "GPU Side",
        "Time", "msd", ['Simulation'], 'linear', 'linear', legend_loc = 'upper left')
    plt.savefig("msd_gpu.png")
    plt.show()
    '''



    # need to force the fit to go through zero
    popt, pcov = curve_fit(lin_fit_through_null, msd_time, msd)

    plt.figure()
    plt.errorbar(msd_time, msd, xerr = None, yerr = 2*d_msd, label="Simulated")
    plt.plot(msd_time, lin_fit_through_null(msd_time, *popt), label="Fit")


    #depends on dimensionality
    #plt.plot(msd_time, lin_fit_through_null(msd_time, 4.), label="Expected")
    plt.plot(msd_time, 4 * msd_time, label="Expected")
    #plt.plot(msd_time, 4 * D_eff * msd_time * v_0**2 / ( 2 * D_r**2 ) * ( 2 * D_r *
        #msd_time  +  np.exp(- 2. * D_r * msd_time) - 1 ), label = "Theory 3d") # for 3d
    plt.annotate("Fit Slope = {}".format(popt), (0.25, 0.5), textcoords="axes fraction")
    plt.xlabel("Time")
    plt.ylabel("MSD")
    plt.legend(loc='upper left')
    plt.show()

if __name__ == '__main__':
    main()
