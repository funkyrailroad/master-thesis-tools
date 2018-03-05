#!/usr/bin/python
'''
TO DO:
    - have theoretical displacement coefficents not hard-coded, formula for the
        dimensional dependence
    - only plot theory curve once

'''


from functions import *
import matplotlib.pyplot as plt
import numpy as np
import sys

def main(filenames):
    for filename in filenames:
        D_par, D_perp, D_rot, v_0 = get_v_0(filename)
        D_eff = ( D_par + 2 * D_perp ) / 3              # effective diffusion constant 
        ocf_time, ocf, d_ocf = get_ocf_data(filename)
        dims = get_dimension(filename)

        # only the first time around is the theoretically expected curve plotted
        if filename == sys.argv[1]:
            # is D_rot correct here? 
            if dims == 2:
                plt.plot(ocf_time, np.exp(- 1. * D_rot * ocf_time), label =
                        "2d Theory $e^{-Dt}$") 
            if dims == 3:
                plt.plot(ocf_time, np.exp(- 2. * D_rot * ocf_time), label =
                        "3d Theory $e^{-2Dt}$") 
    
        plt.plot(ocf_time, ocf, label = filename)
        plt.title('Orientational Correlation Function')
        plt.xlabel('Time')
        plt.ylabel('OCF')
        plt.xscale('log')
        plt.yscale('linear')
        plt.legend(loc='upper right')
    plt.show()

if __name__ == '__main__':
    main(sys.argv[1:])
