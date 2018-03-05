#!/usr/bin/python
'''
TO DO:
    - have theoretical displacement coefficents not hard-coded, formula for the
        dimensional dependence

    - might actually have to print out multiple theory curves if the data are
      for different parameters
    - decide what output I'd like (optimal fit parameters, etc.)
'''


from functions import *
import matplotlib.pyplot as plt
import numpy as np
import sys

def main(filenames):
    for filename in filenames:
        D_par, D_perp, D_rot, v_0 = get_v_0(filename)
        D_eff = ( D_par + 2 * D_perp ) / 3              # effective diffusion constant 
        msd_time, msd, d_msd = get_msd_data(filename)
        dims = get_dimension(filename)

        # only the first time around is the theoretically expected curve
        # plotted
        if filename == sys.argv[1]:
            # is D_eff correct here?
            # there might be an issue with the exponential term, it could
            # depend on the dimensionality
            plt.plot(msd_time, 2 * dims * D_eff * msd_time + v_0**2 / (2 *
                D_rot**2 ) * ( 2 * D_rot * msd_time + np.exp( - 2 * D_rot *
                    msd_time )  - 1 ), label="Theory")

        plt.errorbar(msd_time, msd, xerr = None, yerr = 2 * d_msd, label="Simulated")
        plt.title('Mean Square Displacement')
        plt.xlabel("Time")
        plt.ylabel("MSD")
        plt.xscale('log')
        plt.yscale('log')
        plt.legend(loc='upper left')
    plt.show()

if __name__ == '__main__':
    main(sys.argv[1:])
