'''
TO DO:
    - other functions for getting position, orientation, and velocity


CONCERNS:
    - velocity vector may not have been totally neglected and I should just
      used orientation times the propulsion strength
'''
import copy
import h5py
import math
import numpy as np
import matplotlib.patches as patches
import matplotlib.pyplot as plt

def lin_fit(x, a, b):
    return a*x + b

def lin_fit_through_null(x, a):
    return a*x

def get_dimension(filename):
    H5          = h5py.File(filename, 'r')
    p_full      = H5['particles/all/position']
    position    = np.array(p_full['value'])
    dimension   = position.shape[2]
    H5.close()

    return dimension

def get_position_data(filename):
    # open and read data file
    H5 = h5py.File(filename, 'r')
    
    # particles/* data
    p_full      = H5['particles/all/position']
    p_step      = np.array(p_full['step'])
    p_time      = np.array(p_full['time'])
    position    = np.array(p_full['value'])

    H5.close()
    return p_step, p_time, position

def get_orientation_data(filename):
    # open and read data file
    H5 = h5py.File(filename, 'r')

    o_full      = H5['particles/all/orientation']
    o_step      = np.array(o_full['step'])
    o_time      = np.array(o_full['time'])
    orientation = np.array(o_full['value'])

    H5.close()
    return o_step, o_time, orientation

def get_velocity_data(filename):
    # open and read data file
    H5 = h5py.File(filename, 'r')

    v_full      = H5['particles/all/velocity']
    v_step      = np.array(v_full['step'])
    v_time      = np.array(v_full['time'])
    velocity    = np.array(v_full['value'])

    H5.close()
    return v_step, v_time, velocity

def get_msd_data(filename):
    # open and read data file
    H5 = h5py.File(filename, 'r')
    
    # dynamics/* data (different format than particles/*)    
    msd_full    = H5['dynamics/all/mean_square_displacement']
    msd_count   = np.array(msd_full['count']).flatten() # what is count?
    msd_time    = np.array(msd_full['time']).flatten()
    msd         = np.array(msd_full['value']).flatten()
    d_msd       = np.array(msd_full['error']).flatten()

    H5.close()

    # organize data into a single array so it can be sorted
    sort_column_index = 0 # index of colum by which data should be sorted (0-sorting by time, 1-sorting by value)

    msd_data = np.array((msd_time, msd, d_msd)).T 
    msd_time, msd, d_msd = msd_data[np.argsort(msd_data[:,sort_column_index])].T


    return msd_time, msd, d_msd

def get_ocf_data(filename):
    # open and read data file
    H5 = h5py.File(filename, 'r')
    
    # dynamics/* data (different format than particles/*)    
    ocf_full    = H5['dynamics/all/orientational_autocorrelation']
    ocf_count   = np.array(ocf_full['count']).flatten() # what is count?
    ocf_time    = np.array(ocf_full['time']).flatten()
    ocf         = np.array(ocf_full['value']).flatten()
    d_ocf       = np.array(ocf_full['error']).flatten()

    H5.close()

    # organize data into a single array so it can be sorted
    sort_column_index = 0 # index of colum by which data should be sorted (0-sorting by time, 1-sorting by value)

    ocf_data = np.array((ocf_time, ocf, d_ocf)).T 
    ocf_time, ocf, d_ocf = ocf_data[np.argsort(ocf_data[:,sort_column_index])].T


    return ocf_time, ocf, d_ocf

def get_density_mode_data(filename):
    # open and read data file
    H5 = h5py.File(filename, 'r')
    
    full        = H5['structure/all/density_mode']
    step        = np.array(full['step'])
    time        = np.array(full['time'])
    value       = np.array(full['value'])
    wavevector = np.array(full['wavevector'])

    H5.close()

    value = value[:,:,0] + value[:,:,1]*1j

    return step, time, value, wavevector



def get_v_0(filename):
    '''
    This function is a bit hacky and depends on the format of the .log file
    Is there a better way to get these parameter values? Directly from the .h5
    file?
    '''
    logfile = filename.split('.')[0] + ".log"

    with open(logfile, 'r') as f:
        for line in f:
            if "diffusion constants" in line:
                a = line[line.find("(")+2:line.find(")")].split(',')
                b = [float(i) for i in a]
    return b

def get_box_dims(filename):

    logfile = filename.split('.')[0] + ".log"

    with open(logfile, 'r') as f:
        for line in f:
            if "edge lengths of simulation domain" in line:
                a = line.split(':')[-1].split()
                b = [float(i) for i in a]
    return np.array(b)


def enforce_periodic_boundary_conditions(position, box):
    # Enforce periodic boundary conditions
    # should work for n-dimensions
    for p in position:
        for ii, particle in enumerate(p):
            for jj, axis in enumerate(particle):
                while abs(p[ii,jj]) >= box[jj] / 2.:
                    if p[ii, jj] > 0:
                        p[ii, jj] -= box[jj]
                    else:
                        p[ii, jj] += box[jj]
    return position

def get_cell_index(particle, box, cell_length):
    ''' for square/cubic cells and boxes, and cell indicies referring to a 1D array '''

    side = box[0]
    dims = len(box)
    # ncells is the cells per side
    ncells = math.floor( side / cell_length ) 

    if len(particle) != len(box):
        print "\nERROR: box and particle coordinates have different dimensions.\n"
        exit()
    if ncells != side / cell_length:
        print "\nERROR! the box length is not an exact multiple of the cell length!\n"
        exit()

    index = 0
    for i in range(len(box)):
        next_increment = math.floor( 1 / cell_length * ( particle[i] + 0.5
            * box[i]) ) * ncells ** i 
        index += next_increment

    return int(index)

def create_cell_linked_list(timestep, box, cell_length, dims):
    ''' not yet tested for 3d'''

    if dims == 3:
        print "Warning, create_cell_linked_list function not yet tested for 3d"
        exit()

    nparticles = len(timestep)
    ncells = int(math.floor( box[0] / cell_length ))
    cells = -np.ones( ncells ** dims , dtype = int )
    particles = -np.ones( nparticles , dtype = int )
    for ( p, particle ) in enumerate( timestep ):
        #for 1d cell linked-list
        cell_index           = get_cell_index(particle, box, cell_length)
        particles[p]         = cells[ cell_index ]
        cells[ cell_index ]  = p

    return cells, particles

def particle_and_velocity_binning(cells, particles, dims, orientation, cell_length, ncells):
    avg_v   = np.zeros( ( cells.shape[0], dims) )
    avg_rho = np.zeros( cells.shape )

    for (c, cell)  in enumerate( cells ):
        # next more recent particle added to the box
        nmr_part = cell
        nparticles_in_box = 0
        while nmr_part != -1:
            print nmr_part
            avg_v[c] += orientation[nmr_part]
            nparticles_in_box += 1
            nmr_part = particles[ nmr_part ]

        if nparticles_in_box != 0:
            avg_rho[c] = nparticles_in_box / cell_length ** dims
            avg_v[c]  /= nparticles_in_box

    avg_rho_reshaped =  avg_rho.reshape((ncells, ncells))
    avg_v =  avg_v.reshape((ncells, ncells, dims))

    for y in range(len(avg_rho_reshaped)):
        for x in range(len(avg_rho_reshaped[y])):
            if avg_rho_reshaped[y,x] != avg_rho[x + y * ncells]:
                print "ERROR: the reshaped avg_rho array doesn't match\
                up with the original"
                exit()


    return avg_rho_reshaped, avg_v

def plot_periodic_density_contour_2d(box, avg_rho, cell_length):
    xlist_pad = np.arange( -box[0] / 2. - cell_length - box[0],
            box[0] / 2. + box[0] + cell_length,  cell_length )
    ylist_pad = np.arange( -box[1] / 2. - cell_length - box[1],
            box[1] / 2. + box[1] + cell_length,  cell_length )
    X_pad, Y_pad = np.meshgrid(xlist_pad, ylist_pad)
    X_pad += 0.5 * cell_length
    Y_pad += 0.5 * cell_length

    big_rho = np.concatenate((avg_rho, avg_rho, avg_rho), axis = 0)
    big_rho = np.concatenate((big_rho, big_rho, big_rho), axis = 1)

    big_rho_pad = np.zeros(np.array(big_rho.shape) + 2)
    big_rho_pad[1:-1,1:-1] = big_rho
    plt.contourf(X_pad, Y_pad, big_rho_pad)
    cax = plt.axes([0.90, 0.1, 0.025, 0.8])
    plt.colorbar(cax=cax)



def plot_cells(box, cell_length):
    for i in range(int( box[0] / cell_length )):
        if i == 0:
            continue
        plt.axhline(i *cell_length - box[0] / 2 )
        plt.axvline(i *cell_length - box[0] / 2 )


