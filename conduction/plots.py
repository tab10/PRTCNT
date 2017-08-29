# //////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////// ##  ##  ###  ## ### ### ///////////////////////////// #
# ////////////////////////////// # # # #  #  #   # #  #  ///////////////////////////// #
# ////////////////////////////// ##  ##   #  #   # #  #  ///////////////////////////// #
# ////////////////////////////// #   # #  #  #   # #  #  ///////////////////////////// #
# ////////////////////////////// #   # #  #   ## # #  #  ///////////////////////////// #
# ////////////////////////////// ###  #          ##           # ///////////////////////#
# //////////////////////////////  #      ###     # # # # ### ### ///////////////////// #
# //////////////////////////////  #   #  ###     ##  # # #    # ////////////////////// #
# //////////////////////////////  #   ## # #     # # ### #    ## ///////////////////// #
# //////////////////////////////  #              ## ////////////////////////////////// #
# //////////////////////////////////////////////////////////////////////////////////// #


from __future__ import division
import itertools
import logging
import glob
import matplotlib as mpl

mpl.use('pdf')
import os
# mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy import stats
import scipy as sp

from conduction import backend


# mpl.rcParams['text.usetex'] = True
# option causes problems on Schooner


def histogram_walker_2d_onlat(walker, grid_range, bins):
    """Takes walker instance and histograms how many times location is accessed over the simulation. H not normalized
    (for 1 walker only)
    FUNCTION IS OUTDATED, SHOULDN'T BE USED"""
    logging.info("Histogramming positions")
    walker_pos_xarr = np.zeros(len(walker.pos))
    walker_pos_yarr = np.zeros(len(walker.pos))
    for i in range(len(walker.pos)):
        walker_pos_xarr[i] = walker.pos[i][0]  # converts list of lists to list of arrays
        walker_pos_yarr[i] = walker.pos[i][1]
    H, xedges, yedges = np.histogram2d(walker_pos_xarr, walker_pos_yarr, range=grid_range, bins=bins)
    return H


def histogram_walker_3d_onlat(walker, grid_range, bins):
    """Takes walker instance and histograms how many times location is accessed over the simulation. H not normalized
    (for 1 walker only)
    FUNCTION IS OUTDATED, SHOULDN'T BE USED"""
    H, edges = np.histogramdd(np.asarray(walker.pos), range=grid_range, bins=bins)
    x_edges = edges[0]
    y_edges = edges[1]
    z_edges = edges[2]
    return H, x_edges, y_edges, z_edges


def plot_two_d_random_walk_setup(grid, quiet, save_dir):
    """Plots setup and orientation of nanotubes"""
    grid_size = grid.size
    tube_coords = grid.tube_coords
    tube_radius = grid.tube_radius
    colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y'])
    tube_x = []
    tube_y = []  # x and y coordinates of tubes unzipped
    for i in range(len(tube_coords)):
        tube_x.append(tube_coords[i][0::2])
        tube_y.append(tube_coords[i][1::2])
    if tube_radius == 0:
        logging.info("Plotting setup with no tube excluded volume")
        for i in range(len(tube_x)):
            plt.plot(tube_x[i], tube_y[i], c=next(colors))
            plt.scatter(tube_x[i], tube_y[i], c='black', marker=(5, 1))
    else:
        logging.info("Plotting setup with tube excluded volume")
        for i in range(len(grid.tube_squares)):
            int_x, int_y = list(zip(*grid.tube_squares[i]))
            color = next(colors)
            plt.scatter(int_x[1:-1], int_y[1:-1], c=color)
            plt.scatter(tube_x[i], tube_y[i], c=color, marker=(5, 1))
    plt.xlim(0, grid_size)
    plt.ylim(0, grid_size)
    plt.title('Nanotube locations')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.grid()
    backend.check_for_folder(save_dir)
    plt.savefig('%s/setup.pdf' % save_dir)
    if not quiet:
        plt.show()
    plt.close()


def plot_three_d_random_walk_setup(grid, quiet, save_dir):
    """Plots setup and orientation of nanotubes"""
    grid_size = grid.size
    tube_coords = grid.tube_coords
    tube_radius = grid.tube_radius
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y'])
    tube_x = []
    tube_y = []
    tube_z = []
    for i in range(len(tube_coords)):
        tube_x.append(tube_coords[i][0::3])
        tube_y.append(tube_coords[i][1::3])
        tube_z.append(tube_coords[i][2::3])
    if tube_radius == 0:
        logging.info("Plotting setup with no tube excluded volume")
        for i in range(len(tube_x)):
            ax.plot(tube_x[i], tube_y[i], tube_z[i], c=next(colors))
            ax.scatter(tube_x[i], tube_y[i], tube_z[i], c='black', marker=(5, 1))
    else:
        logging.info("Plotting setup with tube excluded volume")
        for i in range(len(grid.tube_squares)):
            int_x, int_y, int_z = list(zip(*grid.tube_squares[i]))
            color = next(colors)
            ax.scatter(int_x[1:-1], int_y[1:-1], int_z[1:-1], c=color, marker="s")
            ax.scatter(tube_x[i], tube_y[i], tube_z[i], c=color, marker=(5, 1))
    ax.set_xlim(0, grid_size)
    ax.set_ylim(0, grid_size)
    ax.set_zlim(0, grid_size)
    ax.set_title('Nanotube locations')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.grid()
    backend.check_for_folder(save_dir)
    plt.savefig('%s/setup.pdf' % save_dir)
    if not quiet:
        plt.show()
    plt.close()


def plot_check_array_2d(grid, quiet, save_dir, gen_plots):
    """Plots array with bd/vol locations"""
    logging.info("Plotting check array")
    backend.check_for_folder(save_dir)
    if gen_plots:
        if grid.tube_radius == 0:
            plt.pcolor(grid.tube_check_bd.T)  # since pcolor reverses axes
        else:
            plt.pcolor(grid.tube_check_bd_vol.T)  # since pcolor reverses axes
        plt.title("Nanotube setup, 0 nothing, 1 endpoint, -1 interior\n10 reflective, 20 periodic, 30 corner")
        plt.colorbar()
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.xlim(1, grid.size)
        plt.ylim(1, grid.size)
        plt.savefig('%s/setup_array.pdf' % save_dir)
        if not quiet:
            plt.show()
        plt.close()


def plot_colormap_2d(grid, H_tot, quiet, save_dir, gen_plots, title='Temperature density (dimensionless units)',
                     xlab='X', ylab='Y', filename='temp', random_slice=None, bds=False, vmin=None, vmax=None):
    """Plots temperature profile for all walkers
    Can be called anywhere a 2D colormap (of 2D data), basically a histogram, is needed"""
    logging.info("Plotting 2D temperature (histogram)")
    backend.check_for_folder(save_dir)
    # np.savetxt('%s/temp.txt' % save_dir, H_tot, fmt='%.1E')
    cushion = 5
    rand = np.random.randint(0, grid.size)
    if gen_plots:
        temp_profile = H_tot
        if random_slice == 1:
            temp_profile = H_tot[rand][:][:]  # YZ
        if random_slice == 2:
            temp_profile = H_tot[:][rand][:]  # YZ
        if random_slice == 3:
            temp_profile = H_tot[:][:][rand]  # YZ
        plt.title(title)
        # X, Y = np.meshgrid(xedges, yedges)
        if vmin is None:
            vmin = np.min(H_tot)
        if vmax is None:
            vmax = np.max(H_tot)
        plt.pcolor(temp_profile.T, vmin=vmin, vmax=vmax)  # transpose since pcolormesh reverses axes
        plt.xlabel(xlab)
        plt.ylabel(ylab)
        if not bds:
            plt.xlim(1, grid.size)
            plt.ylim(1, grid.size)
        else:
            plt.xlim(0, grid.size + 1)
            plt.ylim(0, grid.size + 1)
        plt.colorbar()
        plt.savefig('%s/%s.pdf' % (save_dir, filename))
        if not quiet:
            plt.show()
        plt.close()
    return temp_profile


def plot_bargraph_3d(grid, H_tot, x_edges, y_edges, quiet, save_dir, gen_plots,
                     title='Temperature density (dimensionless units)', xlab='X', ylab='Y', zlab='Z', filename='temp',
                     random_slice=None):
    """Plots histogram profile for all walkers
    Can be called anywhere a 3D bar-type histogram (of 2D data) is needed"""
    logging.info("Plotting 3D temperature (bar-type) histogram")
    if gen_plots:
        if random_slice:
            cushion = 5
            zax = np.random.randint(cushion, grid.size - cushion)
            H_temp = H_tot[:][:][zax]  # XY
            H_tot = H_temp
        backend.check_for_folder(save_dir)
        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')
        # Construct arrays for the anchor positions of bars
        xpos, ypos = np.meshgrid(x_edges[:-1], y_edges[:-1])  # +0.25 already added
        xpos = xpos.flatten('F')
        ypos = ypos.flatten('F')
        zpos = np.zeros_like(xpos)
        # Construct arrays with the dimensions for bars
        dx = 0.5 * np.ones_like(zpos)
        dy = dx.copy()
        dz = H_tot.flatten()
        # plot
        ax.bar3d(xpos, ypos, zpos, dx, dy, dz, color='b', zsort='average')
        ax.set_title(title)
        ax.set_xlabel(xlab)
        ax.set_ylabel(ylab)
        ax.set_zlabel(zlab)
        ax.set_xlim(1, grid.size)
        ax.set_ylim(1, grid.size)
        plt.savefig('%s/%s.pdf' % (save_dir, filename))
        if not quiet:
            plt.show()
        plt.close()
    return H_tot


def plot_linear_temp(temp_profile, grid_size, quiet, save_dir, plots, cutoff_frac=0.25):
    cutoff_dist = int(cutoff_frac * grid_size)
    # temp_profile sliced to remove
    test_mean = np.mean(temp_profile, axis=1)
    test_std = np.std(temp_profile, axis=1, ddof=1)
    x = np.arange(cutoff_dist, grid_size - cutoff_dist)
    x_full = np.arange(0, grid_size + 1)
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, test_mean[cutoff_dist:grid_size - cutoff_dist])
    gradient_err = np.mean(test_std)
    line = slope * x_full + intercept
    if plots == True:
        plt.errorbar(x_full, test_mean, yerr=test_std)
        plt.plot(x_full, line, color='black')
        plt.title("Temperature density (normalized by time) (dimensionless units)\n"
                  "$R^2=%.4f$, $\\frac{dT(x)}{dx}$=%.4E$\\pm$%.4E" % (r_value ** 2, slope, std_err))
        plt.xlabel('x')
        plt.ylabel('T(x)')
        plt.xlim(0, grid_size)
        plt.savefig('%s/temp_fit.pdf' % save_dir)
        if not quiet:
            plt.show()
        plt.close()
    return slope, std_err


def plot_walker_path_2d_onlat(walker, grid_size, temp, quiet, label, save_dir):
    """Plots path taken by a single walker"""
    logging.info("Plotting walker path")
    pos = walker.pos
    pos_x = np.zeros(len(pos))
    pos_y = np.zeros(len(pos))
    for i in range(len(pos)):
        pos_x[i] = pos[i][0]
        pos_y[i] = pos[i][1]
    plt.plot(pos_x, pos_y, zorder=1)
    plt.scatter(pos_x[0], pos_y[0], c='red', s=200, label='Start', zorder=2)
    plt.scatter(pos_x[-1], pos_y[-1], c='yellow', s=200, label='End', zorder=3)
    plt.legend()
    plt.xlim(0, grid_size)
    plt.ylim(0, grid_size)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.grid()
    plt.title('Walker %d path taken (%s)' % (label, temp))
    backend.check_for_folder(save_dir)
    plt.savefig('%s/%s_walker_%d_path.pdf' % (save_dir, temp, label))
    if not quiet:
        plt.show()
    plt.close()


def plot_walker_path_3d_onlat(walker, grid_size, temp, quiet, label, save_dir):
    """Plots path taken by a single walker 3D"""
    logging.info("Plotting walker path")
    pos = walker.pos
    pos_x = np.zeros(len(pos))
    pos_y = np.zeros(len(pos))
    pos_z = np.zeros(len(pos))
    for i in range(len(pos)):
        pos_x[i] = pos[i][0]
        pos_y[i] = pos[i][1]
        pos_z[i] = pos[i][2]
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    ax.set_xlim(0, grid_size)
    ax.set_ylim(0, grid_size)
    ax.set_zlim(0, grid_size)
    ax.set_title('Walker %d path taken (%s)' % (label, temp))
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.grid()
    ax.plot(pos_x, pos_y, pos_z, label="Path", zorder=1)
    ax.scatter(pos_x[0], pos_y[0], pos_z[0], 'o', c='red', s=200, label="Start", zorder=2)
    ax.scatter(pos_x[-1], pos_y[-1], pos_z[-1], 'o', c='yellow', s=200, label="End", zorder=3)
    plt.legend()
    backend.check_for_folder(save_dir)
    plt.savefig('%s/%s_walker_%d_path.pdf' % (save_dir, temp, label))
    if not quiet:
        plt.show()
    plt.close()


def plot_temp_gradient_2d_onlat(grid, temp_profile, xedges, yedges, quiet, save_dir, gradient_cutoff):
    """Plots temperature gradient for all walkers"""
    logging.info("Plotting temperature gradient")
    # disregard first few x= slices as close to the wall and values have large errors
    # gradient_cutoff = 2 seems to work well. Don't change.
    temp_gradient_x, temp_gradient_y = np.gradient(temp_profile)
    np.savetxt('%s/temp_gradient.txt' % save_dir, temp_gradient_x, fmt='%.1E')
    plt.title('Temperature gradient $\\frac{dT(x)}{dx}$ (dimensionless units)\nExcluding first %d x=bins'
              % gradient_cutoff)
    # Y gradient is irrelevant since we have periodic bd conditions in that direction
    X, Y = np.meshgrid(xedges[gradient_cutoff:], yedges)
    plt.pcolormesh(X, Y, temp_gradient_x[gradient_cutoff:].T)  # transpose since pcolormesh reverses axes
    plt.xlabel('X')
    plt.ylabel('Y')
    # plt.xlim(xedges[0], xedges[-1])
    plt.xlim(0, grid.size + 1)
    plt.ylim(0, grid.size + 1)
    plt.colorbar()
    backend.check_for_folder(save_dir)
    plt.savefig('%s/temp_gradient.pdf' % save_dir)
    if not quiet:
        plt.show()
    plt.close()
    return temp_gradient_x


def plot_check_gradient_noise_floor(temp_gradient_x, quiet, save_dir):
    logging.info("Plotting temperature gradient noise floor")
    conv = []
    size = len(temp_gradient_x)
    for i in range(size):
        conv.append(float(np.mean(temp_gradient_x[i:])))
    d_conv = []
    for i in range(len(conv) - 1):
        temp = np.abs(conv[i + 1] - conv[i])
        d_conv.append(temp)
    noise_floor_vals = d_conv[2:size / 2]  # these are rough bounds that should work
    plt.plot(d_conv)
    noise_floor = np.mean(noise_floor_vals)
    x_floor = list(range(2, 50))
    y_floor = [noise_floor for number in range(2, 50)]
    plt.plot(x_floor, y_floor)
    plt.title("Convergence of $\\frac{dT(x)}{dx}$ calculation sweeping from x to the end\n"
              "Noise floor: %.4E" % noise_floor)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.savefig('%s/temp_gradient_conv.pdf' % save_dir)
    if not quiet:
        plt.show()
    plt.close()


# def plot_distribution_tube_lengths


# def plot_distribution_tube_angles


# def plot_distribution_tube_distances(centers or ends)



def plot_k_convergence(quantity, quiet, save_dir, x_list=None):
    logging.info("Plotting k convergence")
    if x_list is not None:
        plt.plot(x_list, quantity)
        plt.xlabel('Timesteps')
    else:
        plt.plot(quantity)
        plt.xlabel('Total walkers/2')
    plt.title("Convergence of conductivity k")
    plt.ylabel('Conductivity k')
    plt.savefig('%s/k_convergence.pdf' % save_dir)
    if not quiet:
        plt.show()
    plt.close()


def plot_k_convergence_err(quantity, quiet, save_dir, begin_cov_check, x_list=None):
    logging.info("Plotting k convergence error")
    if x_list is not None:
        plt.plot(x_list, quantity)
        plt.xlabel('Timesteps')
    else:
        x = list(range(begin_cov_check, len(quantity) + begin_cov_check))
        plt.plot(x, quantity)
        plt.xlabel('Total walkers/2')
    plt.title("Error in convergence of conductivity k")
    plt.ylabel('Conductivity k error')
    plt.savefig('%s/k_convergence_err.pdf' % save_dir)
    if not quiet:
        plt.show()
    plt.close()


def plot_dt_dx(quantity, quiet, save_dir, x_list=None):
    logging.info("Plotting dt/dx")
    if x_list is not None:
        plt.plot(x_list, quantity)
        plt.xlabel('Timesteps')
    else:
        plt.plot(quantity)
        plt.xlabel('Total walkers/2')
    plt.title("dT(x)/dx")
    plt.ylabel('dT(x)/dx')
    plt.savefig('%s/dt_dx.pdf' % save_dir)
    if not quiet:
        plt.show()
    plt.close()


def plot_heat_flux(quantity, quiet, save_dir, x_list=None):
    logging.info("Plotting heat flux")
    if x_list is not None:
        plt.plot(x_list, quantity)
        plt.xlabel('Timesteps')
    else:
        plt.plot(quantity)
        plt.xlabel('Total walkers/2')
    plt.title("Heat flux")
    plt.ylabel('Heat flux')
    plt.savefig('%s/heat_flux.pdf' % save_dir)
    if not quiet:
        plt.show()
    plt.close()


def plot_k_vs_num_tubes(tube_length, num_configs, grid_size, dim, div_by_k0=True, legend=True, exclude_vals='',
                        tunneling=False, max_tube_num=100000, force_y_int=False, y_max=None,
                        units=False, dec_fill_fract=False, w_err=True):
    """w_err - weighted linear fit based on k error bars from configurations"""
    def fill_fraction_tubes(x, orientation, tunneling, grid_size, dim):
        ######
        random_2d_15_100 = {'0': 0, '10': 2.04, '20': 4.04, '30': 6.15, '40': 8.22, '50': 10.06, '60': 12.21,
                            '70': 14.23,
                        '80': 16.25, '90': 17.97, '100': 20.2, '110': 22.49, '120': 24.37, '130': 26.08,
                        '140': 28.25, '150': 30.34}
        ######
        h_v_2d_15_100 = {'0': 0, '10': 1.63, '20': 3.26, '30': 4.9, '40': 6.53, '50': 8.16, '60': 9.79, '70': 11.43,
                     '80': 13.06, '90': 14.69, '100': 16.32, '110': 17.96, '120': 19.59, '130': 21.22, '140': 22.85,
                     '150': 24.49}
        ###### CORRECTED 8/6/2017 DUE TO NEW RANDOM TUBE GENERATION ALGORITHM ######
        random_3d_15_100 = {'0': 0, '1250': 2.91, '2500': 5.76, '3750': 8.63, '5000': 11.43, '6250': 14.15,
                            '7500': 16.93,
                        '8750': 19.6, '10000': 22.29, '11250': 24.92, '12500': 27.49, '13750': None}
        ######
        h_v_3d_15_100 = {'0': 0, '1250': 2.06, '2500': 4.12, '3750': 6.18, '5000': 8.24, '6250': 10.3, '7500': 12.36,
               '8750': 14.42, '10000': 16.48, '11250': 18.56, '12500': 20.62}
        ######
        random_2d_10_100 = {'0': 0, '14': 1.95, '30': 4.2, '45': 6.29, '58': 8.23, '75': 10.58, '88': 12.32,
                            '102': 14.19,
                        '115': 16.05, '131': 18.48, '143': 19.99, '156': 21.81, '172': 24.14, '191': 26.47,
                        '201': 28.13, '220': 30.67}
        ######
        h_v_2d_10_100 = {'0': 0, '15': 1.63, '29': 3.26, '44': 4.9, '58': 6.53, '73': 8.16, '87': 9.79, '102': 11.43,
                     '116': 13.06, '131': 14.69, '145': 16.32, '160': 17.96, '175': 19.59, '189': 21.22, '204': 22.85}
        ######
        random_3d_10_100 = {'0': 0, '1794': 2.85, '3609': 5.7, '5403': 8.5, '7212': 11.3, '8992': 14.1, '10801': 16.8,
                        '12604': 19.6, '14429': 22.3, '16228': 24.97, '18060': 27.7, '19845': 30.25}
        ######
        h_v_3d_10_100 = {'0': 0, '1813': 2.06, '3640': 4.12, '5453': 6.18, '7266': 8.24, '9080': 10.3, '10906': 12.36,
                     '12719': 14.42, '14533': 16.48, '16376': 18.56, '18190': 20.62}
        ######
        random_2d_20_100 = {'0': 0, '8': 2.02, '15': 4.09, '23': 6.25, '30': 8.02, '39': 10.55, '46': 12.14,
                            '53': 14.28,
                        '60': 15.97, '70': 18.34, '75': 19.73, '82': 21.4, '89': 23.21, '99': 25.97,
                        '108': 28.15, '116': 30.64}
        ######
        h_v_2d_20_100 = {'0': 0, '8': 1.63, '15': 3.26, '23': 4.9, '30': 6.53, '38': 8.16, '46': 9.79, '53': 11.43,
                     '61': 13.06, '69': 14.69, '76': 16.32, '84': 17.96, '91': 19.59, '99': 21.22, '107': 22.85}
        ######
        random_3d_20_100 = {'0': 0, '951': 2.91, '1904': 5.73, '2863': 8.6, '3822': 11.35, '4757': 14.07, '5728': 16.8,
                        '6660': 19.35, '7627': 21.96, '8638': 24.65, '9562': 27.1, '10519': 29.5}
        ######
        h_v_3d_20_100 = {'0': 0, '950': 2.06, '1904': 4.12, '2854': 6.18, '3808': 8.24, '4758': 10.3, '5712': 12.36,
                     '6662': 14.42, '7616': 16.48, '8575': 18.56, '9529': 20.62}
        ######
        h_v_3d_20_150 = {'0': 0, '3476': 2.06, '6953': 4.12, '10429': 6.18, '13905': 8.24, '17381': 10.3,
                         '20858': 12.36,
                         '24334': 14.42, '27810': 16.48, '31320': 18.56, '34796': 20.62}
        ######
        tunnel = 2.0 * float(x) * 100.0 / grid_size ** dim
        if not tunneling:
            if orientation == 'random':
                search_str = 'random_%dd_%d_%d[str(int(x))]' % (dim, tube_length, grid_size)
                fill_fract = eval(search_str)
            elif (orientation == 'horizontal') or (orientation == 'vertical'):
                search_str = 'h_v_%dd_%d_%d[str(int(x))]' % (dim, tube_length, grid_size)
                fill_fract = eval(search_str)
        else:
            fill_fract = tunnel
        return fill_fract

    def lin_fit(x, y, dim):
        '''Fits a linear fit of the form mx+b to the data'''
        dim_dict = {2: 0.5, 3: 1.0 / 3.0}
        fitfunc = lambda params, x: params[0] * x + dim_dict[dim]  # create fitting function of form mx+no_tubes_const
        errfunc = lambda p, x, y: fitfunc(p, x) - y  # create error function for least squares fit

        init_a = 0.5  # find initial value for a (gradient)
        init_p = np.array((init_a))  # bundle initial values in initial parameters

        # calculate best fitting parameters (i.e. m and b) using the error function
        p1, success = sp.optimize.leastsq(errfunc, init_p.copy(), args=(x, y))
        f = fitfunc(p1, x)  # create a fit with those parameters
        return p1, f

    scaling_factor = 87.0  # FOR UNITS
    exclude_vals = list(map(str, exclude_vals))  # array of numbers
    exclude_vals = [x + '_' for x in exclude_vals]
    folds = []  # list of all folder name strings
    zero_folds = []
    orientations = []  # list of all orientations (not unique yet)
    dim = int(dim)
    tube_length = int(tube_length)
    old_plot = 'k_num_tubes_%d_%dD.pdf' % (tube_length, dim)  # let's get rid of the old one!
    if os.path.isfile(old_plot):
        os.remove(old_plot)
    for file in glob.glob("*_*_%d_*" % tube_length):
        checker = file.split('_')[0] + '_'
        config_num = int(file.split('_')[3])
        tube_val = int(file.split('_')[0])
        if (checker not in exclude_vals) and (config_num <= num_configs) and (
            tube_val <= max_tube_num):  # throws out extra config
            folds.append(file)  # all files
            orientations.append(file.split('_')[1])
    for file in glob.glob("0_*_*_*"):
        zero_folds.append(file)
    uni_orientations = list(set(orientations))
    sep_folds = []
    # separate folds by orientation
    for i in range(len(uni_orientations)):
        sep_folds.append([x for x in folds if uni_orientations[i] in x])
        sep_folds[i] = sorted(sep_folds[i])
        sep_folds[i] += zero_folds
    slopes = []
    d_slopes = []  # error on the slope
    y_ints = []
    r_twos = []
    for i in range(len(uni_orientations)):
        uni_tubes = int(len(sep_folds[i]) / num_configs)
        uni_num_tubes = []
        for k in range(uni_tubes):
            uni_num_tubes.append(sep_folds[i][k * num_configs].split('_')[0])
        uni_num_tubes = [float(y) for y in uni_num_tubes]
        all_k_vals = np.zeros(len(sep_folds[i]))
        # all_kapitza_vals = np.zeros(len(sep_folds[i]))
        for j in range(len(sep_folds[i])):
            os.chdir(sep_folds[i][j])
            # kapitza = np.loadtxt('prob_m_cn.txt')
            all_k_vals[j] = np.loadtxt('k.txt')
            os.chdir('..')
        k_vals = []
        k_err = []
        for l in range(len(uni_num_tubes)):
            k_vals.append(np.mean(all_k_vals[l * num_configs:(l + 1) * num_configs]))
            k_err.append(
                np.std(all_k_vals[l * num_configs:(l + 1) * num_configs], ddof=1) / np.sqrt(num_configs))
        fill_fract = []
        for a in range(len(uni_num_tubes)):
            temp_ff = fill_fraction_tubes(uni_num_tubes[a], uni_orientations[i], tunneling, grid_size, dim)
            fill_fract.append(temp_ff)
        # sort data ascending
        fill_fract_temp = np.array(fill_fract)
        k_err_temp = np.array(k_err)
        k_vals_temp = np.array(k_vals)
        idx = np.argsort(fill_fract)
        fill_fract = fill_fract_temp[idx]
        k_err = k_err_temp[idx]
        k_vals = k_vals_temp[idx]
        # sort should be working
        # remove duplicate fill fractions
        unq, unq_idx = np.unique(fill_fract, return_index=True)
        fill_fract = fill_fract[unq_idx]
        k_vals = k_vals[unq_idx]
        k_err = k_err[unq_idx]
        # divide by k0?
        if div_by_k0:
            k0 = k_vals[0]
            k_vals = k_vals / k_vals[0]
        # add units (if requested)
        if units:
            k_vals *= scaling_factor  # adds W/m*K units based on polystyrene matrix as measured by Matt Houck
            k_err *= scaling_factor
        if dec_fill_fract:
            fill_fract *= 0.01  # uses decimal for fill fraction values, more reasonable slopes
        # apply linear fit
        if force_y_int:
            dim_dict = {2: 0.5, 3: 1.0 / 300.0}
            slope, _ = lin_fit(fill_fract, k_vals, dim)
            print(slope)
            raise SystemExit
            x = np.array(fill_fract)
            y = np.array(k_vals)
            # x = x[:, np.newaxis]  # for 0 y intercept
            intercept = dim_dict[dim]
            x = np.vstack([x, np.ones(len(x)) * intercept]).T  # forces set y-int
            a, _, _, _ = np.linalg.lstsq(x, y)
            slope = a[0]
            r_value = a[1]
            x_fit = x
            y_fit = slope * x
        else:
            if w_err:
                k_err_weights = 1.0 / k_err
                p, V = np.polyfit(fill_fract, k_vals, 1, cov=True, w=k_err_weights)
            else:
                p, V = np.polyfit(fill_fract, k_vals, 1, cov=True, w=k_err)
            slope = p[0]
            intercept = p[1]
            d_slope = np.sqrt(V[0][0])
            d_yint = np.sqrt(V[1][1])
            notused_slope, notused_intercept, r_value, p_value, std_err = stats.linregress(fill_fract, k_vals)
            x_fit = np.linspace(min(fill_fract), max(fill_fract), num=50)
            y_fit = slope * x_fit + intercept
            # d_slope = np.abs(slope) * np.sqrt(((1 / r_value ** 2) - 1) / (num_configs - 2))
        slopes.append(slope)
        y_ints.append(intercept)
        r_twos.append(r_value ** 2)
        d_slopes.append(d_slope)
        plt.errorbar(fill_fract, k_vals, yerr=k_err, fmt='o', label=uni_orientations[i])
        fit_label = '%s, slope %.4E, y-int %.4E' % (uni_orientations[i], slope, intercept)
        plt.plot(x_fit, y_fit)  # , label=fit_label)
        g = open("%s_data.txt" % uni_orientations[i], 'w')
        g.write('fill_fract k k_err\n')
        for i in range(len(fill_fract)):
            header = '%.4E %.4E %.4E\n' % (fill_fract[i], k_vals[i], k_err[i])
            g.write(header)
        g.close()
    plt.title(
        'Tubes of length %d in a %dD cube of length %d\n%d configurations' % (
            tube_length, dim, grid_size, num_configs))
    # plt.xlabel('Number of tubes')
    if dec_fill_fract:
        plt.xlabel('Volume fraction')
    else:
        plt.xlabel('Volume fraction %')
    if div_by_k0:
        plt.ylabel('Thermal conductivity $k/k_0$')
    else:
        if units:
            plt.ylabel('Thermal conductivity k (W/(M*K))')
        else:
            plt.ylabel('Thermal conductivity k')
    if legend:
        plt.legend(loc=2)
    plt.tight_layout()
    if y_max:
        plt.ylim((0, y_max))
    else:
        plt.ylim(ymin=0)
    plt.savefig('k_num_tubes_%d_%dD.pdf' % (tube_length, dim))
    f = open("fit.txt", 'w')
    f.write('orientation slope d_slope y_int r_twos\n')
    for i in range(len(uni_orientations)):
        header = '%s %.4E %.4E %.4E %.4E\n' % (uni_orientations[i], slopes[i], d_slopes[i], y_ints[i], r_twos[i])
        f.write(header)
    f.close()
    plt.close()


def plot_all_tube_length(tube_length, num_configs, grid_size, div_by_k0=True, y_max=None, units=False,
                         dec_fill_fract=False, w_err=True):
    # Runs plot analysis in all folders for a tube length
    # div_by_k0 True or False divides the k values by k0 at 0 fill fraction
    os.chdir('tube_length_%d' % tube_length)
    for file in glob.glob("*_*"):
        if 'kapitza' in file:
            dim = file.split('_')[0]
            prob = file.split('_')[2]
            os.chdir(file)
            print('In directory %s' % file)
            plot_k_vs_num_tubes(tube_length, num_configs, grid_size, dim, div_by_k0=div_by_k0, units=units,
                                dec_fill_fract=dec_fill_fract, w_err=w_err)
            os.chdir('..')
        elif ('tunneling' in file) and ('novol' not in file):
            dim = file.split('_')[0]
            os.chdir(file)
            print('In directory %s' % file)
            plot_k_vs_num_tubes(tube_length, num_configs, grid_size, dim, div_by_k0=div_by_k0,
                                units=units, dec_fill_fract=dec_fill_fract, w_err=w_err)
            os.chdir('..')
        elif ('tunneling' in file) and ('novol' in file):
            dim = file.split('_')[0]
            os.chdir(file)
            print('In directory %s' % file)
            plot_k_vs_num_tubes(tube_length, num_configs, grid_size, dim, div_by_k0=div_by_k0,
                                units=units, dec_fill_fract=dec_fill_fract, tunneling=True, w_err=w_err)
            os.chdir('..')
    os.chdir('..')
    print('Done!')


def plot_multi_tube_lengths(tube_lengths_str, type, num_configs=5, mark_size='6_7.5_9',
                            leg_loc=2, leg_size=9, div_by_k0=True, plot_fits=True):
    # plot multiple tube lengths on same plot for each folder
    # one type only, like kapitzaX
    # tube_lengths_str separated with _
    # leg_loc and leg_size define legend location and size
    # errorbar size for the 3 sizes given above
    mark_size_spl = mark_size.split('_')
    color_list = {'10': "r", "15": "g", '20': "b"}
    marker_list = {'horizontal': "^", "vertical": "v", 'random': "o"}
    size_list = {'10': float(mark_size_spl[0]), "15": float(mark_size_spl[1]), '20': float(mark_size_spl[2])}
    tube_lengths_str = tube_lengths_str.split('_')
    for i in range(len(tube_lengths_str)):
        tube_lengths_str[i] = int(tube_lengths_str[i])
        os.chdir('tube_length_%d' % tube_lengths_str[i])
        dirs = []  # directories to check in all tube length folders
        in_dirs = []
        for file in glob.glob("*_*"):
            if os.path.isdir(file) and type == file:
                os.chdir(file)
                dirs.append(file)
                data_files = glob.glob("*_data.txt")
                for j in range(len(data_files)):
                    in_dirs.append(data_files[j])
                    f = open(data_files[j], 'r')
                    lines = f.readlines()[1:]
                    f.close()
                    fill_fract = []
                    k_vals = []
                    k_err = []
                    uni_orientations = []
                    slopes = []
                    d_slopes = []
                    y_ints = []
                    r_twos = []
                    for k in range(len(lines)):
                        fill_fract_t, k_vals_t, k_err_t = lines[k].split(" ")
                        fill_fract.append(float(fill_fract_t))
                        k_vals.append(float(k_vals_t))
                        k_err.append(float(k_err_t[:-1]))  # -1 for the newline removal
                    # import fit
                    g = open('fit.txt', 'r')
                    lines = g.readlines()[1:]
                    for l in range(len(lines)):
                        uni_orientations_t, slopes_t, d_slopes_t, y_ints_t, r_twos_t = lines[l].split(" ")
                        uni_orientations.append(uni_orientations_t)
                        slopes.append(float(slopes_t))
                        d_slopes.append(float(d_slopes_t))
                        y_ints.append(float(y_ints_t))
                        r_twos.append(float(r_twos_t[:-1]))  # -1 for the newline removal
                    g.close()
                    orientation_str = data_files[j].replace('_data.txt', '')
                    model_str = file.replace('_', ' ')
                    tube_l_str = tube_lengths_str[i]
                    for m in range(len(uni_orientations)):
                        if orientation_str == uni_orientations[m]:
                            x_fit = np.linspace(min(fill_fract), max(fill_fract), num=50)
                            y_fit = slopes[m] * x_fit + y_ints[m]
                    # let's plot here
                    legend_label = 'Tube length %d %s %s' % (tube_l_str, model_str, orientation_str)
                    print('Plotting %s' % legend_label)
                    plt.errorbar(fill_fract, k_vals, yerr=k_err, fmt=marker_list[orientation_str],
                                 c=color_list[str(tube_l_str)], label=legend_label,
                                 markersize=size_list[str(tube_l_str)])
                    if plot_fits:
                        plt.plot(x_fit, y_fit, c=color_list[str(tube_l_str)], linewidth=0.25)
                os.chdir('..')
        os.chdir('..')
    plt.legend(loc=leg_loc, prop={'size': leg_size})
    plt.title('Thermal conductivity vs. filling fraction percentage\nTubes of length 10, 15, '
              '20 with different orientations, %d configurations\nModel: %s' % (num_configs, model_str))
    plt.xlabel('Filling fraction %')
    if div_by_k0:
        plt.ylabel('Thermal conductivity $k/k_0$ (dimensionless units)')
    else:
        plt.ylabel('Thermal conductivity k (dimensionless units)')
    plt.ylim(ymin=0)
    plt.xlim(xmax=22)
    plt.tight_layout()
    plt.savefig('%s_%s_plot.pdf' % (tube_lengths_str, type))
    plt.close()


def plot_slopes_bar_graph(tube_lengths_str, type, num_configs=5, mark_size='6_7.5_9',
                          leg_loc=2, leg_size=9, div_by_k0=True, plot_fits=True):
    # plot multiple tube lengths on same plot for each folder
    # one type only, like kapitzaX
    # tube_lengths_str separated with _
    # leg_loc and leg_size define legend location and size
    # errorbar size for the 3 sizes given above
    plot_list = ['random_func', 'random_nofunc', 'horizontal_func', 'horizontal_nofunc', 'vertical_func',
                 'vertical_nofunc']
    mark_size_spl = mark_size.split('_')
    color_list = {'10': "r", "15": "g", '20': "b"}

    size_list = {'10': float(mark_size_spl[0]), "15": float(mark_size_spl[1]), '20': float(mark_size_spl[2])}
    tube_lengths_str = tube_lengths_str.split('_')
    for i in range(len(tube_lengths_str)):
        tube_lengths_str[i] = int(tube_lengths_str[i])
        os.chdir('tube_length_%d' % tube_lengths_str[i])
        dirs = []  # directories to check in all tube length folders
        in_dirs = []
        for file in glob.glob("*_*"):
            if os.path.isdir(file) and type == file:
                os.chdir(file)
                dirs.append(file)
                data_files = glob.glob("*_data.txt")
                for j in range(len(data_files)):
                    in_dirs.append(data_files[j])
                    f = open(data_files[j], 'r')
                    lines = f.readlines()[1:]
                    f.close()
                    fill_fract = []
                    k_vals = []
                    k_err = []
                    uni_orientations = []
                    slopes = []
                    d_slopes = []
                    y_ints = []
                    r_twos = []
                    for k in range(len(lines)):
                        fill_fract_t, k_vals_t, k_err_t = lines[k].split(" ")
                        fill_fract.append(float(fill_fract_t))
                        k_vals.append(float(k_vals_t))
                        k_err.append(float(k_err_t[:-1]))  # -1 for the newline removal
                    # import fit
                    g = open('fit.txt', 'r')
                    lines = g.readlines()[1:]
                    for l in range(len(lines)):
                        uni_orientations_t, slopes_t, d_slopes_t, y_ints_t, r_twos_t = lines[l].split(" ")
                        uni_orientations.append(uni_orientations_t)
                        slopes.append(float(slopes_t))
                        d_slopes.append(float(d_slopes_t))
                        y_ints.append(float(y_ints_t))
                        r_twos.append(float(r_twos_t[:-1]))  # -1 for the newline removal
                    g.close()
                    orientation_str = data_files[j].replace('_data.txt', '')
                    model_str = file.replace('_', ' ')
                    tube_l_str = tube_lengths_str[i]
                    for m in range(len(uni_orientations)):
                        if orientation_str == uni_orientations[m]:
                            x_fit = np.linspace(min(fill_fract), max(fill_fract), num=50)
                            y_fit = slopes[m] * x_fit + y_ints[m]
                    # let's plot here
                    legend_label = 'Tube length %d %s %s' % (tube_l_str, model_str, orientation_str)
                    print('Plotting %s' % legend_label)
                    plt.errorbar(fill_fract, k_vals, yerr=k_err, fmt=marker_list[orientation_str],
                                 c=color_list[str(tube_l_str)], label=legend_label,
                                 markersize=size_list[str(tube_l_str)])
                    if plot_fits:
                        plt.plot(x_fit, y_fit, c=color_list[str(tube_l_str)], linewidth=0.25)
                os.chdir('..')
        os.chdir('..')
    plt.legend(loc=leg_loc, prop={'size': leg_size})
    plt.title('Thermal conductivity vs. filling fraction percentage\nTubes of length 10, 15, '
              '20 with different orientations, %d configurations\nModel: %s' % (num_configs, model_str))
    plt.xlabel('Filling fraction %')
    if div_by_k0:
        plt.ylabel('Thermal conductivity $k/k_0$ (dimensionless units)')
    else:
        plt.ylabel('Thermal conductivity k (dimensionless units)')
    plt.ylim(ymin=0)
    plt.xlim(xmax=22)
    plt.tight_layout()
    plt.savefig('%s_%s_plot.pdf' % (tube_lengths_str, type))
    plt.close()
