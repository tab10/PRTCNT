import itertools
import logging
import glob
import matplotlib as mpl
import os

# mpl.use('Agg')
mpl.use('pdf')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy import stats
import creation
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
            int_x, int_y = zip(*grid.tube_squares[i])
            color = next(colors)
            plt.scatter(int_x[1:-1], int_y[1:-1], c=color)
            plt.scatter(tube_x[i], tube_y[i], c=color, marker=(5, 1))
    plt.xlim(0, grid_size)
    plt.ylim(0, grid_size)
    plt.title('Nanotube locations')
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.grid()
    creation.check_for_folder(save_dir)
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
            int_x, int_y, int_z = zip(*grid.tube_squares[i])
            color = next(colors)
            ax.scatter(int_x[1:-1], int_y[1:-1], int_z[1:-1], c=color)
            ax.scatter(tube_x[i], tube_y[i], tube_z[i], c=color, marker=(5, 1))
    ax.set_xlim(0, grid_size)
    ax.set_ylim(0, grid_size)
    ax.set_zlim(0, grid_size)
    ax.set_title('Nanotube locations')
    ax.set_xlabel('X')
    ax.set_ylabel('Y')
    ax.set_zlabel('Z')
    ax.grid()
    creation.check_for_folder(save_dir)
    plt.savefig('%s/setup.pdf' % save_dir)
    if not quiet:
        plt.show()
    plt.close()


def plot_check_array_2d(grid, quiet, save_dir, gen_plots):
    """Plots array with bd/vol locations"""
    logging.info("Plotting check array")
    creation.check_for_folder(save_dir)
    if gen_plots:
        if grid.tube_radius == 0:
            plt.pcolor(grid.tube_check_bd.T)  # since pcolor reverses axes
        else:
            plt.pcolor(grid.tube_check_bd_vol.T)  # since pcolor reverses axes
        plt.title("Nanotube setup, 0 nothing, 1 endpoint, -1 interior\n10 reflective, 20 periodic, 30 corner")
        plt.colorbar()
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.xlim(0, grid.size + 1)
        plt.ylim(0, grid.size + 1)
        plt.savefig('%s/setup_array.pdf' % save_dir)
        if not quiet:
            plt.show()
        plt.close()


def plot_histogram_walkers_onlat(grid, timesteps, H_tot, xedges, yedges, quiet, save_dir, gen_plots):
    """Plots temperature profile for all walkers"""
    logging.info("Plotting temperature profile")
    # H_tot /= float(timesteps)  # normalization condition
    creation.check_for_folder(save_dir)
    # np.savetxt('%s/temp.txt' % save_dir, H_tot, fmt='%.1E')
    temp_profile = H_tot
    if gen_plots:
        plt.title('Temperature density (dimensionless units)')
        # X, Y = np.meshgrid(xedges, yedges)
        plt.pcolor(temp_profile.T)  # transpose since pcolormesh reverses axes
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.xlim(0, grid.size + 1)
        plt.ylim(0, grid.size + 1)
        plt.colorbar()
        plt.savefig('%s/temp.pdf' % save_dir)
        if not quiet:
            plt.show()
        plt.close()
    return temp_profile


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
    creation.check_for_folder(save_dir)
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
    creation.check_for_folder(save_dir)
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
    creation.check_for_folder(save_dir)
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
    x_floor = range(2, 50)
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
        x = range(begin_cov_check, len(quantity) + begin_cov_check)
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


def plot_k_vs_num_tubes(tube_length, num_configs, grid_size, dim, legend=True, exclude_vals='',
                        tunneling=False, max_tube_num=100000):
    def fill_fraction_tubes(x, orientation, tunneling, grid_size, dim):
        random = {'0': 0, '1250': 1.78, '2500': 3.58, '3750': 5.36, '5000': 7.14, '6250': 8.92, '7500': 10.7,
                  '8750': 12.46, '10000': 14.27, '11250': 16.04, '12500': 17.8, '13750': 19.58}
        h_v = {'0': 0, '1250': 2.06, '2500': 4.12, '3750': 6.18, '5000': 8.24, '6250': 10.3, '7500': 12.36,
               '8750': 14.42, '10000': 16.48, '11250': 18.56, '12500': 20.62}
        tunnel = 2.0 * float(x) / grid_size ** dim
        if not tunneling:
            if orientation == 'random':
                fill_fract = random[str(int(x))]
            else:  # h or v
                fill_fract = h_v[str(int(x))]
        else:
            fill_fract = tunnel
        return fill_fract
    exclude_vals = map(str, exclude_vals)  # array of numbers
    exclude_vals = [x + '_' for x in exclude_vals]
    folds = []  # list of all folder name strings
    zero_folds = []
    orientations = []  # list of all orientations (not unique yet)
    for file in glob.glob("*_*_%d_*" % tube_length):
        checker = file.split('_')[0] + '_'
        config_num = int(file.split('_')[3])
        tube_val = int(file.split('_')[0])
        if (checker not in exclude_vals) and (config_num <= num_configs) and (
            tube_val <= max_tube_num):  # throws out extra config
            folds.append(file)  # all files
            orientations.append(file.split('_')[1])
    print folds
    for file in glob.glob("0_*_*_*"):
        zero_folds.append(file)
    uni_orientations = list(set(orientations))
    sep_folds = []
    # separate folds by orientation
    for i in range(len(uni_orientations)):
        sep_folds.append([x for x in folds if uni_orientations[i] in x])
        sep_folds[i] = sorted(sep_folds[i])
        sep_folds[i] += zero_folds
    for i in range(len(uni_orientations)):
        uni_tubes = len(sep_folds[i]) / num_configs
        uni_num_tubes = []
        for k in range(uni_tubes):
            uni_num_tubes.append(sep_folds[i][k * num_configs].split('_')[0])
        uni_num_tubes = [float(y) for y in uni_num_tubes]
        print uni_num_tubes
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
            k_err.append(np.std(all_k_vals[l * num_configs:(l + 1) * num_configs], ddof=1) / np.sqrt(num_configs))
        # plt.errorbar(uni_num_tubes, k_vals, yerr=k_err, fmt='o', label=uni_orientations[i])
        fill_fract = []
        for a in range(len(uni_num_tubes)):
            temp_ff = fill_fraction_tubes(uni_num_tubes[a], uni_orientations[i], tunneling, grid_size, dim)
            fill_fract.append(temp_ff)
        # apply linear fit
        slope, intercept, r_value, p_value, std_err = stats.linregress(fill_fract, k_vals)
        x_fit = np.linspace(min(fill_fract), max(fill_fract), num=50)
        y_fit = slope * x_fit + intercept
        plt.errorbar(fill_fract, k_vals, yerr=k_err, fmt='o', label=uni_orientations[i])
        fit_label = '%s, slope %.4E, y-int %.4E' % (uni_orientations[i], slope, intercept)
        plt.plot(x_fit, y_fit, label=fit_label)
    plt.title(
        'Tubes of length %d in a %dD cube of length %d\n%d configurations' % (
            tube_length, dim, grid_size, num_configs))
    # plt.xlabel('Number of tubes')
    plt.xlabel('Filling fraction %')
    plt.ylabel('Conductivity k')
    if legend:
        plt.legend(loc=2)
    plt.tight_layout()
    plt.savefig('k_num_tubes_%d_%dD.pdf' % (tube_length, dim))
    plt.close()


def plot_k_kapitza_fill_fract_side_by_side(kapitza_vals, kapitza_num_configs, tube_length, grid_size, dim,
                                           exclude_vals='', max_tube_num=100000):
    # kapitza_vals - a list with vals to use (should be from lowest p to highest)
    # kapitza_num_configs - a list in order with kapitza vals telling the configs to include
    exclude_vals = map(str, exclude_vals)  # array of numbers , ['0','1000','2000']
    exclude_vals = [x + '_' for x in exclude_vals]
    plt.figure()
    num_plots = len(kapitza_vals)
    for z in range(len(kapitza_vals)):
        # go into each directory and get this data
        os.chdir('%dD_kapitza_%s' % (dim, kapitza_vals[z]))
        num_configs = kapitza_num_configs[z]

        folds = []  # list of all folder name strings
        zero_folds = []
        orientations = []  # list of all orientations (not unique yet)
        for file in glob.glob("*_*_%d_*" % tube_length):
            checker = file.split('_')[0] + '_'
            config_num = int(file.split('_')[3])
            tube_val = int(file.split('_')[0])
            if (checker not in exclude_vals) and (config_num <= num_configs) and (
                        tube_val <= max_tube_num):  # throws out extra config
                folds.append(file)  # all files
                orientations.append(file.split('_')[1])
        print folds
        for file in glob.glob("0_*_*_*"):
            zero_folds.append(file)
        uni_orientations = list(set(orientations))
        sep_folds = []
        # separate folds by orientation
        for i in range(len(uni_orientations)):
            sep_folds.append([x for x in folds if uni_orientations[i] in x])
            sep_folds[i] = sorted(sep_folds[i])
            sep_folds[i] += zero_folds
        for i in range(len(uni_orientations)):
            uni_tubes = len(sep_folds[i]) / num_configs
            uni_num_tubes = []
            for k in range(uni_tubes):
                uni_num_tubes.append(sep_folds[i][k * num_configs].split('_')[0])
            uni_num_tubes = [float(y) for y in uni_num_tubes]
            print uni_num_tubes
            all_k_vals = np.zeros(len(sep_folds[i]))
            # all_kapitza_vals = np.zeros(len(sep_folds[i]))
            for j in range(len(sep_folds[i])):
                os.chdir(sep_folds[i][j])
                # kapitza = np.loadtxt('prob_m_cn.txt')
                all_k_vals[j] = np.loadtxt('k.txt')
                os.chdir('..')
            k_vals = []
            k_err = []
            # build filling fraction on x axis
            fill_fract = []  # will be in percentage
            for q in range(len(uni_num_tubes)):
                fill_fract.append((int(uni_num_tubes[q]) * tube_length * 100.0) / grid_size ** dim)
            for l in range(len(uni_num_tubes)):
                k_vals.append(np.mean(all_k_vals[l * num_configs:(l + 1) * num_configs]))
                k_err.append(np.std(all_k_vals[l * num_configs:(l + 1) * num_configs], ddof=1) / np.sqrt(num_configs))
            plt.subplot(num_plots, 1, z)
            plt.errorbar(fill_fract, k_vals, yerr=k_err, fmt='o', label=uni_orientations[i])
        # plt.title('Tubes of length %d in a %dD cube of length %d\n%d configurations' % (tube_length, dim, grid_size, num_configs))
        plt.xlabel('Filling fraction (%)')
        plt.ylabel('Conductivity k')
        # plt.ylim(0.0,)

        os.chdir('..')
    plt.legend()
    plt.tight_layout()
    plt.savefig('kapitza_k_plot.pdf')
    plt.close()
