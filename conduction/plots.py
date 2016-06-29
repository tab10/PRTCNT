import itertools
import logging
import glob
import matplotlib as mpl
import os
# mpl.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy import stats
import creation
# mpl.rcParams['text.usetex'] = True


def histogram_walker_2d_onlat(walker, grid_range, bins):
    """Takes walker instance and histograms how many times location is accessed over the simulation. H not normalized
    (for 1 walker only)"""
    # logging.info("Histogramming positions")
    walker_pos_xarr = np.zeros(len(walker.pos))
    walker_pos_yarr = np.zeros(len(walker.pos))
    for i in range(len(walker.pos)):
        walker_pos_xarr[i] = walker.pos[i][0]  # converts list of lists to list of arrays
        walker_pos_yarr[i] = walker.pos[i][1]
    H, xedges, yedges = np.histogram2d(walker_pos_xarr, walker_pos_yarr, range=grid_range, bins=bins)
    return H, xedges, yedges


def histogram_walker_3d_onlat(walker, grid_range, bins):
    """Takes walker instance and histograms how many times location is accessed over the simulation. H not normalized
    (for 1 walker only)"""
    H, edges = np.histogramdd(np.asarray(walker.pos), range=grid_range, bins=bins)
    x_edges = edges[0]
    y_edges = edges[1]
    z_edges = edges[2]
    return H, x_edges, y_edges, z_edges

def plot_two_d_random_walk_setup(tube_coords, grid_size, quiet, save_dir):
    """Plots setup and orientation of nanotubes"""
    logging.info("Plotting setup")
    colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])
    tube_x = []
    tube_y = []  # x and y coordinates of tubes unzipped
    for i in range(len(tube_coords)):
        tube_x.append(tube_coords[i][0::2])
        tube_y.append(tube_coords[i][1::2])
    for i in range(len(tube_x)):
        plt.plot(tube_x[i], tube_y[i], c=next(colors))
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


def plot_three_d_random_walk_setup(tube_coords, grid_size, quiet, save_dir):
    """Plots setup and orientation of nanotubes"""
    logging.info("Plotting setup")
    fig = plt.figure()
    ax = fig.add_subplot(111, projection='3d')
    colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])
    tube_x = []
    tube_y = []
    tube_z = []
    for i in range(len(tube_coords)):
        tube_x.append(tube_coords[i][0::3])
        tube_y.append(tube_coords[i][1::3])
        tube_z.append(tube_coords[i][2::3])
    for i in range(len(tube_x)):
        ax.plot(tube_x[i], tube_y[i], tube_z[i], c=next(colors))
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
        raise SystemExit
    plt.close()


def plot_histogram_walkers_2d_onlat(timesteps, H_tot, xedges, yedges, quiet, save_dir, gen_plots):
    """Plots temperature profile for all walkers"""
    logging.info("Plotting temperature profile")
    H_tot /= float(timesteps)  # normalization condition
    creation.check_for_folder(save_dir)
    np.savetxt('%s/temp.txt' % save_dir, H_tot, fmt='%.1E')
    temp_profile = H_tot
    if gen_plots:
        plt.title('Temperature density (dimensionless units)')
        X, Y = np.meshgrid(xedges, yedges)
        plt.pcolormesh(X, Y, temp_profile.T)  # transpose since pcolormesh reverses axes
        plt.xlabel('X')
        plt.ylabel('Y')
        plt.xlim(xedges[0], xedges[-1])
        plt.ylim(yedges[0], yedges[-1])
        plt.colorbar()
        plt.savefig('%s/temp.pdf' % save_dir)
        if not quiet:
            plt.show()
        plt.close()
    return temp_profile


def plot_linear_temp(temp_profile, quiet, save_dir, plots):
    test_mean = np.mean(temp_profile[1:], axis=1)
    test_std = np.std(temp_profile[1:], axis=1, ddof=1)
    x = np.arange(1, len(test_mean) + 1)
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, test_mean)
    gradient_err = np.mean(test_std)
    line = slope * x + intercept
    if plots == True:
        plt.errorbar(x, test_mean, yerr=test_std)
        plt.plot(x, line, color='black')
        plt.title("Temperature density (normalized by time) (dimensionless units)\n"
                  "$R^2=%.4f$, $\\frac{dT(x)}{dx}$=%.4E$\\pm$%.4E" % (r_value ** 2, slope, std_err))
        plt.xlabel('x')
        plt.ylabel('T(x)')
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


def plot_temp_gradient_2d_onlat(temp_profile, xedges, yedges, quiet, save_dir, gradient_cutoff):
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
    plt.xlim(xedges[0], xedges[-1])
    plt.ylim(yedges[0], yedges[-1])
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


def plot_k_convergence(quantity, quiet, save_dir):
    logging.info("Plotting k convergence")
    plt.plot(quantity)
    plt.title("Convergence of conductivity k")
    plt.xlabel('Timesteps')
    plt.ylabel('Conductivity k')
    plt.savefig('%s/k_convergence.pdf' % save_dir)
    if not quiet:
        plt.show()
    plt.close()


def plot_k_convergence_err(quantity, quiet, save_dir, begin_cov_check):
    logging.info("Plotting k convergence error")
    x = range(begin_cov_check, len(quantity) + begin_cov_check)
    plt.plot(x, quantity)
    plt.title("Error in convergence of conductivity k")
    plt.xlabel('Timesteps')
    plt.ylabel('Conductivity k error')
    plt.savefig('%s/k_convergence_err.pdf' % save_dir)
    if not quiet:
        plt.show()
    plt.close()


def plot_k_vs_num_tubes(tube_length, num_configs, grid_size, dim, exclude_vals):
    exclude_vals = map(str, exclude_vals)
    exclude_vals = [x + '_' for x in exclude_vals]
    folds = []
    zero_folds = []
    orientations = []
    for file in glob.glob("*_*_%d_*" % tube_length):
        checker = file.split('_')[0] + '_'
        if checker not in exclude_vals:
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
    for i in range(len(uni_orientations)):
        uni_tubes = len(sep_folds[i]) / num_configs
        uni_num_tubes = []
        for k in range(uni_tubes):
            uni_num_tubes.append(sep_folds[i][k * num_configs].split('_')[0])
        all_k_vals = np.zeros(len(sep_folds[i]))
        for j in range(len(sep_folds[i])):
            os.chdir(sep_folds[i][j])
            all_k_vals[j] = np.loadtxt('k.txt')
            os.chdir('..')
        k_vals = []
        k_err = []
        for l in range(len(uni_num_tubes)):
            k_vals.append(np.mean(all_k_vals[l * num_configs:(l + 1) * num_configs]))
            k_err.append(np.std(all_k_vals[l * num_configs:(l + 1) * num_configs], ddof=1) / np.sqrt(num_configs))
        plt.errorbar(uni_num_tubes, k_vals, yerr=k_err, fmt='o', label=uni_orientations[i])
    plt.title(
        'Tubes of length %d in a %dD cube of length %d\n%d configurations' % (tube_length, dim, grid_size, num_configs))
    plt.xlabel('Number of tubes')
    plt.ylabel('Conductivity k')
    plt.legend(loc=2)
    plt.tight_layout()
    plt.savefig('k_num_tubes_%d_%dD.pdf' % (tube_length, dim))
    plt.close()
