import itertools
import matplotlib.pyplot as plt
import matplotlib
import numpy as np
import logging
import creation

matplotlib.rcParams['text.usetex'] = True


def histogram_walker_2d_onlat(walker, grid_range, temp, bins):
    """Takes walker instance and histograms how many times location is accessed over the simulation. H not normalized
    (for 1 walker only)"""
    logging.info("Histogramming positions")
    walker_pos_xarr = np.zeros(len(walker.pos))
    walker_pos_yarr = np.zeros(len(walker.pos))
    for i in range(len(walker.pos)):
        walker_pos_xarr[i] = walker.pos[i][0]  # converts list of lists to list of arrays
        walker_pos_yarr[i] = walker.pos[i][1]
    H, xedges, yedges = np.histogram2d(walker_pos_xarr, walker_pos_yarr, range=grid_range, bins=bins)
    if H == 'cold':
        H = -H  # make walker visit site negative times
    return H, xedges, yedges


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
        plt.plot(tube_x[i], tube_y[i], c=next(colors), linestyle='--')
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


def plot_histogram_walkers_2d_onlat(timesteps, H_tot, xedges, yedges, quiet, save_dir):
    """Plots temperature profile for all walkers"""
    logging.info("Plotting temperature profile")
    H_tot /= float(timesteps)  # normalization condition
    creation.check_for_folder(save_dir)
    np.savetxt('%s/temp.txt' % save_dir, H_tot, fmt='%.1E')
    plt.title('Temperature profile (dimensionless units)')
    X, Y = np.meshgrid(xedges, yedges)
    temp_profile = H_tot
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


def plot_temp_gradient_2d_onlat(temp_profile, xedges, yedges, quiet, save_dir):
    """Plots temperature gradient for all walkers"""
    logging.info("Plotting temperature gradient")
    temp_gradient_x, temp_gradient_y = np.gradient(temp_profile)
    np.savetxt('%s/temp_gradient.txt' % save_dir, temp_gradient_x, fmt='%.1E')
    plt.title('Temperature gradient $\\frac{dT(x)}{dx}$ (dimensionless units)')
    # Y gradient is irrelevant since we have periodic bd conditions in that direction
    X, Y = np.meshgrid(xedges, yedges)
    plt.pcolormesh(X, Y, temp_gradient_x.T)  # transpose since pcolormesh reverses axes
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
    plt.title(
        "Convergence of $\\frac{dT(x)}{dx}$ calculation sweeping from x to the end\nNoise floor: %.4E" % noise_floor)
    plt.xlabel('X')
    plt.ylabel('Y')
    plt.savefig('%s/temp_gradient_conv.pdf' % save_dir)
    if not quiet:
        plt.show()
    plt.close()
