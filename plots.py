import itertools
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import os
import logging


def check_for_folder(folder):
    if not os.path.exists(folder):
        os.mkdir(folder)


def plot_two_d_random_walk_setup(tube_coords, grid_size):
    """Plots setup and orientation of nanotubes"""
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
    plt.grid()
    check_for_folder('plots')
    plt.savefig('plots/setup.pdf')
    plt.show()


def histogram_walker_2d_onlat(walker, grid_range, temp):
    """Takes walker instance and histograms how many times location is accessed over the simulation. H not normalized
    (for 1 walker only)"""
    #print walker.pos[0][0]
    walker_pos_xarr = np.zeros(len(walker.pos))
    walker_pos_yarr = np.zeros(len(walker.pos))
    for i in range(len(walker.pos)):
        walker_pos_xarr[i] = walker.pos[i][0]  # converts list of lists to list of arrays
        walker_pos_yarr[i] = walker.pos[i][1]
    H, xedges, yedges = np.histogram2d(walker_pos_xarr, walker_pos_yarr, range=grid_range)
    if H == 'cold':
        H = -H  # make walker visit site negative times
    return H, xedges, yedges


def plot_histogram_walker_2d_onlat(walker, timesteps, H, xedges, yedges):  # plots for a single walker
    H_norm = H / float(timesteps)  # normalizes by number of timesteps

    fig = plt.figure

    ax = fig.add_subplot(121)
    ax.set_title('pcolormesh: exact bin edges')
    X, Y = np.meshgrid(xedges, yedges)
    ax.pcolormesh(X, Y, H_norm)
    ax.set_aspect('equal')
    plt.colorbar()

    ax = fig.add_subplot(122)
    ax.set_title('NonUniformImage: interpolated')
    im = mpl.image.NonUniformImage(ax, interpolation='bilinear')
    xcenters = xedges[:-1] + 0.5 * (xedges[1:] - xedges[:-1])
    ycenters = yedges[:-1] + 0.5 * (yedges[1:] - yedges[:-1])
    im.set_data(xcenters, ycenters, H)
    ax.images.append(im)
    ax.set_xlim(xedges[0], xedges[-1])
    ax.set_ylim(yedges[0], yedges[-1])
    ax.set_aspect('equal')
    plt.colorbar()
    plt.show()
    return H, xedges, yedges


def plot_histogram_walkers_2d_onlat(walker, timesteps, num_walkers, H_tot, xedges, yedges):  # plots for all walkers
    H_norm = H_tot / (float(num_walkers)*float(timesteps))

    test = sum(sum(H_norm))
    correction = -1.0/test
    H_norm *= correction
    new_h = sum(sum(H_norm))
    logging.info('Plot is normalized')
    fig = plt.figure()
    #if option == 'equidistant':
     #   plt.title('imshow: equidistant')
    im = plt.imshow(H_norm, interpolation='nearest', origin='low',
                        extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]])
        #plt.colorbar()

    #ax = fig.add_subplot(132)
    plt.title('pcolormesh: exact bin edges')
    X, Y = np.meshgrid(xedges, yedges)
    plt.pcolormesh(X, Y, H_norm.T)  # transpose since pcolormesh reverses axes
    plt.xlim(xedges[0], xedges[-1])
    plt.ylim(yedges[0], yedges[-1])
    plt.colorbar()
    plt.show()

    # ax = fig.add_subplot(133)
    # ax.set_title('NonUniformImage: interpolated')
    # im = mpl.image.NonUniformImage(ax, interpolation='bilinear')
    # xcenters = xedges[:-1] + 0.5 * (xedges[1:] - xedges[:-1])
    # ycenters = yedges[:-1] + 0.5 * (yedges[1:] - yedges[:-1])
    # im.set_data(xcenters, ycenters, H_tot)
    # ax.images.append(im)
    # ax.set_xlim(xedges[0], xedges[-1])
    # ax.set_ylim(yedges[0], yedges[-1])
    # ax.set_aspect('equal')
    # plt.colorbar()
    # plt.show()


def plot_walker_path_all(walkers_pos, grid_size):  # plot all paths for all time steps
    colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])
    print len(walkers_pos)
    print len(walkers_pos[-1])
    pos = np.array((len(walkers_pos),len(walkers_pos[-1])))
    print pos.size
    for i in range(len(walkers_pos[-1])):  # one plot for each walker
        coord = []
        for j in range(len(walkers_pos)):  # runs over timesteps
           coord.append(walkers_pos[j][i].pos)
        # coord now holds list of tuples for one walker at every timestep
        print coord
        plt.plot(coord,c=next(colors))
    plt.xlim(-grid_size, grid_size)
    plt.ylim(-grid_size, grid_size)
    plt.grid()
    plt.savefig('plots/walker_path_all.pdf')
    plt.show()

