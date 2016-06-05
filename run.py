import plots
import creation
import randomwalk
import ConfigParser
import os
import logging
import numpy as np
import analysis

# from mpi4py import MPI
# import platform
# import multiprocessing as mp


def ConfigSectionMap(section):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                print("skip: %s" % option)
        except:
            print("exception on %s!" % option)
            dict1[option] = None
    return dict1


def logging_setup(save_dir):
    logging.basicConfig(level=logging.DEBUG,
                        format='%(asctime)s %(name)-12s %(levelname)-8s %(message)s',
                        datefmt='%m-%d %H:%M',
                        filename='%s/log.txt' % save_dir,
                        filemode='w')
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)


def save_walker_loc(walker, save_dir, walker_index, temp):
    header = "timestep [x,y]\n"
    creation.check_for_folder(save_dir)
    f = open("%s/%s_walker_%d_traj.txt" % (save_dir, temp, walker_index + 1), "w")
    f.write(header)
    for i in range(len(walker.pos)):
        f.write(str(i) + " " + str(walker.pos[i]) + "\n")
    f.close()


if __name__ == "__main__":
    Config = ConfigParser.ConfigParser()

    if os.path.exists("config.ini"):
        Config.read("config.ini")
        config_used = 'Using config.ini'
    else:
        Config.read("default.ini")
        config_used = 'Using default.ini'

    dim = Config.getint('config','dim')
    grid_size = Config.getint('config','grid_size')
    orientation = Config.get('config','orientation')
    tube_length = Config.getfloat('config','tube_length')
    num_tubes = Config.getint('config','num_tubes')
    on_lattice = Config.getboolean('config','on_lattice')
    num_walkers = Config.getint('config','num_walkers')
    timesteps = Config.getint('config','timesteps')
    save_dir = Config.get('config', 'save_dir')
    quiet = Config.getboolean('config', 'quiet')
    save_loc_plots = Config.getboolean('config', 'save_loc_plots')
    save_loc_data = Config.getboolean('config', 'save_loc_data')
    #mean_dist_tubes = Config.get('config','mean_dist_tubes')
    #std_dist_tubes = Config.get('config', 'std_dist_tubes')

    # Check if inputs valid
    if num_walkers <= 0:
        logging.error('Invalid number of walkers')
        raise SystemExit

    logging_setup(save_dir)
    logging.info(config_used)

    os.chdir(save_dir)
    plot_save_dir = save_dir + "/plots"
    walker_data_save_dir = save_dir + "/walker_locations"
    walker_plot_save_dir = save_dir + "/walker_plots"

    logging.info("Setting up grid and tubes")
    grid = creation.Grid2D_onlat(grid_size, tube_length, num_tubes, orientation)
    plots.plot_two_d_random_walk_setup(grid.tube_coords, grid.size, quiet, plot_save_dir)

    # logging.info("Initializing MPI")
    # comm = MPI.COMM_WORLD
    # rank = comm.Get_rank()
    # size = comm.Get_size()
    # machine = platform.node()
    # print("Hello MPI from %s %d of %d" % (machine, rank, size))

    grid_range = [[0, grid.size], [0, grid.size]]
    bins = grid.size
    H = np.zeros((grid.size, grid.size))

    for i in range(num_walkers):
        # run hot walker
        logging.info("Start hot walker %d out of %d" % (i+1,num_walkers))
        walker = randomwalk.runrandomwalk_2d_onlat(grid, timesteps, 'hot')
        if save_loc_data:
            save_walker_loc(walker, walker_data_save_dir, i, 'hot')
        if i == 0 & save_loc_plots == False:  # always save one example trajectory plot
            plots.plot_walker_path_2d_onlat(walker, grid_size, 'hot', quiet, i + 1, plot_save_dir)
        elif save_loc_plots:
            plots.plot_walker_path_2d_onlat(walker, grid_size, 'hot', quiet, i + 1, walker_plot_save_dir)
        H_temp, xedges, yedges = plots.histogram_walker_2d_onlat(walker, grid_range, 'hot', bins)
        H += H_temp

        # run cold walker
        logging.info("Start cold walker %d out of %d" % (i+1,num_walkers))
        walker = randomwalk.runrandomwalk_2d_onlat(grid, timesteps, 'cold')
        if save_loc_data:
            save_walker_loc(walker, walker_data_save_dir, i, 'cold')
        if i == 0 & save_loc_plots == False:
            plots.plot_walker_path_2d_onlat(walker, grid_size, 'cold', quiet, i + 1, plot_save_dir)
        elif save_loc_plots:
            plots.plot_walker_path_2d_onlat(walker, grid_size, 'cold', quiet, i + 1, walker_plot_save_dir)
        H_temp, xedges, yedges = plots.histogram_walker_2d_onlat(walker, grid_range, 'cold', bins)
        H -= H_temp


    logging.info("Finished random walks")
    temp_profile = plots.plot_histogram_walkers_2d_onlat(timesteps, H, xedges, yedges, quiet, plot_save_dir)
    temp_gradient_x = plots.plot_temp_gradient_2d_onlat(temp_profile, xedges, yedges, quiet, plot_save_dir)
    plots.plot_check_gradient_noise_floor(temp_gradient_x, quiet, plot_save_dir)
    analysis.calc_conductivity_2d_onlat(num_walkers, temp_gradient_x, grid.size, timesteps)