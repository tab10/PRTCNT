import plots
import creation
import randomwalk
import ConfigParser
import os
import logging
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
                        filename='%s/run.log' % save_dir,
                        filemode='w')
    console = logging.StreamHandler()
    console.setLevel(logging.INFO)
    formatter = logging.Formatter('%(name)-12s: %(levelname)-8s %(message)s')
    console.setFormatter(formatter)
    logging.getLogger('').addHandler(console)


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
    #mean_dist_tubes = Config.get('config','mean_dist_tubes')
    #std_dist_tubes = Config.get('config', 'std_dist_tubes')

    logging_setup(save_dir)
    logging.info(config_used)

    os.chdir(save_dir)
    logging.info("Setting up grid and tubes")
    grid = creation.Grid2D_onlat(grid_size, tube_length, num_tubes, orientation)
    plots.plot_two_d_random_walk_setup(grid.tube_coords, grid.size)

    # logging.info("Initializing MPI")
    # comm = MPI.COMM_WORLD
    # rank = comm.Get_rank()
    # size = comm.Get_size()
    # machine = platform.node()
    # print("Hello MPI from %s %d of %d" % (machine, rank, size))

    grid_range = [[0, grid.size], [0, grid.size]]

    for i in range(num_walkers):
        # run hot walkers
        logging.info("Start hot walker %d out of %d" % (i+1,num_walkers))
        walker = randomwalk.runrandomwalk_2d_onlat(grid, timesteps, 'hot')
        if i == 0:
            H, xedges, yedges = plots.histogram_walker_2d_onlat(walker, grid_range, 'hot')
            # bin_edges should stay consistent since range is fixed
        else:
            H_temp, xedges, yedges = plots.histogram_walker_2d_onlat(walker, grid_range, 'hot')
            H += H_temp

    for i in range(num_walkers):
        # run cold walkers
        logging.info("Start cold walker %d out of %d" % (i+1,num_walkers))
        walker = randomwalk.runrandomwalk_2d_onlat(grid, timesteps, 'cold')
        if i == 0:
            H, xedges, yedges = plots.histogram_walker_2d_onlat(walker, grid_range, 'cold')
            # bin_edges should stay consistent since range is fixed
        else:
            H_temp, xedges, yedges = plots.histogram_walker_2d_onlat(walker, grid_range, 'cold')
            H -= H_temp

    logging.info("Finished random walks")
    plots.plot_histogram_walkers_2d_onlat(walker, timesteps, num_walkers, H, xedges, yedges)

