import ConfigParser
import logging
import os

import creation
import onlat_2d
import onlat_3d

# from mpi4py import MPI


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
    creation.check_for_folder(save_dir)
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
    timesteps = Config.getint('config','timesteps')
    save_dir = Config.get('config', 'save_dir')
    quiet = Config.getboolean('config', 'quiet')
    num_tubes = Config.getint('config', 'num_tubes')
    k_convergence_tolerance = Config.getfloat('config', 'k_convergence_tolerance')
    begin_cov_check = Config.getint('config', 'begin_cov_check')
    k_conv_error_buffer = Config.getint('config', 'k_conv_error_buffer')
    save_loc_plots = Config.getboolean('config', 'save_loc_plots')
    save_loc_data = Config.getboolean('config', 'save_loc_data')
    #mean_dist_tubes = Config.get('config','mean_dist_tubes')
    #std_dist_tubes = Config.get('config', 'std_dist_tubes')

    # Check if inputs valid

    possible_dim = [2, 3]
    if dim not in possible_dim:
        logging.error('Invalid dimension')
        raise SystemExit
    if grid_size < 5:
        logging.error('Invalid grid size')
        raise SystemExit
    if tube_length <= 0:
        logging.error('Invalid tube length')
        raise SystemExit
    if num_tubes < 0:
        logging.error('Invalid number of tubes')
        raise SystemExit
    if timesteps <= 0:
        logging.error('Invalid timesteps')
        raise SystemExit

    os.chdir(save_dir)

    plot_save_dir = creation.get_plot_save_dir(save_dir, num_tubes, orientation, tube_length)
    logging_setup(plot_save_dir)
    logging.info(config_used)

    if (on_lattice == True) & (dim == 2):
        logging.info("Starting 2D on-lattice simulation")
        onlat_2d.sim_2d_onlat(grid_size, tube_length, num_tubes, orientation, timesteps, save_loc_data,
                              quiet, save_loc_plots, save_dir, k_convergence_tolerance, begin_cov_check,
                              k_conv_error_buffer, plot_save_dir)
    elif (on_lattice == True) & (dim == 3):
        tube_diameter = Config.getfloat('config', 'tube_diameter')
        logging.info("Starting 3D on-lattice simulation")
        onlat_3d.sim_3d_onlat(grid_size, tube_length, num_tubes, orientation, timesteps, save_loc_data,
                              quiet, save_loc_plots, save_dir, k_convergence_tolerance, begin_cov_check,
                              k_conv_error_buffer, plot_save_dir, tube_diameter)
    else:
        raise SystemExit
