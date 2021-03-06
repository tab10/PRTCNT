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


# from future import standard_library
import logging
import sys
import importlib
import os

from conduction import backend
from conduction import test_2d
from conduction import test_3d
from conduction import randomwalk_2d
from conduction import randomwalk_3d


#standard_library.install_aliases()

def ConfigSectionMap(section):
    dict1 = {}
    options = Config.options(section)
    for option in options:
        try:
            dict1[option] = Config.get(section, option)
            if dict1[option] == -1:
                print(("skip: %s" % option))
        except:
            print(("exception on %s!" % option))
            dict1[option] = None
    return dict1


def logging_setup(save_dir):
    backend.check_for_folder(save_dir)
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
    backend.check_for_folder(save_dir)
    f = open("%s/%s_walker_%d_traj.txt" % (save_dir, temp, walker_index + 1), "w")
    f.write(header)
    for i in range(len(walker.pos)):
        f.write(str(i) + " " + str(walker.pos[i]) + "\n")
    f.close()


if __name__ == "__main__":

    print("Serial code is currently buggy! Please use parallel code with mpirun -np 2 on at least dual cores.")
    raise SystemExit

    print(sys.version_info[0])
    if sys.version_info[0] < 3:
        importlib.import_module('ConfigParser')
        Config = ConfigParser.ConfigParser()
    else:
        print('ji')
        importlib.import_module('configparser')
        Config = configparser.configparser()

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
    tube_radius = Config.getfloat('config', 'tube_radius')
    num_tubes = Config.getint('config','num_tubes')
    timesteps = Config.getint('config','timesteps')
    save_dir = Config.get('config', 'save_dir')
    quiet = Config.getboolean('config', 'quiet')
    num_tubes = Config.getint('config', 'num_tubes')
    k_convergence_tolerance = Config.getfloat('config', 'k_convergence_tolerance')
    begin_cov_check = Config.getint('config', 'begin_cov_check')
    k_conv_error_buffer = Config.getint('config', 'k_conv_error_buffer')
    save_loc_plots = Config.getboolean('config', 'save_loc_plots')
    save_loc_data = Config.getboolean('config', 'save_loc_data')
    gen_plots = Config.getboolean('config', 'gen_plots')
    kapitza = Config.getboolean('config', 'kapitza')
    prob_m_cn = Config.getfloat('config', 'prob_m_cn')
    run_to_convergence = Config.getboolean('config', 'run_to_convergence')
    num_walkers = Config.getint('config', 'num_walkers')
    disable_func = Config.getboolean('config', 'disable_func')
    rules_test = Config.getboolean('rules_simulation', 'rules_test')
    printout_inc = Config.getint('constant_flux', 'printout_inc')

    #mean_dist_tubes = Config.get('config','mean_dist_tubes')
    #std_dist_tubes = Config.get('config', 'std_dist_tubes')
    # Check if inputs valid

    os.chdir(save_dir)

    plot_save_dir = backend.get_plot_save_dir(save_dir, num_tubes, orientation, tube_length)
    logging_setup(plot_save_dir)
    logging.info(config_used)

    ##### VALUE & COMMON SENSE CHECKS#####
    possible_dim = [2, 3]
    if dim not in possible_dim:
        logging.error('Invalid dimension')
        raise SystemExit
    if kapitza and tube_radius == 0:
        logging.error('Kapitza modeling requires tubes to have a nonzero radius')
        raise SystemExit
    if grid_size < 5:
        logging.error('Invalid grid size')
        raise SystemExit
    if tube_length < 0:
        logging.error('Invalid tube length')
        raise SystemExit
    if num_tubes < 0:
        logging.error('Invalid number of tubes')
        raise SystemExit
    if timesteps <= 0:
        logging.error('Invalid timesteps')
        raise SystemExit
    if kapitza:
        logging.info('Kapitza modeling is ON')
        logging.info('Using prob_m_cn value of %.4f' % prob_m_cn)
    else:
        logging.info('Kapitza modeling is OFF')
        prob_m_cn = 0.0

    logging.info('Simulation parameters:\n%d, P_m-cn=%.2f, %s model, %d total walkers, %d timesteps'
                 % (dim, prob_m_cn, model, num_walkers, timesteps))

    if disable_func:
        logging.info('Functionalization has been disabled, treating ends as volume in rules')

    logging.info('Grid size of %d is being used' % (grid_size + 1))
    ##### #####

    if rules_test:
        logging.info("Starting rules test only random walk. Rules will be checked to ensure they uphold the"
                     "Principle of Detailed Balance.")
        logging.info("k will not be computed, so only simulation setup parameters are used. All others"
                     "are ignored.")
        if dim == 2:
            test_2d.serial_method(grid_size, tube_length, tube_radius, num_tubes, orientation,
                                                 timesteps, quiet, plot_save_dir, gen_plots, kapitza, prob_m_cn,
                                                 num_walkers, printout_inc, k_conv_error_buffer, disable_func)
        elif dim == 3:
            test_3d.serial_method(grid_size, tube_length, tube_radius, num_tubes, orientation, timesteps,
                                  quiet, plot_save_dir,
                                  gen_plots, kapitza, prob_m_cn, num_walkers, printout_inc,
                                  k_conv_error_buffer, disable_func)

    elif not rules_test:
        logging.info("Starting %dD constant flux on-lattice random walk." % dim)
        if dim == 2:
            randomwalk_2d.serial_method(grid_size, tube_length, tube_radius, num_tubes, orientation,
                                        timesteps, quiet, plot_save_dir, gen_plots, kapitza, prob_m_cn,
                                        num_walkers, printout_inc, k_conv_error_buffer, disable_func)
        elif dim == 3:
            randomwalk_3d.serial_method(grid_size, tube_length, tube_radius, num_tubes, orientation, timesteps,
                                        quiet, plot_save_dir,
                                        gen_plots, kapitza, prob_m_cn, num_walkers, printout_inc,
                                        k_conv_error_buffer, disable_func)
    else:
        logging.error('Check inputs')
        raise SystemExit
