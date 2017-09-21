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


# from __future__ import absolute_import
import logging
import os
import argparse
from mpi4py import MPI
import ast

from conduction import backend
from conduction import test_3d
from conduction import test_2d
from conduction import randomwalk_3d
from conduction import randomwalk_2d

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

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    parser = argparse.ArgumentParser(description='.')

    parser.add_argument('--dim', type=int, required=True, help='Dimensionality of simulation.')
    parser.add_argument('--grid_size', type=int, default=99, help='Size of grid of use. TRUE SIZE USED IS VALUE + 1, '
                                                                  'TO COMPARE WITH ANALYTICAL.')
    parser.add_argument('--tube_length', type=float, default=15, help='Length of nanotubes.')
    parser.add_argument('--num_tubes', type=int, default=0,
                        help='How many tubes are there for random walker to use.')
    parser.add_argument('--orientation', type=str, default='horizontal', help='Orientation of nanotubes in medium. '
                                                                       'random, horizontal, vertical, or'
                                                                       ' angle in DEGREES.')
    parser.add_argument('--timesteps', type=int, default=25000, help='How many steps to run each walker for. '
                                                                     'Should be (grid_size+1)**2 to have even '
                                                                     'temperature distribution.')
    parser.add_argument('--k_convergence_tolerance', type=float, default=1E-05, help='Simulation runs until '
                                                                                 'std. dev. of time fluctuations in k drop below this value.')
    parser.add_argument('--begin_cov_check', type=int, default=100, help='Start checking for convergence '
                                                                         'after this many walkers.')
    parser.add_argument('--k_conv_error_buffer', type=int, default=25, help='Include the last X values of time '
                                                                            'in the std. dev. of k calculation. '
                                                                            'Reduced automatically depending on '
                                                                            '# of cores.')
    parser.add_argument('--gen_plots', type=str, default='True', help='Gives the option to not generate any plots. '
                                                                      'Useful on the supercomputer.')
    parser.add_argument('--save_dir', type=str, default=os.getcwd(), help='Path for plots and data and config.ini.')
    parser.add_argument('--save_loc_plots', type=str, default='False', help='Save location plots for all walkers.')
    parser.add_argument('--quiet', type=str, default='True', help='Do not show various plots throughout simulation.')
    parser.add_argument('--model', type=str, required=True, help='Simulation model type. kapitza, tunneling_w_vol, '
                                                                 'tunneling_wo_vol')
    parser.add_argument('--prob_m_cn', type=float, default=0.5, help='Probability a walker will enter the CNT. '
                                                                     'Only used in kapitza models.')
    parser.add_argument('--run_to_convergence', type=str, default='False', help='True does this or False runs '
                                                                              'for number of walkers.')
    parser.add_argument('--restart', type=str, default='False', help='Looks in previous directory for H to extend or '
                                                                    'restart simulation.')
    parser.add_argument('--num_walkers', type=int, default=50000, help='Total walkers to use for simulaton. '
                                                                      'Only used if convergence is false.')
    parser.add_argument('--disable_func', type=str, default='False',
                        help='Turn off functionalization of the tube ends.')
    parser.add_argument('--printout_inc', type=int, default=50, help='deltaT increment for printing out conductivity '
                                                                     'info for constant flux simulations. Should be '
                                                                     'somewhat large because histogramming has to be done every time.')
    parser.add_argument('--rules_test', type=str, default=False, help='Starts a rules test only simulation. '
                                                                      'This checks that the simulation will obey'
                                                                      'detailed balance. Available with serial'
                                                                      'or MPI options as in the primary '
                                                                      'program.')

    args = parser.parse_args()

    comm.Barrier()

    for k in args.__dict__:  # convert string "BOOLEANS" to actual booleans
        tester = args.__dict__[k]
        if (tester == 'True' or tester == 'False') and isinstance(tester, str):
            new_test = ast.literal_eval(tester)
            args.__dict__[k] = new_test

    dim = args.dim
    grid_size = args.grid_size
    orientation = args.orientation
    tube_length = args.tube_length
    num_tubes = args.num_tubes
    timesteps = args.timesteps
    save_dir = args.save_dir
    quiet = args.quiet
    num_tubes = args.num_tubes
    k_convergence_tolerance = args.k_convergence_tolerance
    begin_cov_check = args.begin_cov_check
    k_conv_error_buffer = args.k_conv_error_buffer
    save_loc_plots = args.save_loc_plots
    gen_plots = args.gen_plots
    prob_m_cn = args.prob_m_cn
    run_to_convergence = args.run_to_convergence
    num_walkers = args.num_walkers
    printout_inc = args.printout_inc
    restart = args.restart
    disable_func = args.disable_func
    rules_test = args.rules_test
    model = args.model

    os.chdir(save_dir)

    if rank == 0:
        plot_save_dir = backend.get_plot_save_dir(save_dir, num_tubes, orientation, tube_length, restart)
        logging_setup(plot_save_dir)
    else:
        plot_save_dir = None

    plot_save_dir = comm.bcast(plot_save_dir, root=0)

    ##### VALUE & COMMON SENSE CHECKS#####
    possible_dim = [2, 3]
    if model == 'kapitza':
        tube_radius = 0.5
        kapitza = True
        inert_vol = False
        logging.info('Simulation model: kapitza. CNTs have volume, functionalized/non-functionalized ends, and '
                     'kapitza resistance')
        logging.info('Using prob_m_cn value of %.4f' % prob_m_cn)
    elif model == 'tunneling_w_vol':
        tube_radius = 0.5
        prob_m_cn = 0.0
        kapitza = False
        inert_vol = True
        disable_func = False
        logging.info('Simulation model: tunneling with volume. CNTs have (non-functioning/excluded) volume and '
                     'functionalized ends. CNTs CANNOT cross in space. Limit of infinite kapitza resistance.')
    elif model == 'tunneling_wo_vol':
        tube_radius = 0.0
        prob_m_cn = 0.0
        kapitza = False
        inert_vol = False
        disable_func = False
        logging.info('Simulation model: tunneling without volume. CNTs only have '
                     'functionalized ends. CNTs CANNOT cross in space. Models the tunneling case when'
                     'there is no excluded volume effect. Limit of infinite kapitza resistance.')
    else:
        logging.error('Incorrect simulation model specified.')
        raise SystemExit
    if disable_func:
        logging.info('Functionalization of ends DISABLED')
    else:
        logging.info('Functionalization of ends ENABLED')
    if dim not in possible_dim:
        logging.error('Invalid dimension')
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
    if run_to_convergence:
        logging.info('Simulation will run to convergence')
    else:
        logging.info('Simulation will run to %d walkers' % num_walkers)
        if begin_cov_check >= num_walkers:
            logging.warning('begin_cov_check is less than or equal to num_walkers, forcing 3*num_walkers')
            num_walkers *= 3
    logging.info('Grid size of %d is being used' % (grid_size + 1))
    if disable_func:
        logging.info('Functionalization has been disabled, treating ends as volume in rules.')
    ##### #####

    comm.Barrier()

    k_conv_error_buffer /= size  # accounts for extra processors
    if k_conv_error_buffer < 5:  # imposes minimum steps to check convergence at
        k_conv_error_buffer = 5
    if rank == 0:
        logging.info('Using last %d k values to check for convergence' % k_conv_error_buffer)

    #  all processes have control now
    if rules_test:
        logging.info("Starting rules test only random walk.\nRules will be checked to ensure they uphold the"
                     "Principle of Detailed Balance.\nTo do this, we don't use Fourier's Law as our walkers are"
                     " all positive and no heat flux is generated.\nDifferences include: Walkers can start from "
                     "anywhere in the box,\nALL boundaries are periodic, Walkers are all positive,\nALL visited"
                     " positions are histogrammed as opposed to keeping just 1")
        if dim == 2:
            test_2d.parallel_method(grid_size, tube_length, tube_radius, num_tubes, orientation,
                                    timesteps, quiet, plot_save_dir, gen_plots, kapitza, prob_m_cn,
                                    num_walkers, disable_func, rank, size, rules_test, restart, inert_vol)
        elif dim == 3:
            test_3d.parallel_method(grid_size, tube_length, tube_radius, num_tubes, orientation, timesteps, quiet,
                                    plot_save_dir, gen_plots, kapitza, prob_m_cn, num_walkers, disable_func, rank,
                                    size, rules_test, restart, inert_vol)
    else:
        logging.info("Starting %dD constant flux on-lattice random walk." % dim)
        if dim == 2:
            randomwalk_2d.parallel_method(grid_size, tube_length, tube_radius, num_tubes, orientation, timesteps,
                                          quiet, plot_save_dir, gen_plots, kapitza, prob_m_cn,
                                          num_walkers, printout_inc, k_conv_error_buffer, disable_func, rank, size,
                                          rules_test, restart, inert_vol, save_loc_plots)
        elif dim == 3:
            randomwalk_3d.parallel_method(grid_size, tube_length, tube_radius, num_tubes, orientation, timesteps,
                                          quiet, plot_save_dir, gen_plots, kapitza, prob_m_cn, num_walkers,
                                          printout_inc,
                                          k_conv_error_buffer, disable_func, rank, size, rules_test, restart, inert_vol,
                                          save_loc_plots)
