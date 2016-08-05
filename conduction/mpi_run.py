import logging
import os
import argparse
import creation
import onlat_2d
import onlat_3d
from mpi4py import MPI


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

    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    parser = argparse.ArgumentParser(description='.')

    parser.add_argument('--dim', type=int, required=True, help='Dimensionality of simulation.')
    parser.add_argument('--grid_size', type=int, default=100, help='Size of square grid of use.')
    parser.add_argument('--tube_length', type=float, required=True, help='Length of nanotubes.')
    parser.add_argument('--num_tubes', type=int, required=True,
                        help='How many tubes are there for random walker to use.')
    parser.add_argument('--tube_radius', type=float, default=0.5,
                        help='Radius of tubes. Only used if kapitza is True.')
    parser.add_argument('--orientation', type=str, required=True, help='Orientation of nanotubes in medium. '
                                                                       'random, horizontal, vertical, or'
                                                                       ' angle in DEGREES.')
    parser.add_argument('--timesteps', type=int, default=10000, help='How many steps to run each walker for.')
    parser.add_argument('--k_convergence_tolerance', type=float, default=1E-05, help='Simulation runs until '
                                                                                 'std. dev. of time fluctuations in k drop below this value.')
    parser.add_argument('--begin_cov_check', type=int, default=100, help='Start checking for convergence '
                                                                         'after this many walkers.')
    parser.add_argument('--k_conv_error_buffer', type=int, default=25, help='Include the last X values of time '
                                                                            'in the std. dev. of k calculation. '
                                                                            'Reduced automatically depending on '
                                                                            '# of cores.')
    parser.add_argument('--gen_plots', type=bool, default=False, help='Gives the option to not generate any plots. '
                                                                      'Useful on the supercomputer.')
    parser.add_argument('--save_dir', type=str, default=os.getcwd(), help='Path for plots and data and config.ini.')
    parser.add_argument('--save_loc_plots', type=bool, default=False, help='Save location plots for all walkers.')
    parser.add_argument('--save_loc_data', type=bool, default=False, help='Save location data for all walkers.')
    parser.add_argument('--quiet', type=bool, default=True, help='Show various plots throughout simulation.')
    parser.add_argument('--on_lattice', type=bool, default=True, help='True for on lattice random walk.')
    parser.add_argument('--kapitza', type=bool, default=False, help='Adds kapitza resistance '
                                                                    'to simulation, see readme.md.')
    args = parser.parse_args()

    comm.Barrier()

    dim = args.dim
    grid_size = args.grid_size
    orientation = args.orientation
    tube_length = args.tube_length
    tube_radius = args.tube_radius
    num_tubes = args.num_tubes
    on_lattice = args.on_lattice
    timesteps = args.timesteps
    save_dir = args.save_dir
    quiet = args.quiet
    num_tubes = args.num_tubes
    k_convergence_tolerance = args.k_convergence_tolerance
    begin_cov_check = args.begin_cov_check
    k_conv_error_buffer = args.k_conv_error_buffer
    save_loc_plots = args.save_loc_plots
    save_loc_data = args.save_loc_data
    gen_plots = args.gen_plots
    kapitza = args.kapitza

    # Check if inputs valid
    possible_dim = [2, 3]
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

    os.chdir(save_dir)

    if rank == 0:
        plot_save_dir = creation.get_plot_save_dir(save_dir, num_tubes, orientation, tube_length)
        logging_setup(plot_save_dir)
    else:
        plot_save_dir = None

    plot_save_dir = comm.bcast(plot_save_dir, root=0)

    if rank == 0:
        logging.info('Using convergence value of %.4E' % k_convergence_tolerance)

    comm.Barrier()

    k_conv_error_buffer /= size  # accounts for extra processors
    if k_conv_error_buffer < 5:  # imposes minimum steps to check convergence at
        k_conv_error_buffer = 5
    if rank == 0:
        logging.info('Using last %d k values to check for convergence' % k_conv_error_buffer)

    #  all processes have control now
    if on_lattice and (dim == 2):
        logging.info("Starting MPI 2D on-lattice simulation")
        comm.Barrier()
        onlat_2d.sim_2d_onlat_MPI(grid_size, tube_length, tube_radius, num_tubes, orientation, timesteps, save_loc_data,
                                  quiet, save_loc_plots, save_dir, k_convergence_tolerance, begin_cov_check,
                                  k_conv_error_buffer, plot_save_dir, gen_plots, kapitza, rank, size)
    elif on_lattice and (dim == 3):
        logging.info("Starting MPI 3D on-lattice simulation")
        comm.Barrier()
        onlat_3d.sim_3d_onlat_MPI(grid_size, tube_length, tube_radius, num_tubes, orientation, timesteps, save_loc_data,
                                  quiet, save_loc_plots, save_dir, k_convergence_tolerance, begin_cov_check,
                                  k_conv_error_buffer, plot_save_dir, tube_radius, gen_plots, kapitza, rank, size)
    else:
        print 'Off lattice not implemented yet'
        raise SystemExit
