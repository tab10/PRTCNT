import logging
import numpy as np
import analysis
import creation
import plots
import randomwalk
import run
from mpi4py import MPI


def sim_3d_onlat(grid_size, tube_length, num_tubes, orientation, timesteps, save_loc_data,
                 quiet, save_loc_plots, save_dir, k_convergence_tolerance, begin_cov_check,
                 k_conv_error_buffer, plot_save_dir, tube_radius, gen_plots):
    walker_data_save_dir = save_dir + "/walker_locations"
    walker_plot_save_dir = save_dir + "/walker_plots"

    grid = creation.Grid3D_onlat(grid_size, tube_length, tube_radius, num_tubes, orientation)
    if gen_plots:
        plots.plot_three_d_random_walk_setup(grid.tube_coords, grid.size, quiet, plot_save_dir)
    fill_fract = tube_length * float(num_tubes) / grid.size ** 3
    logging.info("Filling fraction is %.2f %%" % (fill_fract * 100.0))

    grid_range = [[0, grid.size], [0, grid.size], [0, grid.size]]
    bins = grid.size
    H = np.zeros((grid.size, grid.size, grid.size))

    i = 0
    k_list = []
    k_convergence_err_list = []
    k_convergence_err = 1.0

    while k_convergence_err > k_convergence_tolerance:
        # for i in range(num_walkers):

        # run hot walker
        # logging.info("Start hot walker %d" % (i+1))
        walker = randomwalk.runrandomwalk_3d_onlat(grid, timesteps, 'hot')
        if save_loc_data:
            run.save_walker_loc(walker, walker_data_save_dir, i, 'hot')
        if i == 0 & save_loc_plots == False & gen_plots == True:  # always save one example trajectory plot
            plots.plot_walker_path_3d_onlat(walker, grid_size, 'hot', quiet, i + 1, plot_save_dir)
        elif save_loc_plots & gen_plots == True:
            plots.plot_walker_path_3d_onlat(walker, grid_size, 'hot', quiet, i + 1, walker_plot_save_dir)
        H_temp, x_edges, y_edges, z_edges = plots.histogram_walker_3d_onlat(walker, grid_range, bins)
        H += H_temp

        # run cold walker
        # logging.info("Start cold walker %d" % (i+1))
        walker = randomwalk.runrandomwalk_3d_onlat(grid, timesteps, 'cold')
        if save_loc_data:
            run.save_walker_loc(walker, walker_data_save_dir, i, 'cold')
        if i == 0 & save_loc_plots == False & gen_plots == True:
            plots.plot_walker_path_3d_onlat(walker, grid_size, 'cold', quiet, i + 1, plot_save_dir)
        elif save_loc_plots & gen_plots == True:
            plots.plot_walker_path_3d_onlat(walker, grid_size, 'cold', quiet, i + 1, walker_plot_save_dir)
        H_temp, x_edges, y_edges, z_edges = plots.histogram_walker_3d_onlat(walker, grid_range, bins)
        H -= H_temp

        i += 1

        dt_dx, heat_flux, dt_dx_err, k, k_err, r2, temp_profile_sum = analysis.check_convergence_3d_onlat(H, i * 2,
                                                                                                          grid.size,
                                                                                                          timesteps)
        k_list.append(k)
        logging.info("%d: R squared: %.4f, k: %.4E" % (i, r2, k))
        if i > begin_cov_check:
            k_convergence_err = np.std(np.array(k_list[-k_conv_error_buffer:]), ddof=1)
            k_convergence_val = np.mean(np.array(k_list[-k_conv_error_buffer:]))
            k_convergence_err_list.append(k_convergence_err)
            logging.info("k: %.4E" % k_convergence_val)
            logging.info("k error: %.4E" % k_convergence_err)

    logging.info("Simulation has converged with %d total walkers" % (i * 2))
    logging.info("Finished random walks")
    if gen_plots:
        plots.plot_k_convergence(k_list, quiet, plot_save_dir)
        plots.plot_k_convergence_err(k_convergence_err_list, quiet, plot_save_dir, begin_cov_check)
    gradient_avg, gradient_std = plots.plot_linear_temp(temp_profile_sum, quiet, plot_save_dir, gen_plots)
    analysis.final_conductivity_3d_onlat(i * 2, grid.size, timesteps, gradient_avg,
                                         gradient_std, k_convergence_err, num_tubes, save_dir,
                                         k_convergence_val, gradient_cutoff=2)
    logging.info("Complete")


def sim_3d_onlat_MPI(grid_size, tube_length, num_tubes, orientation, timesteps, save_loc_data,
                     quiet, save_loc_plots, save_dir, k_convergence_tolerance, begin_cov_check,
                     k_conv_error_buffer, plot_save_dir, tube_radius, gen_plots, rank, size):
    comm = MPI.COMM_WORLD
    walker_data_save_dir = save_dir + "/walker_locations"
    walker_plot_save_dir = save_dir + "/walker_plots"

    if rank == 0:
        grid = creation.Grid3D_onlat(grid_size, tube_length, tube_radius, num_tubes, orientation)
        if gen_plots:
            plots.plot_three_d_random_walk_setup(grid.tube_coords, grid.size, quiet, plot_save_dir)
        fill_fract = tube_length * float(num_tubes) / grid.size ** 3
        logging.info("Filling fraction is %.2f %%" % (fill_fract * 100.0))
    else:
        grid = None

    grid = comm.bcast(grid, root=0)

    grid_range = [[0, grid.size], [0, grid.size], [0, grid.size]]
    bins = grid.size
    H = np.zeros((grid.size, grid.size, grid.size))
    tot_H = np.zeros((grid.size, grid.size, grid.size))

    i = 0

    k_list = []
    k_convergence_err_list = []
    k_convergence_err = 1.0

    start = MPI.Wtime()

    while k_convergence_err > k_convergence_tolerance:
        # for i in range(num_walkers):

        # run hot walker
        # logging.info("Start hot walker %d" % (i+1))
        walker = randomwalk.runrandomwalk_3d_onlat(grid, timesteps, 'hot')

        if rank == 0:
            if save_loc_data:
                run.save_walker_loc(walker, walker_data_save_dir, i, 'hot')
            if i == 0 & save_loc_plots == False & gen_plots == True:  # always save one example trajectory plot
                plots.plot_walker_path_3d_onlat(walker, grid_size, 'hot', quiet, i + 1, plot_save_dir)
            elif save_loc_plots & gen_plots == True:
                plots.plot_walker_path_3d_onlat(walker, grid_size, 'hot', quiet, i + 1, walker_plot_save_dir)

        H_temp, x_edges, y_edges, z_edges = plots.histogram_walker_3d_onlat(walker, grid_range, bins)
        H += H_temp

        # run cold walker
        # logging.info("Start cold walker %d" % (i+1))
        walker = randomwalk.runrandomwalk_3d_onlat(grid, timesteps, 'cold')

        if rank == 0:
            if save_loc_data:
                run.save_walker_loc(walker, walker_data_save_dir, i, 'cold')
            if i == 0 & save_loc_plots == False & gen_plots == True:
                plots.plot_walker_path_3d_onlat(walker, grid_size, 'cold', quiet, i + 1, plot_save_dir)
            elif save_loc_plots & gen_plots == True:
                plots.plot_walker_path_3d_onlat(walker, grid_size, 'cold', quiet, i + 1, walker_plot_save_dir)

        H_temp, x_edges, y_edges, z_edges = plots.histogram_walker_3d_onlat(walker, grid_range, bins)
        H -= H_temp

        comm.Reduce(H, tot_H, op=MPI.SUM, root=0)
        # H is updated on every core for every i independently
        # tot_H is the total across all cores

        comm.Barrier()
        i += 1  # i starts at 0

        if rank == 0:
            dt_dx, heat_flux, dt_dx_err, k, k_err, r2, temp_profile_sum = analysis.check_convergence_3d_onlat(tot_H,
                                                                                                              i * 2 * size,
                                                                                                              grid.size,
                                                                                                              timesteps)
            k_list.append(k)
            logging.info("%d: R squared: %.4f, k: %.4E" % (i * size, r2, k))

        comm.Barrier()

        if (i * size) > begin_cov_check:
            if rank == 0:
                k_convergence_err = np.std(np.array(k_list[-k_conv_error_buffer:]), ddof=1)
                k_convergence_val = np.mean(np.array(k_list[-k_conv_error_buffer:]))
                k_convergence_err_list.append(k_convergence_err)
                logging.info("k: %.4E" % k_convergence_val)
                logging.info("k error: %.4E" % k_convergence_err)
            else:
                k_convergence_err = None
                k_convergence_val = None

            k_convergence_err = comm.bcast(k_convergence_err, root=0)
            k_convergence_val = comm.bcast(k_convergence_val, root=0)

        comm.Barrier()

    if rank == 0:
        end = MPI.Wtime()
        logging.info("Simulation has converged with %d total walkers" % (i * 2 * size))
        logging.info("Finished random walks")
        logging.info("Using %d cores, parallel simulation time was %.4f s" % (size, end - start))
        walk_sec = (i * 2 * size) / (end - start)
        logging.info("Crunched %.4f walkers/second" % walk_sec)
        if gen_plots:
            plots.plot_k_convergence(k_list, quiet, plot_save_dir)
            plots.plot_k_convergence_err(k_convergence_err_list, quiet, plot_save_dir, begin_cov_check)
        gradient_avg, gradient_std = plots.plot_linear_temp(temp_profile_sum, quiet, plot_save_dir, gen_plots)
        analysis.final_conductivity_3d_onlat(i * 2 * size, grid.size, timesteps, gradient_avg, gradient_std,
                                             k_convergence_err, num_tubes, plot_save_dir, k_convergence_val,
                                             gradient_cutoff=2)
        logging.info("Complete")
