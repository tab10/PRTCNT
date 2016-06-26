import logging

import numpy as np
import time
import analysis
import creation
import plots
import randomwalk
import run
from mpi4py import MPI


def sim_2d_onlat(grid_size, tube_length, num_tubes, orientation, timesteps, save_loc_data,
                 quiet, save_loc_plots, save_dir, k_convergence_tolerance, begin_cov_check,
                 k_conv_error_buffer, plot_save_dir, gen_plots):
    walker_data_save_dir = save_dir + "/walker_locations"
    walker_plot_save_dir = save_dir + "/walker_plots"

    logging.info("Setting up grid and tubes")
    grid = creation.Grid2D_onlat(grid_size, tube_length, num_tubes, orientation)
    if gen_plots:
        plots.plot_two_d_random_walk_setup(grid.tube_coords, grid.size, quiet, plot_save_dir)
    fill_fract = analysis.filling_fraction(grid.tube_coords, grid.size)
    logging.info("Filling fraction is %.2f" % fill_fract)

    grid_range = [[0, grid.size], [0, grid.size]]
    bins = grid.size
    H = np.zeros((grid.size, grid.size))

    i = 0
    k_list = []
    k_convergence_err_list = []
    k_convergence_err = 1.0

    start = time.clock()

    while k_convergence_err > k_convergence_tolerance:
        # for i in range(num_walkers):

        # run hot walker
        # logging.info("Start hot walker %d" % (i+1))
        walker = randomwalk.runrandomwalk_2d_onlat(grid, timesteps, 'hot')
        if save_loc_data:
            run.save_walker_loc(walker, walker_data_save_dir, i, 'hot')
        if i == 0 & save_loc_plots == False & gen_plots == True:  # always save one example trajectory plot
            plots.plot_walker_path_2d_onlat(walker, grid_size, 'hot', quiet, i + 1, plot_save_dir)
        elif save_loc_plots & gen_plots == True:
            plots.plot_walker_path_2d_onlat(walker, grid_size, 'hot', quiet, i + 1, walker_plot_save_dir)
        H_temp, xedges, yedges = plots.histogram_walker_2d_onlat(walker, grid_range, bins)
        H += H_temp

        # run cold walker
        # logging.info("Start cold walker %d" % (i+1))
        walker = randomwalk.runrandomwalk_2d_onlat(grid, timesteps, 'cold')
        if save_loc_data:
            run.save_walker_loc(walker, walker_data_save_dir, i, 'cold')
        if i == 0 & save_loc_plots == False & gen_plots == True:
            plots.plot_walker_path_2d_onlat(walker, grid_size, 'cold', quiet, i + 1, plot_save_dir)
        elif save_loc_plots & gen_plots == True:
            plots.plot_walker_path_2d_onlat(walker, grid_size, 'cold', quiet, i + 1, walker_plot_save_dir)
        H_temp, xedges, yedges = plots.histogram_walker_2d_onlat(walker, grid_range, bins)
        H -= H_temp

        i += 1

        dt_dx, heat_flux, dt_dx_err, k, k_err, r2 = analysis.check_convergence_2d_onlat(H, i * 2, grid.size, timesteps)
        k_list.append(k)
        logging.info("%d: R squared: %.4f, k: %.4E" % (i * 2, r2, k))
        if i > begin_cov_check:
            k_convergence_err = np.std(np.array(k_list[-k_conv_error_buffer:]), ddof=1)
            k_convergence_val = np.mean(np.array(k_list[-k_conv_error_buffer:]))
            k_convergence_err_list.append(k_convergence_err)
            logging.info("k: %.4E" % k_convergence_val)
            logging.info("k error: %.4E" % k_convergence_err)

    end = time.clock()
    logging.info("Simulation has converged with %d total walkers" % (i * 2))
    logging.info("Finished random walks")
    logging.info("Serial simulation time was %.4f s" % (end - start))
    walk_sec = (i * 2) / (end - start)
    logging.info("Crunched %.4f walkers/second" % walk_sec)
    temp_profile = plots.plot_histogram_walkers_2d_onlat(timesteps, H, xedges, yedges, quiet, plot_save_dir, gen_plots)
    if gen_plots == True:
        plots.plot_k_convergence(k_list, quiet, plot_save_dir)
        plots.plot_k_convergence_err(k_convergence_err_list, quiet, plot_save_dir, begin_cov_check)
        temp_gradient_x = plots.plot_temp_gradient_2d_onlat(temp_profile, xedges, yedges, quiet,
                                                        plot_save_dir, gradient_cutoff=2)
    gradient_avg, gradient_std = plots.plot_linear_temp(temp_profile, quiet, plot_save_dir, gen_plots)
    analysis.final_conductivity_2d_onlat(i * 2, grid.size, timesteps, gradient_avg, gradient_std,
                                         k_convergence_err, num_tubes, plot_save_dir, k_convergence_val,
                                         gradient_cutoff=2)
    logging.info("Complete")


def sim_2d_onlat_MPI(grid_size, tube_length, num_tubes, orientation, timesteps, save_loc_data,
                     quiet, save_loc_plots, save_dir, k_convergence_tolerance, begin_cov_check,
                     k_conv_error_buffer, plot_save_dir, gen_plots, rank, size):

    comm = MPI.COMM_WORLD
    walker_data_save_dir = save_dir + "/walker_locations"
    walker_plot_save_dir = save_dir + "/walker_plots"
    if size < 12:
        logging.info("Less than 12 cores detected. Using serial convergence algorithm")
    elif size >= 12:
        logging.info("12 or more cores detected. Using parallel convergence algorithm")
    if rank == 0:
        logging.info("Setting up grid and tubes")
        grid = creation.Grid2D_onlat(grid_size, tube_length, num_tubes, orientation)
        if gen_plots:
            plots.plot_two_d_random_walk_setup(grid.tube_coords, grid.size, quiet, plot_save_dir)
        fill_fract = analysis.filling_fraction(grid.tube_coords, grid.size)
        logging.info("Filling fraction is %.2f" % fill_fract)
    else:
        grid = None

    grid = comm.bcast(grid, root=0)

    grid_range = [[0, grid.size], [0, grid.size]]
    bins = grid.size
    H = np.zeros((grid.size, grid.size))

    i = 0

    k_list = []
    k_convergence_err_list = []
    k_convergence_err = 1.0
    k_core = np.zeros(size)
    r2_core = np.zeros(size)

    start = MPI.Wtime()

    while k_convergence_err > k_convergence_tolerance:
        # for i in range(num_walkers):

        # run hot walker
        # logging.info("Start hot walker %d" % (i+1))
        walker = randomwalk.runrandomwalk_2d_onlat(grid, timesteps, 'hot')

        if rank == 0:
            if save_loc_data:
                run.save_walker_loc(walker, walker_data_save_dir, i, 'hot')
            if i == 0 & save_loc_plots == False & gen_plots == True:  # always save one example trajectory plot
                plots.plot_walker_path_2d_onlat(walker, grid_size, 'hot', quiet, i + 1, plot_save_dir)
            elif save_loc_plots & gen_plots == True:
                plots.plot_walker_path_2d_onlat(walker, grid_size, 'hot', quiet, i + 1, walker_plot_save_dir)

        H_temp, xedges, yedges = plots.histogram_walker_2d_onlat(walker, grid_range, bins)
        H += H_temp

        # run cold walker
        # logging.info("Start cold walker %d" % (i+1))
        walker = randomwalk.runrandomwalk_2d_onlat(grid, timesteps, 'cold')

        if rank == 0:
            if save_loc_data:
                run.save_walker_loc(walker, walker_data_save_dir, i, 'cold')
            if i == 0 & save_loc_plots == False & gen_plots == True:
                plots.plot_walker_path_2d_onlat(walker, grid_size, 'cold', quiet, i + 1, plot_save_dir)
            elif save_loc_plots & gen_plots == True:
                plots.plot_walker_path_2d_onlat(walker, grid_size, 'cold', quiet, i + 1, walker_plot_save_dir)

        H_temp, xedges, yedges = plots.histogram_walker_2d_onlat(walker, grid_range, bins)
        H -= H_temp

        tot_H = comm.allreduce(H, op=MPI.SUM)
        # H is updated on every core for every i independently
        # tot_H is the total across all cores
        i += 1

        if size >= 12:
            dt_dx, heat_flux, dt_dx_err, k, k_err, r2 = analysis.check_convergence_2d_onlat(H, i * 2 * (rank + 1),
                                                                                            grid.size, timesteps)
            logging.info("%d: R squared: %.4f, k: %.4E" % (i * 2 * (rank + 1), r2, k))
            if rank != 0:
                comm.send(k, dest=0, tag=10)
            else:
                k_list.append(k)
                for j in range(1, size):
                    temp_k = comm.recv(source=j, tag=10)
                    k_list.append(temp_k)
            print size(k_list)
                    #
                    #
                    # if rank != 0:
                    #     comm.send(k, dest=0, tag=10)
                    #     comm.send(r2, dest=0, tag=15)
                    # else:
                    #     for j in range(1, size):
                    #         k_core[j] = comm.recv(source=j, tag=10)
                    #         r2_core[j] = comm.recv(source=j, tag=15)
                    #     k_core[0] = k
                    #     r2_core[0] = r2
                    #     print k_core
                    #     print r2_core
                    #     k_error = np.std(k_core, ddof=1)
                    #     k_mean = np.mean(k_core)
                    #     r2_error = np.std(r2_core, ddof=1)
                    #     r2_mean = np.mean(r2_core)
                    #     logging.info("%d: Avg. across cores R squared: %.4f, k: %.4E" % (i*2*size, r2_mean, k_mean))


        elif size < 12:
            if rank == 0:
                dt_dx, heat_flux, dt_dx_err, k, k_err, r2 = analysis.check_convergence_2d_onlat(tot_H, i * 2 * size,
                                                                                                grid.size, timesteps)
                k_list.append(k)
                logging.info("%d: R squared: %.4f, k: %.4E" % (i * size, r2, k))

        comm.Barrier()

        if (i * 2 * size) > begin_cov_check:
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
        temp_profile = plots.plot_histogram_walkers_2d_onlat(timesteps, tot_H, xedges, yedges, quiet, plot_save_dir,
                                                             gen_plots)
        if gen_plots:
            plots.plot_k_convergence(k_list, quiet, plot_save_dir)
            plots.plot_k_convergence_err(k_convergence_err_list, quiet, plot_save_dir, begin_cov_check)
            temp_gradient_x = plots.plot_temp_gradient_2d_onlat(temp_profile, xedges, yedges, quiet,
                                                            plot_save_dir, gradient_cutoff=2)
        gradient_avg, gradient_std = plots.plot_linear_temp(temp_profile, quiet, plot_save_dir, gen_plots)
        analysis.final_conductivity_2d_onlat(i * 2 * size, grid.size, timesteps, gradient_avg, gradient_std,
                                             k_convergence_err, num_tubes, plot_save_dir, k_convergence_val,
                                             gradient_cutoff=2)
        logging.info("Complete")
