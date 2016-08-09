import logging
import numpy as np
import time
import analysis
import creation
import plots
import randomwalk
import run
from mpi4py import MPI


def sim_2d_onlat(grid_size, tube_length, tube_radius, num_tubes, orientation, timesteps, save_loc_data,
                 quiet, save_loc_plots, save_dir, k_convergence_tolerance, begin_cov_check,
                 k_conv_error_buffer, plot_save_dir, gen_plots, kapitza, prob_m_cn, run_to_convergence, num_walkers):
    walker_data_save_dir = plot_save_dir + "/walker_locations"
    walker_plot_save_dir = plot_save_dir + "/walker_plots"

    grid = creation.Grid2D_onlat(grid_size, tube_length, num_tubes, orientation, tube_radius)
    if gen_plots:
        plots.plot_two_d_random_walk_setup(grid, quiet, plot_save_dir)
        plots.plot_check_array_2d(grid, quiet, plot_save_dir, gen_plots)

    grid_range = [[0, grid.size + 1], [0, grid.size + 1]]
    bins = grid.size + 1
    H = np.zeros((grid.size + 1, grid.size + 1))

    i = 0
    k_list = []
    k_convergence_err_list = []
    k_convergence_err = 1.0

    start = time.clock()
    if run_to_convergence:
        while k_convergence_err > k_convergence_tolerance:
            walker, H, xedges, yedges, i = randomwalk_routine_2d_serial(grid, grid_range, timesteps, save_loc_data,
                                                                        quiet,
                                                                        save_loc_plots, bins, plot_save_dir,
                                                                        walker_plot_save_dir, walker_data_save_dir,
                                                                        gen_plots, kapitza, prob_m_cn, i, H)
            dt_dx, heat_flux, dt_dx_err, k, k_err, r2 = analysis.check_convergence_2d_onlat(H, i * 2,
                                                                                            grid.size, timesteps)
            k_list.append(k)
            logging.info("%d: R squared: %.4f, k: %.4E" % (i, r2, k))
            if i > begin_cov_check:
                k_convergence_err = np.std(np.array(k_list[-k_conv_error_buffer:]), ddof=1)
                k_convergence_val = np.mean(np.array(k_list[-k_conv_error_buffer:]))
                k_convergence_err_list.append(k_convergence_err)
                logging.info("k: %.4E" % k_convergence_val)
                logging.info("k error: %.4E" % k_convergence_err)
    else:
        for i in range(num_walkers / 2):
            walker, H, xedges, yedges, i = randomwalk_routine_2d_serial(grid, grid_range, timesteps, save_loc_data,
                                                                        quiet,
                                                                        save_loc_plots, bins, plot_save_dir,
                                                                        walker_plot_save_dir, walker_data_save_dir,
                                                                        gen_plots, kapitza, prob_m_cn, i, H)
            dt_dx, heat_flux, dt_dx_err, k, k_err, r2 = analysis.check_convergence_2d_onlat(H, i * 2,
                                                                                            grid.size, timesteps)
            k_list.append(k)
            logging.info("%d: R squared: %.4f, k: %.4E" % (i, r2, k))
            k_list.append(k)
            logging.info("%d: R squared: %.4f, k: %.4E" % (i, r2, k))
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
    temp_profile = plots.plot_histogram_walkers_2d_onlat(grid, timesteps, H, xedges, yedges, quiet, plot_save_dir,
                                                         gen_plots)
    if gen_plots == True:
        plots.plot_k_convergence(k_list, quiet, plot_save_dir, begin_cov_check)
        plots.plot_k_convergence_err(k_convergence_err_list, quiet, plot_save_dir, begin_cov_check)
        temp_gradient_x = plots.plot_temp_gradient_2d_onlat(grid, temp_profile, xedges, yedges, quiet,
                                                            plot_save_dir, gradient_cutoff=0)
    gradient_avg, gradient_std = plots.plot_linear_temp(temp_profile, grid_size, quiet, plot_save_dir,
                                                        gen_plots)
    analysis.final_conductivity_2d_onlat(i * 2, grid.size, timesteps, gradient_avg, gradient_std,
                                         k_convergence_err, num_tubes, plot_save_dir, k_convergence_val,
                                         gradient_cutoff=0)
    max_temp = np.max(temp_profile)
    min_temp = np.min(temp_profile)
    diff_temp = np.abs(max_temp) - np.abs(min_temp)
    logging.info("Max temp is %d, min temp is %d, with difference %d" % (max_temp, min_temp, diff_temp))
    logging.info("Complete")


def sim_2d_onlat_MPI(grid_size, tube_length, tube_radius, num_tubes, orientation, timesteps, save_loc_data,
                     quiet, save_loc_plots, save_dir, k_convergence_tolerance, begin_cov_check,
                     k_conv_error_buffer, plot_save_dir, gen_plots, kapitza, prob_m_cn, rank, size,
                     run_to_convergence, num_walkers):

    comm = MPI.COMM_WORLD
    walker_data_save_dir = plot_save_dir + "/walker_locations"
    walker_plot_save_dir = plot_save_dir + "/walker_plots"

    if rank == 0:
        grid = creation.Grid2D_onlat(grid_size, tube_length, num_tubes, orientation, tube_radius)
        if gen_plots:
            plots.plot_two_d_random_walk_setup(grid, quiet, plot_save_dir)
            plots.plot_check_array_2d(grid, quiet, plot_save_dir, gen_plots)
    else:
        grid = None

    grid = comm.bcast(grid, root=0)

    grid_range = [[0, grid.size + 1], [0, grid.size + 1]]
    bins = grid.size + 1
    H = np.zeros((grid.size + 1, grid.size + 1))
    tot_H = np.zeros((grid.size + 1, grid.size + 1))

    i = 0

    k_list = []
    k_convergence_err_list = []
    k_convergence_err = 1.0

    start = MPI.Wtime()

    if run_to_convergence:
        while k_convergence_err > k_convergence_tolerance:
            walker, H, xedges, yedges, i, tot_H = randomwalk_routine_2d_MPI(grid, grid_range, timesteps,
                                                                            save_loc_data, quiet, save_loc_plots, bins,
                                                                            plot_save_dir, walker_plot_save_dir,
                                                                            walker_data_save_dir, gen_plots, kapitza,
                                                                            prob_m_cn, i, H, rank, comm, tot_H)
            if rank == 0:
                dt_dx, heat_flux, dt_dx_err, k, k_err, r2 = analysis.check_convergence_2d_onlat(tot_H, i * 2 * size,
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
    else:
        for i in range(num_walkers / (2 * size)):  # rounds down total walkers slightly
            walker, H, xedges, yedges, i, tot_H = randomwalk_routine_2d_MPI(grid, grid_range, timesteps,
                                                                            save_loc_data, quiet, save_loc_plots, bins,
                                                                            plot_save_dir, walker_plot_save_dir,
                                                                            walker_data_save_dir, gen_plots, kapitza,
                                                                            prob_m_cn, i, H, rank, comm, tot_H)
            if rank == 0:
                dt_dx, heat_flux, dt_dx_err, k, k_err, r2 = analysis.check_convergence_2d_onlat(tot_H, i * 2 * size,
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
        temp_profile = plots.plot_histogram_walkers_2d_onlat(grid, timesteps, tot_H, xedges, yedges, quiet,
                                                             plot_save_dir,
                                                             gen_plots)
        if gen_plots:
            plots.plot_k_convergence(k_list, quiet, plot_save_dir, begin_cov_check)
            plots.plot_k_convergence_err(k_convergence_err_list, quiet, plot_save_dir, begin_cov_check)
            temp_gradient_x = plots.plot_temp_gradient_2d_onlat(grid, temp_profile, xedges, yedges, quiet,
                                                                plot_save_dir, gradient_cutoff=0)
        gradient_avg, gradient_std = plots.plot_linear_temp(temp_profile, grid_size, quiet, plot_save_dir,
                                                            gen_plots)
        analysis.final_conductivity_2d_onlat(i * 2 * size, grid.size, timesteps, gradient_avg, gradient_std,
                                             k_convergence_err, num_tubes, plot_save_dir, k_convergence_val,
                                             gradient_cutoff=0)
        max_temp = np.max(temp_profile)
        min_temp = np.min(temp_profile)
        diff_temp = np.abs(max_temp) - np.abs(min_temp)
        logging.info("Max temp is %d, min temp is %d, with difference %d" % (max_temp, min_temp, diff_temp))
        logging.info("Complete")


def randomwalk_routine_2d_serial(grid, grid_range, timesteps, save_loc_data, quiet, save_loc_plots, bins, plot_save_dir,
                                 walker_plot_save_dir, walker_data_save_dir, gen_plots, kapitza, prob_m_cn,
                                 i, H):
    # run hot walker
    # logging.info("Start hot walker %d" % (i+1))
    walker = randomwalk.runrandomwalk_2d_onlat(grid, timesteps, 'hot', kapitza, prob_m_cn)
    if save_loc_data:
        run.save_walker_loc(walker, walker_data_save_dir, i, 'hot')
    if i == 0 & save_loc_plots == False & gen_plots == True:  # always save one example trajectory plot
        plots.plot_walker_path_2d_onlat(walker, grid.size, 'hot', quiet, i + 1, plot_save_dir)
    elif save_loc_plots & gen_plots == True:
        plots.plot_walker_path_2d_onlat(walker, grid.size, 'hot', quiet, i + 1, walker_plot_save_dir)
    H_temp, xedges, yedges = plots.histogram_walker_2d_onlat(walker, grid_range, bins)
    H += H_temp

    # run cold walker
    # logging.info("Start cold walker %d" % (i+1))
    walker = randomwalk.runrandomwalk_2d_onlat(grid, timesteps, 'cold', kapitza, prob_m_cn)
    if save_loc_data:
        run.save_walker_loc(walker, walker_data_save_dir, i, 'cold')
    if i == 0 & save_loc_plots == False & gen_plots == True:
        plots.plot_walker_path_2d_onlat(walker, grid.size, 'cold', quiet, i + 1, plot_save_dir)
    elif save_loc_plots & gen_plots == True:
        plots.plot_walker_path_2d_onlat(walker, grid.size, 'cold', quiet, i + 1, walker_plot_save_dir)
    H_temp, xedges, yedges = plots.histogram_walker_2d_onlat(walker, grid_range, bins)
    H -= H_temp

    i += 1

    return walker, H, xedges, yedges, i


def randomwalk_routine_2d_MPI(grid, grid_range, timesteps, save_loc_data, quiet, save_loc_plots, bins, plot_save_dir,
                              walker_plot_save_dir, walker_data_save_dir, gen_plots, kapitza, prob_m_cn,
                              i, H, rank, comm, tot_H):
    # run hot walker
    # logging.info("Start hot walker %d" % (i+1))
    walker = randomwalk.runrandomwalk_2d_onlat(grid, timesteps, 'hot', kapitza, prob_m_cn)

    if rank == 0:
        if save_loc_data:
            run.save_walker_loc(walker, walker_data_save_dir, i, 'hot')
        if i == 0 & save_loc_plots == False & gen_plots == True:  # always save one example trajectory plot
            plots.plot_walker_path_2d_onlat(walker, grid.size, 'hot', quiet, i + 1, plot_save_dir)
        elif save_loc_plots & gen_plots == True:
            plots.plot_walker_path_2d_onlat(walker, grid.size, 'hot', quiet, i + 1, walker_plot_save_dir)

    H_temp, xedges, yedges = plots.histogram_walker_2d_onlat(walker, grid_range, bins)
    H += H_temp

    # run cold walker
    # logging.info("Start cold walker %d" % (i+1))
    walker = randomwalk.runrandomwalk_2d_onlat(grid, timesteps, 'cold', kapitza, prob_m_cn)

    if rank == 0:
        if save_loc_data:
            run.save_walker_loc(walker, walker_data_save_dir, i, 'cold')
        if i == 0 & save_loc_plots == False & gen_plots == True:
            plots.plot_walker_path_2d_onlat(walker, grid.size, 'cold', quiet, i + 1, plot_save_dir)
        elif save_loc_plots & gen_plots == True:
            plots.plot_walker_path_2d_onlat(walker, grid.size, 'cold', quiet, i + 1, walker_plot_save_dir)

    H_temp, xedges, yedges = plots.histogram_walker_2d_onlat(walker, grid_range, bins)
    H -= H_temp

    comm.Reduce(H, tot_H, op=MPI.SUM, root=0)
    # H is updated on every core for every i independently
    # tot_H is the total across all cores

    comm.Barrier()
    i += 1  # i starts at 0

    return walker, H, xedges, yedges, i, tot_H
