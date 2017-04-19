from __future__ import division
from builtins import range
from past.utils import old_div
import logging
import numpy as np
import time
from mpi4py import MPI
from conduction import *


def serial_method(grid_size, tube_length, tube_radius, num_tubes, orientation, timesteps, save_loc_data,
                  quiet, save_loc_plots, save_dir, k_convergence_tolerance, begin_cov_check,
                  k_conv_error_buffer, plot_save_dir, gen_plots, kapitza, prob_m_cn, run_to_convergence,
                  num_walkers):
    walker_data_save_dir = plot_save_dir + "/walker_locations"
    walker_plot_save_dir = plot_save_dir + "/walker_plots"

    grid = creation.Grid2D_onlat(grid_size, tube_length, num_tubes, orientation, tube_radius)
    if gen_plots:
        plots.plot_two_d_random_walk_setup(grid, quiet, plot_save_dir)
        # plots.plot_check_array_2d(grid, quiet, plot_save_dir, gen_plots)

    grid_range = [[0, grid.size + 1], [0, grid.size + 1]]
    bins = grid.size + 1
    H = np.zeros((grid.size + 1, grid.size + 1), dtype=int)

    i = 0
    k_list = []
    dt_dx_list = []
    heat_flux_list = []
    k_convergence_err_list = []
    k_convergence_err = 1.0
    xedges = list(range(0, bins))
    yedges = list(range(0, bins))

    start = time.clock()

    if run_to_convergence:
        while k_convergence_err > k_convergence_tolerance:
            H = randomwalk_routine_2d_serial(grid, timesteps, save_loc_data, quiet, save_loc_plots,
                                             plot_save_dir, walker_plot_save_dir, walker_data_save_dir,
                                             gen_plots, kapitza, prob_m_cn, i, H)

            i += 1

            dt_dx, heat_flux, dt_dx_err, k, k_err, r2 = analysis.check_convergence_2d_onlat(H, i * 2,
                                                                                            grid.size, timesteps)
            k_list.append(k)
            dt_dx_list.append(dt_dx)
            heat_flux_list.append(heat_flux)
            logging.info("%d: R squared: %.4f, k: %.4E, dT/dx: %.4E" % (i, r2, k, dt_dx))
            if i > begin_cov_check:
                k_convergence_err = np.std(np.array(k_list[-k_conv_error_buffer:]), ddof=1)
                k_convergence_val = np.mean(np.array(k_list[-k_conv_error_buffer:]))
                k_convergence_err_list.append(k_convergence_err)
                logging.info("k: %.4E" % k_convergence_val)
                logging.info("k error: %.4E" % k_convergence_err)
    else:
        for i in range(old_div(num_walkers, 2)):
            H = randomwalk_routine_2d_serial(grid, timesteps, save_loc_data, quiet, save_loc_plots, plot_save_dir,
                                             walker_plot_save_dir, walker_data_save_dir,
                                             gen_plots, kapitza, prob_m_cn, i, H)

            dt_dx, heat_flux, dt_dx_err, k, k_err, r2 = analysis.check_convergence_2d_onlat(H, i * 2,
                                                                                            grid.size, timesteps)
            k_list.append(k)
            dt_dx_list.append(dt_dx)
            heat_flux_list.append(heat_flux)
            logging.info("%d: R squared: %.4f, k: %.4E, dT/dx: %.4E" % (i, r2, k, dt_dx))
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
    walk_sec = old_div((i * 2), (end - start))
    logging.info("Crunched %.4f walkers/second" % walk_sec)
    temp_profile = plots.plot_colormap_2d(grid, H, quiet, plot_save_dir, gen_plots)
    if gen_plots == True:
        plots.plot_k_convergence(k_list, quiet, plot_save_dir)
        plots.plot_k_convergence_err(k_convergence_err_list, quiet, plot_save_dir, begin_cov_check)
        plots.plot_dt_dx(dt_dx_list, quiet, plot_save_dir)
        plots.plot_heat_flux(heat_flux_list, quiet, plot_save_dir)
        temp_gradient_x = plots.plot_temp_gradient_2d_onlat(grid, temp_profile, xedges, yedges, quiet,
                                                            plot_save_dir, gradient_cutoff=0)
    gradient_avg, gradient_std = plots.plot_linear_temp(temp_profile, grid_size, quiet, plot_save_dir,
                                                        gen_plots)
    analysis.final_conductivity_2d_onlat(i * 2, grid.size, timesteps, gradient_avg, gradient_std,
                                         k_convergence_err, num_tubes, plot_save_dir, k_convergence_val,
                                         prob_m_cn, gradient_cutoff=0)
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
            #plots.plot_check_array_2d(grid, quiet, plot_save_dir, gen_plots)
    else:
        grid = None

    grid = comm.bcast(grid, root=0)

    grid_range = [[0, grid.size + 1], [0, grid.size + 1]]
    bins = grid.size + 1
    H = np.zeros((grid.size + 1, grid.size + 1), dtype=int)
    tot_H = np.zeros((grid.size + 1, grid.size + 1), dtype=int)

    i = 0

    k_list = []
    dt_dx_list = []
    heat_flux_list = []
    k_convergence_err_list = []
    k_convergence_err = 1.0
    xedges = list(range(0, bins))
    yedges = list(range(0, bins))

    start = MPI.Wtime()

    if run_to_convergence:
        while k_convergence_err > k_convergence_tolerance:
            tot_H = randomwalk_routine_2d_MPI(grid, timesteps, save_loc_data, quiet, save_loc_plots, plot_save_dir,
                                              walker_plot_save_dir, walker_data_save_dir, gen_plots, kapitza,
                                              prob_m_cn, i, H, rank, comm, tot_H)

            i += 1

            if rank == 0:
                dt_dx, heat_flux, dt_dx_err, k, k_err, r2 = analysis.check_convergence_2d_onlat(tot_H, i * 2 * size,
                                                                                                grid.size,
                                                                                                timesteps)
                k_list.append(k)
                dt_dx_list.append(dt_dx)
                heat_flux_list.append(heat_flux)
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
        for i in range(old_div(num_walkers, (2 * size))):  # rounds down total walkers slightly
            tot_H = randomwalk_routine_2d_MPI(grid, timesteps, save_loc_data, quiet, save_loc_plots, plot_save_dir,
                                              walker_plot_save_dir, walker_data_save_dir, gen_plots, kapitza,
                                              prob_m_cn, i, H, rank, comm, tot_H)
            if rank == 0:
                dt_dx, heat_flux, dt_dx_err, k, k_err, r2 = analysis.check_convergence_2d_onlat(tot_H, i * 2 * size,
                                                                                                grid.size,
                                                                                                timesteps)
                k_list.append(k)
                dt_dx_list.append(dt_dx)
                heat_flux_list.append(heat_flux)
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
        walk_sec = old_div((i * 2 * size), (end - start))
        logging.info("Crunched %.4f walkers/second" % walk_sec)
        temp_profile = plots.plot_colormap_2d(grid, tot_H, quiet, plot_save_dir, gen_plots)
        if gen_plots:
            plots.plot_k_convergence(k_list, quiet, plot_save_dir)
            plots.plot_k_convergence_err(k_convergence_err_list, quiet, plot_save_dir, begin_cov_check)
            plots.plot_dt_dx(dt_dx_list, quiet, plot_save_dir)
            plots.plot_heat_flux(heat_flux_list, quiet, plot_save_dir)
            temp_gradient_x = plots.plot_temp_gradient_2d_onlat(grid, temp_profile, xedges, yedges, quiet,
                                                                plot_save_dir, gradient_cutoff=0)
        gradient_avg, gradient_std = plots.plot_linear_temp(temp_profile, grid_size, quiet, plot_save_dir,
                                                            gen_plots)
        analysis.final_conductivity_2d_onlat(i * 2 * size, grid.size, timesteps, gradient_avg, gradient_std,
                                             k_convergence_err, num_tubes, plot_save_dir, k_convergence_val,
                                             prob_m_cn, gradient_cutoff=0)
        max_temp = np.max(temp_profile)
        min_temp = np.min(temp_profile)
        diff_temp = np.abs(max_temp) - np.abs(min_temp)
        logging.info("Max temp is %d, min temp is %d, with difference %d" % (max_temp, min_temp, diff_temp))
        logging.info("Complete")


def randomwalk_routine_2d_serial(grid, timesteps, save_loc_data, quiet, save_loc_plots, plot_save_dir,
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

    H[walker.pos[-1][0], walker.pos[-1][1]] += 1

    # run cold walker
    # logging.info("Start cold walker %d" % (i+1))
    walker = randomwalk.runrandomwalk_2d_onlat(grid, timesteps, 'cold', kapitza, prob_m_cn)
    if save_loc_data:
        run.save_walker_loc(walker, walker_data_save_dir, i, 'cold')
    if i == 0 & save_loc_plots == False & gen_plots == True:
        plots.plot_walker_path_2d_onlat(walker, grid.size, 'cold', quiet, i + 1, plot_save_dir)
    elif save_loc_plots & gen_plots == True:
        plots.plot_walker_path_2d_onlat(walker, grid.size, 'cold', quiet, i + 1, walker_plot_save_dir)
    H[walker.pos[-1][0], walker.pos[-1][1]] -= 1

    return H


def randomwalk_routine_2d_MPI(grid, timesteps, save_loc_data, quiet, save_loc_plots, plot_save_dir,
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

    H[walker.pos[-1][0], walker.pos[-1][1]] += 1

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

    H[walker.pos[-1][0], walker.pos[-1][1]] -= 1

    comm.Reduce(H, tot_H, op=MPI.SUM, root=0)
    # H is updated on every core for every i independently
    # tot_H is the total across all cores

    comm.Barrier()

    return tot_H
