import logging
import numpy as np
import time
import analysis
import creation
import plots
import randomwalk
import run
from scipy import stats
from mpi4py import MPI


def sim_2d_onlat_constant_flux(grid_size, tube_length, tube_radius, num_tubes, orientation, timesteps, save_loc_data,
                               quiet, save_loc_plots, save_dir, k_convergence_tolerance, begin_cov_check,
                               k_conv_error_buffer, plot_save_dir, gen_plots, kapitza, prob_m_cn, run_to_convergence,
                               num_walkers, method, printout_inc):
    walker_data_save_dir = plot_save_dir + "/walker_locations"
    walker_plot_save_dir = plot_save_dir + "/walker_plots"

    grid = creation.Grid2D_onlat(grid_size, tube_length, num_tubes, orientation, tube_radius)
    if gen_plots:
        plots.plot_two_d_random_walk_setup(grid, quiet, plot_save_dir)
        plots.plot_check_array_2d(grid, quiet, plot_save_dir, gen_plots)

    grid_range = [[0, grid.size + 1], [0, grid.size + 1]]
    bins = grid.size + 1
    H = np.zeros((grid.size + 1, grid.size + 1))

    start = time.clock()

    hot_walker_master, cold_walker_master, H, xedges, yedges, k_list, dt_dx_list, heat_flux_list, timestep_list = \
        randomwalk_routine_2d_serial(grid, grid_range, timesteps, save_loc_data, quiet, save_loc_plots, bins,
                                     plot_save_dir, walker_plot_save_dir, walker_data_save_dir, gen_plots, kapitza,
                                     prob_m_cn, H, num_walkers,
                                     printout_inc)

    dt_dx, heat_flux, dt_dx_err, k, k_err, r2 = analysis.check_convergence_2d_onlat(H, num_walkers,
                                                                                    grid.size, timesteps)
    logging.info("%d walkers: R squared: %.4f, k: %.4E, heat flux: %.4E" % (num_walkers, r2, k, heat_flux))

    end = time.clock()
    logging.info("Constant flux simulation has completed")
    logging.info("Serial simulation time was %.4f s" % (end - start))

    temp_profile = plots.plot_histogram_walkers_2d_onlat(grid, timesteps, H, xedges, yedges, quiet, plot_save_dir,
                                                         gen_plots)
    if gen_plots:
        plots.plot_k_convergence(k_list, quiet, plot_save_dir, timestep_list)
        plots.plot_dt_dx(dt_dx_list, quiet, plot_save_dir, timestep_list)
        plots.plot_heat_flux(heat_flux_list, quiet, plot_save_dir, timestep_list)
        temp_gradient_x = plots.plot_temp_gradient_2d_onlat(grid, temp_profile, xedges, yedges, quiet,
                                                            plot_save_dir, gradient_cutoff=0)
        gradient_avg, gradient_std = plots.plot_linear_temp(temp_profile, grid_size, quiet, plot_save_dir,
                                                            gen_plots)
    logging.info("Complete")


def randomwalk_routine_2d_serial(grid, grid_range, tot_time, save_loc_data, quiet, save_loc_plots, bins, plot_save_dir,
                                 walker_plot_save_dir, walker_data_save_dir, gen_plots, kapitza, prob_m_cn,
                                 H, tot_walkers, printout_inc):
    def histogram_walker_list(hot_walker_master, cold_walker_master, H, grid_range, bins):
        H = np.zeros((grid.size + 1, grid.size + 1))  # resets H every time function is called, IMPORTANT
        for i in range(len(hot_walker_master)):
            hot_temp = hot_walker_master[i]
            H_temp, xedges, yedges = plots.histogram_walker_2d_onlat(hot_temp, grid_range, bins)
            H += H_temp

            cold_temp = cold_walker_master[i]
            H_temp, xedges, yedges = plots.histogram_walker_2d_onlat(cold_temp, grid_range, bins)
            H -= H_temp
        return H, xedges, yedges

    # these will hold the walker objects, need to access walker.pos for the positions
    hot_walker_master = []
    cold_walker_master = []
    k_list = []
    dt_dx_list = []
    heat_flux_list = []
    timestep_list = []  # x axis for plots

    d_add = tot_time / (tot_walkers / 2.0)  # how often to add hot/cold walkers

    for i in range(tot_time):
        if (i % printout_inc) == 0 and (i > (5 * printout_inc)):
            cur_num_walkers = len(hot_walker_master) * 2
            Htemp, xedges, yedges = histogram_walker_list(hot_walker_master, cold_walker_master, H, grid_range, bins)
            dt_dx, heat_flux, dt_dx_err, k, k_err, r2 = analysis.check_convergence_2d_onlat(Htemp, cur_num_walkers,
                                                                                            grid.size, i)
            k_list.append(k)
            dt_dx_list.append(dt_dx)
            heat_flux_list.append(heat_flux)
            timestep_list.append(i)

            logging.info("Timestep: %d, %d walkers, R2: %.4f, k: %.4E, heat flux: %.4E, dT(x)/dx: %.4E"
                         % (i, cur_num_walkers, r2, k, heat_flux, dt_dx))
        if (i % d_add) == 0:
            trigger = 1
            # let's add the 2 new walkers
            hot_walker_temp = creation.Walker2D_onlat(grid.size, 'hot')
            cold_walker_temp = creation.Walker2D_onlat(grid.size, 'cold')
            hot_walker_master.append(hot_walker_temp)
            cold_walker_master.append(cold_walker_temp)
        else:
            trigger = 0
        # let's update all the positions of the activated walkers
        for j in range(len(hot_walker_master) - trigger):  # except the new ones
            hot_temp = hot_walker_master[j]
            current_hot_updated = randomwalk.apply_moves_2d(hot_temp, kapitza, grid, prob_m_cn)
            hot_walker_master[j] = current_hot_updated

            cold_temp = cold_walker_master[j]
            current_cold_updated = randomwalk.apply_moves_2d(cold_temp, kapitza, grid, prob_m_cn)
            cold_walker_master[j] = current_cold_updated

    # let's histogram everything
    logging.info('Finished random walks, histogramming...')

    H, xedges, yedges = histogram_walker_list(hot_walker_master, cold_walker_master, H, grid_range, bins)

    return hot_walker_master, cold_walker_master, H, xedges, yedges, k_list, dt_dx_list, heat_flux_list, timestep_list
