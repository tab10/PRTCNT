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


"""test_3d.py - Tim Burt 3/24/2017
PRTCNT package

This file is designed to test the random walk rules to make sure they satisfy detailed balance, or that the
system has microscopic reversibility and locations are equally probable as they must be in equilibrium"""

from __future__ import division
from builtins import range
from past.utils import old_div
import logging
import numpy as np
import time
from mpi4py import MPI
from conduction import *


def serial_method(grid_size, tube_length, tube_radius, num_tubes, orientation, tot_time, quiet, plot_save_dir,
                  gen_plots, kapitza, prob_m_cn, tot_walkers, printout_inc, k_conv_error_buffer, disable_func,
                  rules_test, bound):
    def histogram_walker_list(H, walker_master):
        # DOES NOT resets H every time function is called, IMPORTANT
        for i in range(len(walker_master)):
            walker_temp = walker_master[i]
            H[walker_temp.pos[-1][0], walker_temp.pos[-1][1], walker_temp.pos[-1][2]] += 1
        # DEBUG
        # print np.sum(H), np.max(H), np.min(H)
        return H

    grid = creation.Grid3D_onlat(grid_size, tube_length, num_tubes, orientation, tube_radius, False, plot_save_dir,
                                 disable_func, rules_test)
    if gen_plots:
        plots.plot_three_d_random_walk_setup(grid, quiet, plot_save_dir)

    grid_range = [[0, grid.size + 1], [0, grid.size + 1], [0, grid.size + 1]]
    bins = grid.size + 1

    start = time.clock()

    # these will hold the walker objects, need to access walker.pos for the positions
    H = np.zeros((grid.size + 1, grid.size + 1, grid.size + 1))
    walker_master = []
    k_list = []
    k_err_list = []
    dt_dx_list = []
    heat_flux_list = []
    timestep_list = []  # x axis for plots
    xedges = list(range(0, bins))
    yedges = list(range(0, bins))
    zedges = list(range(0, bins))
    start_k_err_check = old_div(tot_time, 2)
    prev_type = 0  # assume walker was in matrix before 1st step. this tells the previous type of
    # cell the walker was on

    # d_add - how often to add a walker
    d_add = tot_time / tot_walkers  # as a float
    if d_add.is_integer() and d_add >= 1:
        d_add = int(tot_time, tot_walkers)
        walker_frac_trigger = 0  # add a walker every d_add timesteps
    elif d_add < 1:  # this is a fractional number < 1, implies more than 1 walker should be added every timestep
        d_add = int(old_div(1.0, d_add))
        walker_frac_trigger = 1
    else:  # change num_walkers or timesteps
        logging.error('Choose tot_time / (tot_walkers) so that it is integer or less than 1')
        raise SystemExit
    if walker_frac_trigger == 1:
        logging.info('Adding %d walker(s) every timestep, this might not converge' % d_add)
    elif walker_frac_trigger == 0:
        logging.info('Adding 1 walker(s) every %d timesteps' % d_add)

    for i in range(tot_time):
        if (i % printout_inc) == 0 and (i >= (2 * printout_inc)):
            cur_num_walkers = len(walker_master)
            Htemp = histogram_walker_list(H, walker_master)
            dt_dx, heat_flux, dt_dx_err, k, k_err, r2, temp_profile_sum = analysis.check_convergence_3d_onlat(Htemp,
                                                                                                              cur_num_walkers,
                                                                                                              grid.size,
                                                                                                              i)
            k_list.append(k)
            dt_dx_list.append(dt_dx)
            heat_flux_list.append(heat_flux)
            timestep_list.append(i)
            if i >= start_k_err_check:
                k_err = np.std(k_list[-k_conv_error_buffer:], ddof=1)
                k_err_list.append(k_err)
                logging.info("Timestep: %d, %d walkers, R2: %.4f, k: %.4E, k err: %.4E, heat flux: %.4E, dT(x)/dx: %.4E"
                             % (i, cur_num_walkers, r2, k, k_err, heat_flux, dt_dx))
            else:
                logging.info("Timestep: %d, %d walkers, R2: %.4f, k: %.4E, heat flux: %.4E, dT(x)/dx: %.4E"
                             % (i, cur_num_walkers, r2, k, heat_flux, dt_dx))
        if walker_frac_trigger == 0:
            if (i % d_add) == 0:
                trigger = 1
                # let's add the 1 new walker
                walker_temp = creation.Walker3D_onlat(grid.size, 'hot', rules_test)
                walker_master.append(walker_temp)
            else:
                trigger = 0
        elif walker_frac_trigger == 1:
            trigger = d_add  # trigger never 0 in this case
            # let's add the 1 new walkers, several times
            for k in range(d_add):
                walker_temp = creation.Walker3D_onlat(grid.size, 'hot', rules_test)
                walker_master.append(walker_temp)
        # let's update all the positions of the activated walkers
        for j in range(len(walker_master) - trigger):  # except the new ones
            temp = walker_master[j]
            current_updated, prev_type_hot = rules_3d.apply_moves_3d(temp, kapitza, grid, prob_m_cn, prev_type_hot)(
                walker, kapitza, grid, prob_m_cn, inside_cnt, bound)
            current_updated.erase_prev_pos()
            walker_master[j] = current_updated

    # let's histogram everything
    logging.info('Finished random walks, histogramming...')

    H = histogram_walker_list(hot_walker_master, cold_walker_master)

    dt_dx, heat_flux, dt_dx_err, k, k_err, r2, temp_profile_sum = analysis.check_convergence_3d_onlat(H, tot_walkers,
                                                                                                      grid.size,
                                                                                                      tot_time)
    logging.info("%d walkers: R squared: %.4f, k: %.4E, heat flux: %.4E, dT(x)/dx: %.4E" % (tot_walkers, r2, k,
                                                                                            heat_flux, dt_dx))

    end = time.clock()
    logging.info("Constant flux simulation has completed")
    logging.info("Serial simulation time was %.4f s" % (end - start))

    temp_profile = plots.plot_colormap_2d(grid, tot_time, temp_profile_sum, xedges, yedges, quiet,
                                          plot_save_dir,
                                          gen_plots)
    if gen_plots:
        plots.plot_k_convergence(k_list, quiet, plot_save_dir, timestep_list)
        plots.plot_k_convergence_err(k_list, quiet, plot_save_dir, start_k_err_check, timestep_list)
        plots.plot_dt_dx(dt_dx_list, quiet, plot_save_dir, timestep_list)
        plots.plot_heat_flux(heat_flux_list, quiet, plot_save_dir, timestep_list)
        temp_gradient_x = plots.plot_temp_gradient_2d_onlat(grid, temp_profile, xedges, yedges, quiet,
                                                            plot_save_dir, gradient_cutoff=0)
        gradient_avg, gradient_std = plots.plot_linear_temp(temp_profile, grid_size, quiet, plot_save_dir,
                                                            gen_plots)
    logging.info("Complete")


def parallel_method(grid_size, tube_length, tube_radius, num_tubes, orientation, tot_time, quiet, plot_save_dir,
                    gen_plots, kapitza, prob_m_cn, tot_walkers, disable_func, rank, size, rules_test, restart):
    comm = MPI.COMM_WORLD

    # serial tube generation
    if rank == 0:
        grid = creation.Grid3D_onlat(grid_size, tube_length, num_tubes, orientation, tube_radius, False, plot_save_dir,
                                     disable_func, rules_test)

    comm.Barrier()

    if rank == 0:
        if gen_plots:
            plots.plot_three_d_random_walk_setup(grid, quiet, plot_save_dir)
    else:
        grid = None

    comm.Barrier()
    grid = comm.bcast(grid, root=0)

    grid_range = [[0, grid.size + 1], [0, grid.size + 1], [0, grid.size + 1]]
    bins = grid.size + 1

    start = MPI.Wtime()

    # these will hold the walker objects, need to access walker.pos for the positions

    walker_master_pos = []
    timestep_list = []  # x axis for plots
    x_edges = list(range(0, bins))
    y_edges = list(range(0, bins))
    z_edges = list(range(0, bins))
    # start_k_err_check = old_div(tot_time, 2)

    # d_add - how often to add a hot/cold walker pair
    # d_add = old_div(tot_time, (old_div(tot_walkers, 2.0)))  # as a float
    # if d_add.is_integer() and d_add >= 1:
    #     d_add = int(tot_time / tot_walkers)
    #     walker_frac_trigger = 0  # add a pair every d_add timesteps
    # elif d_add < 1:  # this is a fractional number < 1, implies more than 1 walker pair should be added every timestep
    #     d_add = old_div(1.0, d_add)
    #     if not d_add.is_integer():
    #         logging.error('Choose tot_time / tot_walkers so that it is integer or less than 1 and whole')
    #         raise SystemExit
    #     else:
    #         d_add = int(d_add)
    #     walker_frac_trigger = 1
    # else:  # change num_walkers or timesteps
    #     logging.error('Choose tot_time / tot_walkers so that it is integer or less than 1')
    #     raise SystemExit
    # if walker_frac_trigger == 1:
    #     logging.info('Adding %d walkers pair(s) every timestep' % d_add)
    #     walkers_per_core_whole = int(np.floor(old_div(tot_walkers, (size * d_add))))
    # elif walker_frac_trigger == 0:
    #     logging.info(
    #         'Adding 1 walker(s) every %d timesteps. Likely will not have enough walkers.' % d_add)
    # walkers_per_core_whole = int(np.floor(old_div(tot_walkers, (size))))

    walkers_per_core_whole = int(np.floor(old_div(tot_walkers, (size))))

    comm.Barrier()

    walkers_per_core_remain = int(tot_walkers % size)
    if walkers_per_core_remain != 0:
        logging.error('Algorithm cannot currently handle a remainder between tot_walkers and tot_cores')
        raise SystemExit

    comm.Barrier()

    H_local = np.zeros((grid.size + 1, grid.size + 1, grid.size + 1), dtype=int)
    for i in range(walkers_per_core_whole):
        H_master = np.zeros((grid.size + 1, grid.size + 1, grid.size + 1), dtype=int)  # should be reset every iteration
        # if walker_frac_trigger == 0:
        #     core_time = ((i * size) + rank) * d_add
        #     cur_num_walkers = i * size
        #     walkers_per_timestep = 1
        # elif walker_frac_trigger == 1:
        #     core_time = i * size + rank
        #     cur_num_walkers = i * size * d_add
        #     walkers_per_timestep = d_add
        cur_num_walkers = i * size
        walker_temp = rules_3d.runrandomwalk_3d_onlat(grid, tot_time, 'hot', kapitza, prob_m_cn, grid.bound,
                                                          rules_test)  # 'hot' since +1
        # histogram ALL positions
        for j in range(len(walker_temp.pos)):
            H_local[walker_temp.pos[j][0], walker_temp.pos[j][1], walker_temp.pos[j][2]] += 1
        # send to core 0
        # as long as size is somewhat small, this barrier won't slow things down much
        # H_local now has all walker positions.
        comm.Barrier()
        comm.Reduce(H_local, H_master, op=MPI.SUM, root=0)
        # analysis
        if rank == 0 and (i > 0):
            # print np.count_nonzero(H_master)
            timestep_list.append(tot_time)
            logging.info("Parallel iteration %d out of %d, %d walkers"
                         % (i, walkers_per_core_whole, cur_num_walkers))
        comm.Barrier()

    comm.Barrier()  # make sure whole walks are done

    if rank == 0:
        logging.info('Finished random walk rules test. Analyzing....')
        mean_temp, mean_temp_norm, std_temp, std_temp_norm = analysis.rules_test_analysis(H_local, cur_num_walkers)
        logging.info('Histogram: mean-%.4E, std-%.4E' % (mean_temp, std_temp))
        logging.info('Histogram normalized: mean-%.4E, std-%.4E' % (mean_temp_norm, std_temp_norm))
        # plot 2d cuts in each of the 3 planes, all directions periodic
        temp_profile_yz = np.sum(H_local, axis=0)
        temp_profile_xz = np.sum(H_local, axis=1)
        temp_profile_xy = np.sum(H_local, axis=2)
        plots.plot_colormap_2d(grid, temp_profile_xy, quiet, plot_save_dir, gen_plots,
                               title='Temperature density (dimensionless units)',
                               xlab='X', ylab='Y', filename='temp_xy')
        # plots.plot_bargraph_3d(grid, temp_profile_xy, x_edges, y_edges, quiet, save_dir, gen_plots,
        #                 title='Temperature density (dimensionless units)', xlab='X', ylab='Y', zlab='Z',
        #                 filename='temp')
        end = MPI.Wtime()
        logging.info("Constant flux simulation has completed")
        logging.info("Using %d cores, parallel simulation time was %.4f min" % (size, old_div((end - start), 60.0)))
        walk_sec = old_div(tot_walkers, (end - start))
        logging.info("Crunched %.4f walkers/second" % walk_sec)
        # temp_profile = plots.plot_colormap_2d(grid, tot_time, temp_profile_sum, xedges, yedges, quiet,
        #                                      plot_save_dir, gen_plots)
        # if gen_plots:
        #     plots.plot_k_convergence(k_list, quiet, plot_save_dir, timestep_list)
        #     plots.plot_k_convergence_err(k_list, quiet, plot_save_dir, start_k_err_check, timestep_list)
        #     plots.plot_dt_dx(dt_dx_list, quiet, plot_save_dir, timestep_list)
        #     plots.plot_heat_flux(heat_flux_list, quiet, plot_save_dir, timestep_list)
        #     temp_gradient_x = plots.plot_temp_gradient_2d_onlat(grid, temp_profile, xedges, yedges, quiet,
        #                                                         plot_save_dir, gradient_cutoff=0)
        #     gradient_avg, gradient_std = plots.plot_linear_temp(temp_profile, grid_size, quiet, plot_save_dir,
        #                                                         gen_plots)
        logging.info("Complete")
