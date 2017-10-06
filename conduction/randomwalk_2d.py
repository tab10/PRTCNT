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


from __future__ import division
import logging
import numpy as np
from mpi4py import MPI

from conduction import creation_2d
from conduction import plots
from conduction import rules_2d
from conduction import analysis


def parallel_method(grid_size, tube_length, tube_radius, num_tubes, orientation, tot_time, quiet, plot_save_dir,
                    gen_plots, kapitza, prob_m_cn, tot_walkers, printout_inc, k_conv_error_buffer, disable_func, rank,
                    size, rules_test, restart, inert_vol, save_loc_plots):
    comm = MPI.COMM_WORLD

    # serial tube generation
    if rank == 0:
        grid = creation_2d.Grid2D_onlat(grid_size, tube_length, num_tubes, orientation, tube_radius, False,
                                        plot_save_dir,
                                        disable_func, rules_test, inert_vol)

    comm.Barrier()

    if rank == 0:
        if gen_plots:
            plots.plot_two_d_random_walk_setup(grid, quiet, plot_save_dir, inert_vol)
            #plots.plot_check_array_2d(grid, quiet, plot_save_dir, gen_plots)
    else:
        grid = None

    comm.Barrier()
    grid = comm.bcast(grid, root=0)

    grid_range = [[0, grid.size + 1], [0, grid.size + 1]]
    bins = grid.size + 1

    start = MPI.Wtime()

    # these will hold the walker objects, need to access walker.pos for the positions

    hot_walker_master_pos = []
    cold_walker_master_pos = []
    k_list = []
    k_err_list = []
    dt_dx_list = []
    heat_flux_list = []
    timestep_list = []  # x axis for plots
    xedges = list(range(0, bins))
    yedges = list(range(0, bins))
    start_k_err_check = tot_time / 2

    # d_add - how often to add a hot/cold walker pair
    d_add = tot_time / (tot_walkers / 2.0)  # as a float
    # print d_add
    if d_add.is_integer() and d_add >= 1:
        d_add = int(tot_time / (tot_walkers / 2.0))
        walker_frac_trigger = 0  # add a pair every d_add timesteps
    elif d_add < 1:  # this is a fractional number < 1, implies more than 1 walker pair should be added every timestep
        d_add = 1.0 / d_add
        if not d_add.is_integer():
            logging.error('Choose tot_time / (tot_walkers / 2.0) so that it is integer or less than 1 and whole')
            raise SystemExit
        else:
            d_add = int(d_add)
        walker_frac_trigger = 1
    else:  # change num_walkers or timesteps
        logging.error('Choose tot_time / (tot_walkers / 2.0) so that it is integer or less than 1 and whole')
        raise SystemExit
    if walker_frac_trigger == 1:
        logging.info('Adding %d hot/cold walker pair(s) every timestep' % d_add)
        walkers_per_core_whole = int(np.floor(tot_walkers / (2.0 * size * d_add)))
    elif walker_frac_trigger == 0:
        logging.info(
            'Adding 1 hot/cold walker pair(s) every %d timesteps. Likely will not have enough walkers.' % d_add)
        walkers_per_core_whole = int(np.floor(tot_walkers / (2.0 * size)))

    comm.Barrier()

    walkers_per_core_remain = int(tot_walkers % size)
    if walkers_per_core_remain != 0:
        logging.error('Algorithm cannot currently handle a remainder between tot_walkers and tot_cores')
        raise SystemExit

    H_local = np.zeros((grid.size + 1, grid.size + 1), dtype=int)

    comm.Barrier()

    for i in range(walkers_per_core_whole):
        H_master = np.zeros((grid.size + 1, grid.size + 1), dtype=int)  # should be reset every iteration
        if walker_frac_trigger == 0:
            core_time = ((i * size) + rank) * d_add
            cur_num_walkers = 2 * i * size
            walkers_per_timestep = 1
        elif walker_frac_trigger == 1:
            core_time = i * size + rank
            cur_num_walkers = 2 * i * size * d_add
            walkers_per_timestep = d_add
        for j in range(walkers_per_timestep):
            # print '%d on core %d' % (core_time, rank)
            # run trajectories for that long
            hot_temp = rules_2d.runrandomwalk_2d_onlat(grid, core_time, 'hot', kapitza, prob_m_cn, grid.bound,
                                                       rules_test)
            cold_temp = rules_2d.runrandomwalk_2d_onlat(grid, core_time, 'cold', kapitza, prob_m_cn, grid.bound,
                                                        rules_test)
            # plot walker path if desired
            if save_loc_plots:
                plots.plot_walker_path_2d_onlat(hot_temp, grid_size, 'hot', quiet, i, plot_save_dir)
                plots.plot_walker_path_2d_onlat(cold_temp, grid_size, 'cold', quiet, i, plot_save_dir)
            # get last position of walker
            hot_temp_pos = hot_temp.pos[-1]
            cold_temp_pos = cold_temp.pos[-1]
            # histogram
            H_local[hot_temp_pos[0], hot_temp_pos[1]] += 1
            H_local[cold_temp_pos[0], cold_temp_pos[1]] -= 1
            # send to core 0
            # as long as size is somewhat small, this barrier won't slow things down much and ensures a correct k value
            comm.Barrier()
        comm.Barrier()
        comm.Reduce(H_local, H_master, op=MPI.SUM, root=0)
        # analysis
        if rank == 0 and (i > 0):
            # print np.count_nonzero(H_master)
            dt_dx, heat_flux, dt_dx_err, k, k_err, r2 = analysis.check_convergence_2d_onlat(H_master, tot_walkers,
                                                                                            grid.size, tot_time)
            # since final k is based on core 0 calculations, heat flux will slide a little since
            # core 0 will run slower, and this gives a more accurate result
            # np.savetxt("%s/H.txt" % plot_save_dir, H_master, fmt='%d')  # write histo to file
            k_list.append(k)
            dt_dx_list.append(dt_dx)
            heat_flux_list.append(heat_flux)
            timestep_list.append(core_time)
            logging.info("Parallel iteration %d out of %d, timestep %d, %d walkers, R2: %.4f, "
                         "k: %.4E, heat flux: %.4E, dT(x)/dx: %.4E"
                         % (i, walkers_per_core_whole, core_time, cur_num_walkers, r2, k, heat_flux, dt_dx))
        comm.Barrier()

    comm.Barrier()  # make sure whole walks are done

    if rank == 0:
        logging.info('Finished random walks, histogramming...')
        analysis.final_conductivity_onlat(plot_save_dir, prob_m_cn, dt_dx_list, k_list, k_conv_error_buffer)
        end = MPI.Wtime()
        logging.info("Constant flux simulation has completed")
        logging.info("Using %d cores, parallel simulation time was %.4f min" % (size, (end - start) / 60.0))
        walk_sec = tot_walkers / (end - start)
        logging.info("Crunched %.4f walkers/second" % walk_sec)
        temp_profile = plots.plot_colormap_2d(grid, H_master, quiet, plot_save_dir, gen_plots)
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
