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


"""test_2d.py - Tim Burt 3/24/2017
CONDUCTION package

This file is designed to test the random walk rules to make sure they satisfy detailed balance, or that the
system has microscopic reversibility and locations are equally probable as they must be in equilibrium"""

from __future__ import division
import logging
import numpy as np
import time
from mpi4py import MPI

from conduction import creation_2d
from conduction import plots
from conduction import analysis
from conduction import rules_2d


def parallel_method(grid_size, tube_length, tube_radius, num_tubes, orientation, tot_time, quiet, plot_save_dir,
                    gen_plots, kapitza, prob_m_cn, tot_walkers, disable_func, rank, size, rules_test, restart,
                    inert_vol):

    comm = MPI.COMM_WORLD

    # serial tube generation
    if rank == 0:
        grid = creation_2d.Grid2D_onlat(grid_size, tube_length, num_tubes, orientation, tube_radius, False,
                                        plot_save_dir,
                                        disable_func, rules_test, inert_vol)

    comm.Barrier()

    if rank == 0:
        if gen_plots:
            plots.plot_two_d_random_walk_setup(grid, quiet, plot_save_dir)
            plots.plot_colormap_2d(grid, grid.tube_check_bd_vol, quiet, plot_save_dir, gen_plots, title='CNT Type',
                                   filename='type', bds=False,
                                   vmin=-2, vmax=2)
            # plots.plot_colormap_2d(grid, grid.tube_check_index, quiet, plot_save_dir, gen_plots, title='CNT Index',
            #                       filename='index', bds=False)
            # plots.plot_colormap_2d(grid, grid.p_cn_m, quiet, plot_save_dir, gen_plots,
            #                       title='P_CN_M at each pixel',
            #                       xlab='X', ylab='Y', filename='P_cn_m')
            # plots.plot_colormap_2d(grid, grid.tube_bds_lkup, quiet, plot_save_dir, gen_plots,
            #                        title='CNT Boundaries',
            #                        xlab='X', ylab='Y', filename='tube_bds', bds=True)
            #plots.plot_check_array_2d(grid, quiet, plot_save_dir, gen_plots)
    else:
        grid = None

    comm.Barrier()
    grid = comm.bcast(grid, root=0)

    bins = grid.size + 1

    start = MPI.Wtime()

    # these will hold the walker objects, need to access walker.pos for the positions

    x_edges = np.asarray(range(0, bins)) + 0.25
    y_edges = np.asarray(range(0, bins)) + 0.25

    walkers_per_core_whole = int(np.floor(tot_walkers / size))
    walkers_per_core_remain = int(tot_walkers % size)  # run a last iteration on only this many cores

    comm.Barrier()

    H_local = np.zeros((grid.size + 1, grid.size + 1), dtype=float)
    if walkers_per_core_remain == 0:
        iterations = float(walkers_per_core_whole)
    else:
        iterations = float(walkers_per_core_whole + 1)

    for i in range(walkers_per_core_whole):
        cur_num_walkers = i * size
        H_master = np.zeros((grid.size + 1, grid.size + 1), dtype=float)  # should be reset every iteration
        walker_temp = rules_2d.runrandomwalk_2d_onlat(grid, tot_time, 'hot', kapitza, prob_m_cn, grid.bound,
                                                      rules_test)  # 'hot' since +1
        # histogram ALL positions
        for j in range(len(walker_temp.pos)):
            H_local[walker_temp.pos[j][0], walker_temp.pos[j][1]] += 1.0

        comm.Barrier()
        comm.Reduce(H_local, H_master, op=MPI.SUM, root=0)

        if rank == 0 and (i > 0):
            per_complete = float(i) * 100.0 / iterations
            """DEBUG - if temp_profile_norm_sum approaches 1, the histogram reductions are working correctly"""
            # temp_profile_norm = np.divide(H_master, (float(cur_num_walkers) * float(tot_time)))
            # temp_profile_norm_sum = np.sum(temp_profile_norm)
            # print(temp_profile_norm_sum)
            """END DEBUG"""
            logging.info("Parallel iteration %d out of %d, %d walkers, %d percent complete"
                         % (i, iterations, cur_num_walkers, per_complete))
        comm.Barrier()

    comm.Barrier()  # now remainder set

    if rank < walkers_per_core_remain:
        cur_num_walkers = walkers_per_core_whole * size + walkers_per_core_remain
        H_master = np.zeros((grid.size + 1, grid.size + 1), dtype=float)  # should be reset every iteration
        walker_temp = rules_2d.runrandomwalk_2d_onlat(grid, tot_time, 'hot', kapitza, prob_m_cn, grid.bound,
                                                      rules_test)  # 'hot' since +1
        # histogram ALL positions
        for j in range(len(walker_temp.pos)):
            H_local[walker_temp.pos[j][0], walker_temp.pos[j][1]] += 1.0

    comm.Barrier()
    comm.Reduce(H_local, H_master, op=MPI.SUM, root=0)

    # analysis
    if rank == 0:
        logging.info("Parallel iteration %d out of %d, %d walkers, 100 percent complete"
                     % (iterations, iterations, cur_num_walkers))

    comm.Barrier()  # make sure whole walks are done

    if rank == 0:
        logging.info('Finished random walk rules test. Analyzing....')
        mean_temp, mean_temp_norm, std_temp, std_temp_norm, temp_profile_norm = analysis.rules_test_analysis \
            (H_master, tot_walkers, tot_time)
        logging.info('Histogram: mean %.4E, std %.4E' % (mean_temp, std_temp))
        logging.info('Histogram normalized: mean %.4E, std %.4E' % (mean_temp_norm, std_temp_norm))
        # plots
        temp_profile = H_master
        plots.plot_colormap_2d(grid, temp_profile, quiet, plot_save_dir, gen_plots,
                               title='Number of times visited',
                               xlab='X', ylab='Y', filename='H_xy')
        plots.plot_colormap_2d(grid, temp_profile_norm, quiet, plot_save_dir, gen_plots,
                               title='Probability of being visited',
                               xlab='X', ylab='Y', filename='H_xy_norm')
        end = MPI.Wtime()
        logging.info("Rules test has completed. Please see results to verify if rules obey P.D.B.")
        logging.info("Using %d cores, parallel simulation time was %.4f min" % (size, (end - start) / 60.0))
        walk_sec = tot_walkers / (end - start)
        logging.info("Crunched %.4f walkers/second" % walk_sec)
        logging.info("Complete")
