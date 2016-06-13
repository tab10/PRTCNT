import logging

import numpy as np

import analysis
import creation
import plots
import randomwalk
import run


def sim_3d_onlat(grid_size, tube_length, num_tubes, orientation, timesteps, tube_radius, save_loc_data,
                 quiet, save_loc_plots, save_dir, plot_save_dir):
    plot_save_dir = save_dir + "/plots"
    walker_data_save_dir = save_dir + "/walker_locations"
    walker_plot_save_dir = save_dir + "/walker_plots"

    logging.info("Setting up grid and tubes")
    grid = creation.Grid3D_onlat(grid_size, tube_length, tube_radius, num_tubes, orientation)
    plots.plot_three_d_random_walk_setup(grid.tube_coords, grid.size, quiet, plot_save_dir)

    raise SystemExit

    grid_range = [[0, grid.size], [0, grid.size]]
    bins = grid.size
    H = np.zeros((grid.size, grid.size))

    for i in range(num_walkers):
        # run hot walker
        logging.info("Start hot walker %d out of %d" % (i + 1, num_walkers))
        walker = randomwalk.runrandomwalk_2d_onlat(grid, timesteps, 'hot')
        if save_loc_data:
            run.save_walker_loc(walker, walker_data_save_dir, i, 'hot')
        if i == 0 & save_loc_plots == False:  # always save one example trajectory plot
            plots.plot_walker_path_2d_onlat(walker, grid_size, 'hot', quiet, i + 1, plot_save_dir)
        elif save_loc_plots:
            plots.plot_walker_path_2d_onlat(walker, grid_size, 'hot', quiet, i + 1, walker_plot_save_dir)
        H_temp, xedges, yedges = plots.histogram_walker_2d_onlat(walker, grid_range, 'hot', bins)
        H += H_temp

        # run cold walker
        logging.info("Start cold walker %d out of %d" % (i + 1, num_walkers))
        walker = randomwalk.runrandomwalk_2d_onlat(grid, timesteps, 'cold')
        if save_loc_data:
            run.save_walker_loc(walker, walker_data_save_dir, i, 'cold')
        if i == 0 & save_loc_plots == False:
            plots.plot_walker_path_2d_onlat(walker, grid_size, 'cold', quiet, i + 1, plot_save_dir)
        elif save_loc_plots:
            plots.plot_walker_path_2d_onlat(walker, grid_size, 'cold', quiet, i + 1, walker_plot_save_dir)
        H_temp, xedges, yedges = plots.histogram_walker_2d_onlat(walker, grid_range, 'cold', bins)
        H -= H_temp

    logging.info("Finished random walks")
    temp_profile = plots.plot_histogram_walkers_2d_onlat(timesteps, H, xedges, yedges, quiet, plot_save_dir)
    temp_gradient_x, gradient_cutoff = plots.plot_temp_gradient_2d_onlat(temp_profile, xedges, yedges, quiet,
                                                                         plot_save_dir)
    plots.plot_check_gradient_noise_floor(temp_gradient_x, quiet, plot_save_dir)
    analysis.calc_conductivity_2d_onlat(num_walkers, temp_gradient_x, grid.size, timesteps, gradient_cutoff)
    logging.info("Complete")
