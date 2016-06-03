import numpy as np
import argparse
import plots
import setup
import randomwalk


def yesno():
    response = raw_input('    Continue? (y/n) ')
    if len(response) == 0:  # [CR] returns true
        return True
    elif response[0] == 'n' or response[0] == 'N':
        return False
    else:  # Default
        return True

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='.')
    parser.add_argument('--grid_size', type=int, default=100, help='Size of square grid of use')
    parser.add_argument('--tube_length', type=float, default=20, help='Length of nanotube')
    parser.add_argument('--num_tubes', type=int, default=100, help='How many tubes are there for random walker to use')
    parser.add_argument('--mean_dist_tubes', type=int, default=40,
                        help='Mean distance between centers of tubes')
    parser.add_argument('--exact_sol', type=bool, default=False, help='Exact integer solutions used when calculating '
                                                                      'positions of nanotubes. Otherwise x and y values'
                                                                      ' rounded.')
    parser.add_argument('--orientation', type=str, default='random', help='Orientation of nanotubes in medium. random, '
                                                                          'horizontal, vertical, or angle in DEGREES')
    parser.add_argument('--walker_start_dir', type=str, default='left', help='Direction positive walkers come in from. '
                                                                         'top, bottom, left, right for 2D or '
                                                                         'Miller Index for 3D')
    parser.add_argument('--num_walkers', type=int, default=50, help='How many total walkers (pos and neg) are used')
    parser.add_argument('--timesteps', type=int, default=50, help='How many steps to run the simulation for')
    parser.add_argument('--bd_condition', type=str, default='exit', help='What happens when walkers hit the wall.'
                                                                               'reflect - isothermal,'
                                                                            'exit - constant flux')
    args = parser.parse_args()

    tube_length = args.tube_length
    grid_size = args.grid_size
    num_tubes = args.num_tubes
    exact_sol = args.exact_sol
    orientation = args.orientation
    num_walkers = args.num_walkers
    mean_dist_tubes = np.abs(args.mean_dist_tubes)
    walker_start_dir = args.walker_start_dir
    bd_condition = args.bd_condition
    timesteps = args.timesteps

    grid = setup.Grid2D(grid_size, tube_length, num_tubes, mean_dist_tubes, exact_sol, orientation)
    plots.plot_two_d_random_walk_setup(grid.tube_coords, grid_size)
    walkers_pos = setup.setupwalkers2d(num_walkers, grid_size, walker_start_dir) # list of walker object list
    #  at every timestep
    randomwalk.runrandomwalk2d(walkers_pos, grid_size, timesteps, bd_condition, walker_start_dir)

