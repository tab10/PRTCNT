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
                                                                          'horizontal, or vertical.')
    args = parser.parse_args()

    tube_length = args.tube_length
    grid_size = args.grid_size
    num_tubes = args.num_tubes
    exact_sol = args.exact_sol
    orientation = args.orientation
    mean_dist_tubes = np.abs(args.mean_dist_tubes)

    grid = setup.Grid2D(grid_size, tube_length, num_tubes, mean_dist_tubes)
    plots.plot_two_d_random_walk_setup(grid.tube_coords, grid_size)
