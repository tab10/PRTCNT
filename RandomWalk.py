import numpy as np
import matplotlib.pyplot as plt
import argparse
import itertools


def setup_two_d_random_walk(grid_size, tube_length, num_tubes, mean_dist_tubes):
    if tube_length > grid_size:
        print 'Nanotube is too large for grid.'
        raise SystemExit
    tube_coords = []
    tube_centers = []
    if num_tubes > 0:  # tubes exist
        for i in range(num_tubes):  # currently no mean dist used, ADD LATER?
            tube_centers.append((np.random.randint(-grid_size + tube_length, grid_size - tube_length),
                                   np.random.randint(-grid_size + tube_length,
                                                     grid_size - tube_length)))  # reduced to ensure tubes stay in grid
            tube_coords.append(
                (generate_2d_tube(tube_centers[i][0], tube_centers[i][1], tube_length, grid_size)))
    return tube_coords, tube_centers


def two_d_random_walk(grid_size, tube_coords):
    ### Start the random walk
    step = 0
    moves = [(0, 1), (1, 0), (0, -1), (-1, 0)]
    x = []
    y = []
    x.append(x_i)
    y.append(y_i)
    while (x[-1] != x_f) or (y[-1] != y_f):
        step += 1
        dx, dy = moves[np.random.randint(0, 3)]
        x.append(x[-1] + dx)
        y.append(y[-1] + dy)
        # bd conditions
        if x[-1] > grid_size:
            x[-1] = -grid_size
        elif x[-1] < -grid_size:
            x[-1] = np.abs(grid_size)
        if y[-1] > grid_size:
            y[-1] = -grid_size
        elif y[-1] < -grid_size:
            y[-1] = np.abs(grid_size)
    print 'Done in %d steps' % step
    return x, y


def plot_two_d_random_walk_setup(tube_coords, grid_size):
    tube_x = []
    tube_y = []  # x and y coordinates of tubes unzipped
    for i in range(len(tube_coords)):
        tube_x.append(tube_coords[i][0::2])
        tube_y.append(tube_coords[i][1::2])
    colors = itertools.cycle(['b', 'g', 'r', 'c', 'm', 'y', 'k'])
    for i in range(len(tube_x)):
        plt.plot(tube_x[i], tube_y[i], c=next(colors), linestyle='--')
    plt.xlim(-grid_size, grid_size)
    plt.ylim(-grid_size, grid_size)
    plt.grid()
    plt.savefig('plots/setup.pdf')
    plt.show()


def int_on_circle(
        radius):  # finds all integer solutions on the circumference of a circle centered at origin for a given radius
    maxlegdist = int(np.floor(radius))
    sols = []
    for i in range(maxlegdist + 1):
        test = np.sqrt(radius ** 2 - i ** 2)
        sol_test = test.is_integer()
        if sol_test == True:
            sols.append((i, int(test)))  # stores all x,y integer solutions as tuples
            sols.append((i, -int(test)))
            sols.append((-i, int(test)))
            sols.append((-i, -int(test)))
    return sols


def generate_2d_tube(x_c, y_c, radius, grid_size,
                     int_sols=False):  # integer solutions finds exact solutions only for endpoints with given radius, otherwise rounded, giving more possibilities
    x_l = -grid_size - 1
    x_r = grid_size + 1
    y_l = -grid_size - 1
    y_r = grid_size + 1  # sets points outside grid so that while statement will continue until points inside grid
    while (x_l < -grid_size) or (x_r > grid_size) or (y_l < -grid_size) or (y_r > grid_size):
        if int_sols == True:
            sols = int_on_circle(radius)
            num_sols = len(sols)
            choice = np.random.randint(0, num_sols)  # random index of solution to use to make tube
            x_l = x_c - sols[choice][0]
            x_r = x_c + sols[choice][0]
            y_l = y_c - sols[choice][1]
            y_r = y_c + sols[choice][1]
        else:
            angle = 2.0 * np.pi * np.random.random_sample()  # randomly chosen
            x_vect = radius * np.cos(angle)
            y_vect = np.sqrt(radius ** 2 - x_vect ** 2)
            x_l = x_c - int(round(x_vect))
            x_r = x_c + int(round(x_vect))
            y_l = y_c - int(round(y_vect))
            y_r = y_c + int(round(y_vect))
    return x_l, y_l, x_r, y_r


def two_d_path_plot(x, y):
    plt.plot(x, y)
    plt.show()


def taxicab_dist(x0, y0, x1, y1):
    dist = np.abs(x1 - x0) + np.abs(y1 - y0)
    return dist


def radius(x, y):
    radius = np.sqrt(x ** 2 + y ** 2)
    return radius


def euc_dist(x0, y0, x1, y1):
    dist = np.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2)
    return dist


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
    parser.add_argument('--num_tubes', type=int, default=20, help='How many tubes are there for random walker to use')
    parser.add_argument('--mean_dist_tubes', type=int, default=40,
                        help='Mean distance between centers of tubes')
    args = parser.parse_args()

    tube_length = args.tube_length
    grid_size = args.grid_size
    num_tubes = args.num_tubes
    mean_dist_tubes = np.abs(args.mean_dist_tubes)

    tube_coords, tube_centers = setup_two_d_random_walk(grid_size, tube_length, num_tubes, mean_dist_tubes)
    plot_two_d_random_walk_setup(tube_coords, grid_size)
