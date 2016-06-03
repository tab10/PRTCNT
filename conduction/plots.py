import itertools
import matplotlib.pyplot as plt

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


def two_d_path_plot(x, y):
    plt.plot(x, y)
    plt.show()
