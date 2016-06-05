import numpy as np
import logging
import os


def check_for_folder(folder):
    if not os.path.exists(folder):
        os.mkdir(folder)


class Grid2D_onlat(object):
    def __init__(self, grid_size, tube_length, num_tubes, orientation):
        """Grid in first quadrant only for convenience"""
        self.size = grid_size
        if tube_length > grid_size:
            logging.error('Nanotube is too large for grid')
            raise SystemExit
        self.tube_coords = []
        self.tube_centers = []
        if num_tubes > 0:  # tubes exist
            for i in range(num_tubes):  # currently no mean dist used, ADD LATER?
                x_l, y_l, x_r, y_r, x_c, y_c = self.generate_2d_tube(tube_length, orientation)
                self.tube_centers.append([x_c, y_c])
                self.tube_coords.append([x_l, y_l, x_r, y_r])

    def generate_2d_tube(self, radius, orientation):
        """Finds appropriate angles within one degree that can be chosen from for random, should be good enough"""
        # first generate a left endpoint anywhere in the box
        x_l = np.random.randint(0, self.size + 1)
        y_l = np.random.randint(0, self.size + 1)
        good_angles = []
        if orientation == 'random':
            angle_range = range(0,360)
        elif orientation == 'vertical':
            angle_range = [90, 270]
        elif orientation == 'horizontal':
            angle_range = [0, 180]
        for i in angle_range:
            x_test = radius * np.cos(np.deg2rad(i)) + x_l
            y_test = radius * np.sin(np.deg2rad(i)) + y_l
            if (x_test >= 0) & (x_test <= self.size) & (y_test >= 0) & (y_test <= self.size):
                # angle and x_l choice inside box
                good_angles.append(i)
        angle = np.random.choice(good_angles)  # in degrees
        x_r = radius * np.cos(np.deg2rad(angle)) + x_l
        y_r = radius * np.sin(np.deg2rad(angle)) + y_l
        x_c = (radius/2.0) * np.cos(np.deg2rad(angle)) + x_l
        y_c = (radius/2.0) * np.sin(np.deg2rad(angle)) + y_l
        return x_l, y_l, x_r, y_r, x_c, y_c


def taxicab_dist(x0, y0, x1, y1):
    dist = np.abs(x1 - x0) + np.abs(y1 - y0)
    return dist


def radius(x, y):
    radius = np.sqrt(x ** 2 + y ** 2)
    return radius


def euc_dist(x0, y0, x1, y1):
    dist = np.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2)
    return dist


class Walker2D_onlat(object):
    def __init__(self, grid_size, temp):
        if temp == 'hot':
            start_x = 0
        elif temp == 'cold':
            start_x = grid_size
        else:
            logging.error('Invalid walker temperature')
            raise SystemExit
        start_y = np.random.randint(0, grid_size+1)
        start = [start_x, start_y]
        self.pos = [start]

    def add_dpos(self, newpos):  # add incremented position
        arr1 = np.array(self.pos[-1])
        arr2 = np.array(newpos)
        #print arr1, arr2
        self.pos.append(list(arr1 + arr2))

    def add_pos(self, newpos):  # add new position
        self.pos.append(list(newpos))

    def replace_pos(self, newpos):  # replace current position
        self.pos[-1] = list(newpos)

    def replace_dpos(self, newpos):  # replace incremented position
        arr1 = np.array(self.pos[-1])
        arr2 = np.array(newpos)
        self.pos[-1] = list(arr1 + arr2)


def int_on_circle(radius):  # finds all integer solutions on the circumference of a circle
    # centered at origin for a given radius
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