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
        """Finds appropriate angles within one degree that can be chosen from for random, should be good enough.
        This method generates better tubes then a setup similar to the 3D method"""
        # first generate a left endpoint anywhere in the box
        x_l = np.random.randint(0, self.size + 1)
        y_l = np.random.randint(0, self.size + 1)
        good_angles = []
        if orientation == 'random':
            angle_range = range(0, 360)
        elif orientation == 'vertical':
            angle_range = [90, 270]
        elif orientation == 'horizontal':
            angle_range = [0, 180]
        else:
            logging.error("Invalid orientation specified")
            raise SystemExit
        for i in angle_range:  # round ensures endpoints stay on-grid
            x_test = round(radius * np.cos(np.deg2rad(i)) + x_l)
            y_test = round(radius * np.sin(np.deg2rad(i)) + y_l)
            if (x_test >= 0) & (x_test <= self.size) & (y_test >= 0) & (y_test <= self.size):
                # angle and x_l choice inside box
                good_angles.append(i)
        if good_angles == []:
            logging.error("Check box size and/or tube specs. No tubes can fit in the box.")
            raise SystemExit
        angle = np.random.choice(good_angles)  # in degrees
        x_r = round(radius * np.cos(np.deg2rad(angle)) + x_l)
        y_r = round(radius * np.sin(np.deg2rad(angle)) + y_l)
        if not ((x_l >= 0) & (x_r <= self.size) & (y_l >= 0) & (y_r <= self.size)):
            logging.error("Point found outside box.")
            raise SystemExit
        x_c = round(radius * np.cos(np.deg2rad(angle)) + x_l) / 2
        y_c = round(radius * np.sin(np.deg2rad(angle)) + y_l) / 2
        return x_l, y_l, x_r, y_r, x_c, y_c

    @staticmethod
    def taxicab_dist(x0, y0, x1, y1):
        dist = np.abs(x1 - x0) + np.abs(y1 - y0)
        return dist

    @staticmethod
    def radius(x, y):
        radius = np.sqrt(x ** 2 + y ** 2)
        return radius

    @staticmethod
    def euc_dist(x0, y0, x1, y1):
        dist = np.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2)
        return dist


class Grid3D_onlat(object):
    def __init__(self, grid_size, tube_length, tube_radius, num_tubes, orientation):
        """Grid in first quadrant only for convenience"""
        self.size = grid_size
        if tube_length > grid_size:
            logging.error('Nanotube is too large for grid')
            raise SystemExit
        self.tube_coords = []
        self.tube_centers = []
        if num_tubes > 0:  # tubes exist
            for i in range(num_tubes):  # currently no mean dist used, ADD LATER?
                x_l, y_l, z_l, x_r, y_r, z_r, x_c, y_c, z_c = self.generate_3d_tube(tube_length, tube_radius,
                                                                                    orientation)
                self.tube_centers.append([x_c, y_c, z_c])
                self.tube_coords.append([x_l, y_l, z_l, x_r, y_r, z_r])

    def generate_3d_tube(self, radius, tube_radius, orientation):
        """Finds appropriate angles within one degree that can be chosen from for random, should be good enough"""
        good_theta_angles = []
        good_phi_angles = []
        if orientation == 'random':
            x_l = np.random.randint(radius, self.size + 1 - radius)
            y_l = np.random.randint(radius, self.size + 1 - radius)
            z_l = np.random.randint(radius, self.size + 1 - radius)
            theta_angle_range = range(0, 360)  # y-z plane
            phi_angle_range = range(0, 360)  # x-y plane
            theta_angle = np.random.choice(theta_angle_range)
            phi_angle = np.random.choice(phi_angle_range)
        elif orientation == 'vertical':
            x_l = np.random.randint(0, self.size + 1)
            y_l = np.random.randint(0, self.size + 1)
            z_l = np.random.randint(radius, self.size + 1 - radius)
            phi_angle = 0
            theta_angle_range = [180]
            theta_angle = np.random.choice(theta_angle_range)
        elif orientation == 'horizontal':
            x_l = np.random.randint(radius, self.size + 1 - radius)
            y_l = np.random.randint(0, self.size + 1)
            z_l = np.random.randint(0, self.size + 1)
            phi_angle = 0
            theta_angle_range = [90, 270]
            theta_angle = np.random.choice(theta_angle_range)
        else:
            logging.error("Invalid orientation specified")
            raise SystemExit
        x_r = round(radius * np.sin(np.deg2rad(theta_angle)) * np.cos(np.deg2rad(phi_angle)) + x_l)
        y_r = round(radius * np.sin(np.deg2rad(theta_angle)) * np.sin(np.deg2rad(phi_angle)) + y_l)
        z_r = round(radius * np.cos(np.deg2rad(phi_angle)) + z_l)
        if not ((x_l >= 0) & (x_r <= self.size) & (y_l >= 0) & (y_r <= self.size) &
                    (z_l >= 0) & (z_r <= self.size)):
            logging.error("Point found outside box.")
            raise SystemExit
        x_c = round(radius * np.sin(np.deg2rad(theta_angle)) * np.cos(np.deg2rad(phi_angle)) + x_l) / 2
        y_c = round(radius * np.sin(np.deg2rad(theta_angle)) * np.sin(np.deg2rad(phi_angle)) + y_l) / 2
        z_c = round(radius * np.cos(np.deg2rad(phi_angle)) + z_l) / 2
        return x_l, y_l, z_l, x_r, y_r, z_r, x_c, y_c, z_c

    @staticmethod
    def taxicab_dist(x0, y0, z0, x1, y1, z1):
        dist = np.abs(x1 - x0) + np.abs(y1 - y0) + np.abs(z1 - z0)
        return dist

    @staticmethod
    def radius(x, y, z):
        radius = np.sqrt(x ** 2 + y ** 2 + z ** 2)
        return radius

    @staticmethod
    def euc_dist(x0, y0, z0, x1, y1, z1):
        dist = np.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2 + (z1 - z0) ** 2)
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
        start_y = np.random.randint(0, grid_size + 1)
        start = [start_x, start_y]
        self.pos = [start]

    def add_dpos(self, newpos):  # add incremented position
        arr1 = np.array(self.pos[-1])
        arr2 = np.array(newpos)
        # print arr1, arr2
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
