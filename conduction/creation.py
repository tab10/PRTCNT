import numpy as np
import logging
import os
import glob
from mpi4py import MPI


def check_for_folder(folder):
    if not os.path.exists(folder):
        os.mkdir(folder)


def get_plot_save_dir(folder, num_tubes, orientation, tube_length):
    ends = []
    for file in glob.glob("%d_%s_%d_*" % (num_tubes, orientation, tube_length)):
        name_temp = file.split('_')
        ends.append(name_temp[3])
    if not ends:  # empty list is false
        maxnum = 1
    else:
        maxnum = max(map(int, ends)) + 1
    plot_save_dir = "%d_%s_%d_%d" % (num_tubes, orientation, tube_length, maxnum)
    os.mkdir(plot_save_dir)
    return plot_save_dir


class Grid2D_onlat(object):
    def __init__(self, grid_size, tube_length, num_tubes, orientation):
        """Grid in first quadrant only for convenience"""
        logging.info("Setting up grid and tubes serially")
        self.size = grid_size
        if tube_length > grid_size:
            logging.error('Nanotube is too large for grid')
            raise SystemExit
        self.tube_coords = []
        self.tube_coords_l = []
        self.tube_coords_r = []
        self.tube_centers = []
        counter = 0  # counts num of non-unique tubes replaced
        if num_tubes > 0:  # tubes exist
            for i in range(num_tubes):  # currently no mean dist used, ADD LATER?
                x_l, y_l, x_r, y_r, x_c, y_c = self.generate_2d_tube(tube_length, orientation)
                self.tube_centers.append([x_c, y_c])
                self.tube_coords.append([x_l, y_l, x_r, y_r])
                self.tube_coords_l.append([x_l, y_l])
                self.tube_coords_r.append([x_r, y_r])
                if i >= 1:
                    uni_flag = self.check_tube_unique()  # ensures no endpoints, left or right, are in the same spot
                    while not uni_flag:
                        counter += 1
                        self.tube_centers.pop()
                        self.tube_coords.pop()
                        self.tube_coords_l.pop()
                        self.tube_coords_r.pop()
                        x_l, y_l, x_r, y_r, x_c, y_c = self.generate_2d_tube(tube_length, orientation)
                        self.tube_centers.append([x_c, y_c])
                        self.tube_coords.append([x_l, y_l, x_r, y_r])
                        self.tube_coords_l.append([x_l, y_l])
                        self.tube_coords_r.append([x_r, y_r])
                        uni_flag = self.check_tube_unique()
            logging.info("Corrected %d overlapping tube endpoints" % counter)
        self.tube_check_l, self.tube_check_r = self.generate_tube_check_array_2d()

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
        # print x_l, y_l, x_r, y_r
        # print self.euc_dist(x_l, y_l, x_r, y_r)
        return x_l, y_l, x_r, y_r, x_c, y_c

    def generate_tube_check_array_2d(self):
        tube_check_l = np.zeros((self.size + 1, self.size + 1), dtype=int)
        tube_check_r = np.zeros((self.size + 1, self.size + 1), dtype=int)
        for i in range(len(self.tube_coords)):
            tube_check_l[self.tube_coords[i][0], self.tube_coords[i][1]] = i
            tube_check_r[self.tube_coords[i][2], self.tube_coords[i][3]] = i
            # holds index of tube_coords, if a walker on that position has a nonzero value in this array,
            # pull the right or left tube endpoint (array positions are at left and right endpoints respectively)
        # np.set_printoptions(threshold=np.inf)
        # print tube_check_l
        return tube_check_l, tube_check_r

    def check_tube_unique(self):
        uni_flag = None
        current = self.tube_coords[-1]
        cur_l = [current[0], current[1]]
        cur_r = [current[2], current[3]]
        # separate
        tube_l = []
        tube_r = []
        for i in range(len(self.tube_coords) - 1):  # -1 accounts for not including the current tube
            tube_l.append([self.tube_coords[i][0], self.tube_coords[i][1]])
            tube_r.append([self.tube_coords[i][2], self.tube_coords[i][3]])
        if (cur_l in tube_l) or (cur_l in tube_r) or (cur_r in tube_l) or (cur_r in tube_r):
            uni_flag = False
        else:
            uni_flag = True
        # ensures no endpoints, left or right, are in the same spot
        return uni_flag

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
        self.tube_coords_l = []
        self.tube_coords_r = []
        self.tube_centers = []
        counter = 0  # counts num of non-unique tubes replaced
        if num_tubes > 0:  # tubes exist
            for i in range(num_tubes):  # currently no mean dist used, ADD LATER?
                x_l, y_l, z_l, x_r, y_r, z_r, x_c, y_c, z_c = self.generate_3d_tube(tube_length, tube_radius,
                                                                                    orientation)
                self.tube_centers.append([x_c, y_c, z_c])
                self.tube_coords.append([x_l, y_l, z_l, x_r, y_r, z_r])
                self.tube_coords_l.append([x_l, y_l, z_l])
                self.tube_coords_r.append([x_r, y_r, z_r])
                if i >= 1:
                    uni_flag = self.check_tube_unique()  # ensures no endpoints, left or right, are in the same spot
                    # print uni_flag
                    while not uni_flag:
                        counter += 1
                        self.tube_centers.pop()
                        self.tube_coords.pop()
                        self.tube_coords_l.pop()
                        self.tube_coords_r.pop()
                        x_l, y_l, z_l, x_r, y_r, z_r, x_c, y_c, z_c = self.generate_3d_tube(tube_length, tube_radius,
                                                                                            orientation)
                        self.tube_centers.append([x_c, y_c, z_c])
                        self.tube_coords.append([x_l, y_l, z_l, x_r, y_r, z_r])
                        self.tube_coords_l.append([x_l, y_l, z_l])
                        self.tube_coords_r.append([x_r, y_r, z_r])
                        uni_flag = self.check_tube_unique()
            logging.info("Corrected %d overlapping tube endpoints" % counter)
        self.tube_check_l, self.tube_check_r = self.generate_tube_check_array_3d()

    def generate_3d_tube(self, radius, tube_radius, orientation):
        """Finds appropriate angles within one degree that can be chosen from for random, should be good enough"""
        good_theta_angles = []
        good_phi_angles = []
        inside_box = False  # keeps track of if tube is in box
        while inside_box == False:
            # let's not put tube ends on the edges
            x_l = np.random.randint(1, self.size)
            y_l = np.random.randint(1, self.size)
            z_l = np.random.randint(1, self.size)
            if orientation == 'random':
                theta_angle_range = range(0, 360)  # y-z plane
                phi_angle_range = range(0, 360)  # x-y plane
                theta_angle = np.random.choice(theta_angle_range)
                phi_angle = np.random.choice(phi_angle_range)
            elif orientation == 'vertical':
                phi_angle = 0
                theta_angle_range = [0, 180]
                theta_angle = np.random.choice(theta_angle_range)
            elif orientation == 'horizontal':
                phi_angle = 0
                theta_angle_range = [90, 270]
                theta_angle = np.random.choice(theta_angle_range)
            else:
                logging.error("Invalid orientation specified")
                raise SystemExit
            x_r = round(self.coord(radius, theta_angle, phi_angle)[0] + x_l)
            y_r = round(self.coord(radius, theta_angle, phi_angle)[1] + y_l)
            z_r = round(self.coord(radius, theta_angle, phi_angle)[2] + z_l)
            if (x_r > 0) and (x_r < self.size) and (y_r > 0) and (y_r < self.size) and (z_r > 0) \
                    and (z_r < self.size):
                inside_box = True
        x_c = round(self.coord(radius, theta_angle, phi_angle)[0] + x_l) / 2
        y_c = round(self.coord(radius, theta_angle, phi_angle)[1] + y_l) / 2
        z_c = round(self.coord(radius, theta_angle, phi_angle)[2] + z_l) / 2
        return x_l, y_l, z_l, x_r, y_r, z_r, x_c, y_c, z_c

    def generate_tube_check_array_3d(self):
        tube_check_l = np.zeros((self.size + 1, self.size + 1, self.size + 1), dtype=int)
        tube_check_r = np.zeros((self.size + 1, self.size + 1, self.size + 1), dtype=int)
        for i in range(len(self.tube_coords)):
            tube_check_l[int(self.tube_coords[i][0]), int(self.tube_coords[i][1]), int(self.tube_coords[i][2])] = i
            tube_check_r[int(self.tube_coords[i][3]), int(self.tube_coords[i][4]), int(self.tube_coords[i][5])] = i
            # holds index of tube_coords, if a walker on that position has a nonzero value in this array,
            # pull the right or left tube endpoint (array positions are at left and right endpoints respectively)
        return tube_check_l, tube_check_r

    def check_tube_unique(self):
        uni_flag = None
        current = self.tube_coords[-1]
        cur_l = [current[0], current[1], current[2]]
        cur_r = [current[3], current[4], current[5]]
        # separate
        tube_l = []
        tube_r = []
        for i in range(len(self.tube_coords) - 1):  # -1 accounts for not including the current tube
            tube_l.append([self.tube_coords[i][0], self.tube_coords[i][1], self.tube_coords[i][2]])
            tube_r.append([self.tube_coords[i][3], self.tube_coords[i][4], self.tube_coords[i][5]])
        if (cur_l in tube_l) or (cur_l in tube_r) or (cur_r in tube_l) or (cur_r in tube_r):
            uni_flag = False
        else:
            uni_flag = True
        # ensures no endpoints, left or right, are in the same spot
        return uni_flag

    @staticmethod
    def coord(radius, theta_angle, phi_angle):
        # convention - theta from + z axis, phi from + x axis
        x = radius * np.sin(np.deg2rad(theta_angle)) * np.cos(np.deg2rad(phi_angle))
        y = radius * np.sin(np.deg2rad(theta_angle)) * np.sin(np.deg2rad(phi_angle))
        z = radius * np.cos(np.deg2rad(theta_angle))
        return x, y, z

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


class Walker3D_onlat(object):
    def __init__(self, grid_size, temp):
        if temp == 'hot':
            start_x = 0
        elif temp == 'cold':
            start_x = grid_size
        else:
            logging.error('Invalid walker temperature')
            raise SystemExit
        start_y = np.random.randint(0, grid_size + 1)
        start_z = np.random.randint(0, grid_size + 1)
        start = [start_x, start_y, start_z]
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
