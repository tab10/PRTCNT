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
import numpy as np
import logging
from conduction import *


class Grid2D_onlat(object):
    def __init__(self, grid_size, tube_length, num_tubes, orientation, tube_radius, parallel, plot_save_dir,
                 disable_func, rank=None,
                 size=None):
        """Grid in first quadrant only for convenience"""
        # serial implementation
        logging.info("Setting up grid and tubes serially")
        self.size = grid_size
        if tube_length > grid_size:
            logging.error('Nanotube is too large for grid')
            raise SystemExit
        self.tube_coords = []
        self.tube_coords_l = []
        self.tube_coords_r = []
        self.tube_centers = []
        self.theta = []
        self.tube_radius = tube_radius
        self.setup_tube_vol_check_array_2d()
        counter = 0  # counts num of non-unique tubes replaced
        status_counter = 0
        if tube_radius == 0:
            logging.info("Zero tube radius given. Tubes will have no volume.")
            if disable_func:
                logging.info("Ignoring disabling functionalization since tubes are volumeless.")
            fill_fract = 2.0 * float(num_tubes) / grid_size ** 2
            logging.info("Filling fraction is %.2f %%" % (fill_fract * 100.0))
            backend.save_fill_frac(plot_save_dir, fill_fract)
            if num_tubes > 0:  # tubes exist
                for i in range(num_tubes):  # currently no mean dist used, ADD LATER?
                    if (i % 50) == 0:
                        status_counter += 50
                        logging.info('Generating tube %d...' % (status_counter - 50))
                    x_l, y_l, x_r, y_r, x_c, y_c, theta = self.generate_2d_tube(tube_length, orientation,
                                                                                tube_radius)
                    self.tube_centers.append([x_c, y_c])
                    self.tube_coords.append([x_l, y_l, x_r, y_r])
                    self.tube_coords_l.append([x_l, y_l])
                    self.tube_coords_r.append([x_r, y_r])
                    self.theta.append(theta)
                    if i == 0:
                        self.add_tube_vol_check_array_2d([x_l, y_l, x_r, y_r], None, disable_func)
                    if i >= 1:
                        uni_flag = self.check_tube_unique_2d_arraymethod([x_l, y_l, x_r, y_r])
                        # uni_flag = self.check_tube_unique(self.tube_coords, False)
                        #  ensures no endpoints, left or right, are in the same spot
                        while not uni_flag:
                            counter += 1
                            self.tube_centers.pop()
                            self.tube_coords.pop()
                            self.tube_coords_l.pop()
                            self.tube_coords_r.pop()
                            self.theta.pop()
                            x_l, y_l, x_r, y_r, x_c, y_c, theta = self.generate_2d_tube(tube_length, orientation,
                                                                                        tube_radius)
                            self.tube_centers.append([x_c, y_c])
                            self.tube_coords.append([x_l, y_l, x_r, y_r])
                            self.tube_coords_l.append([x_l, y_l])
                            self.tube_coords_r.append([x_r, y_r])
                            self.theta.append(theta)
                            uni_flag = self.check_tube_unique_2d_arraymethod([x_l, y_l, x_r, y_r])
                            # uni_flag = self.check_tube_unique(self.tube_coords, False)
                        self.add_tube_vol_check_array_2d([x_l, y_l, x_r, y_r], None, disable_func)
                logging.info("Corrected %d overlapping tube endpoints" % counter)
            self.tube_check_l, self.tube_check_r, self.tube_check_bd = self.generate_tube_check_array_2d()
        else:
            logging.info("Non-zero tube radius given. Tubes will have excluded volume.")
            l_d = tube_length / (2.0 * tube_radius)
            logging.info("L/D is %.4f." % l_d)
            self.tube_squares = []  # grid squares that a tube passes through, for every tube
            self.setup_tube_vol_check_array_2d()
            if num_tubes > 0:  # tubes exist
                for i in range(num_tubes):  # currently no mean dist used, ADD LATER?
                    if (i % 50) == 0:
                        status_counter += 50
                        logging.info('Generating tube %d...' % (status_counter - 50))
                    x_l, y_l, x_r, y_r, x_c, y_c, theta = self.generate_2d_tube(tube_length, orientation,
                                                                                tube_radius)
                    tube_squares = self.find_squares([x_l, y_l], [x_r, y_r], tube_radius)
                    self.tube_centers.append([x_c, y_c])
                    self.tube_coords.append([x_l, y_l, x_r, y_r])
                    self.tube_coords_l.append([x_l, y_l])
                    self.tube_coords_r.append([x_r, y_r])
                    self.tube_squares.append(tube_squares)
                    self.theta.append(theta)
                    if i == 0:
                        self.add_tube_vol_check_array_2d([x_l, y_l, x_r, y_r], tube_squares, disable_func)
                    if i >= 1:
                        uni_flag = self.check_tube_and_vol_unique_2d_arraymethod(tube_squares)
                        # uni_flag = self.check_tube_and_vol_unique(self.tube_squares, False)
                        while not uni_flag:
                            counter += 1
                            self.tube_centers.pop()
                            self.tube_coords.pop()
                            self.tube_coords_l.pop()
                            self.tube_coords_r.pop()
                            self.tube_squares.pop()
                            self.theta.pop()
                            x_l, y_l, x_r, y_r, x_c, y_c, theta = self.generate_2d_tube(tube_length, orientation,
                                                                                        tube_radius)
                            tube_squares = self.find_squares([x_l, y_l], [x_r, y_r], tube_radius)
                            self.theta.append(theta)
                            self.tube_centers.append([x_c, y_c])
                            self.tube_coords.append([x_l, y_l, x_r, y_r])
                            self.tube_coords_l.append([x_l, y_l])
                            self.tube_coords_r.append([x_r, y_r])
                            self.tube_squares.append(tube_squares)
                            # uni_flag = self.check_tube_and_vol_unique(self.tube_squares, False)
                            uni_flag = self.check_tube_and_vol_unique_2d_arraymethod(tube_squares)
                        # uni_flag must have been true to get here, so write it to the array
                        self.add_tube_vol_check_array_2d([x_l, y_l, x_r, y_r], tube_squares, disable_func)
                logging.info("Corrected %d overlapping tube endpoints and/or volume points" % counter)
            # get number of squares filled
            cube_count = 0  # each cube has area 1
            for i in range(len(self.tube_squares)):
                cube_count += len(self.tube_squares[i])
            fill_fract = float(cube_count) * 2.0 * tube_radius / grid_size ** 2
            # each cube has area 1, times the tube radius (important if not 1)
            logging.info("Filling fraction is %.2f %%" % (fill_fract * 100.0))
            backend.save_fill_frac(plot_save_dir, fill_fract)
            self.tube_check_l, self.tube_check_r, self.tube_check_bd = self.generate_tube_check_array_2d()
            self.tube_check_bd_vol, self.tube_check_index = self.generate_vol_check_array_2d(disable_func)
        self.avg_tube_len, self.std_tube_len, self.tube_lengths = self.check_tube_lengths()
        logging.info("Actual tube length avg+std: %.4f +- %.4f" % (self.avg_tube_len, self.std_tube_len))


    def generate_2d_tube(self, radius, orientation, tube_radius):
        """Finds appropriate angles within one degree that can be chosen from for random, should be good enough.
        This method generates better tubes then a setup similar to the 3D method"""
        # first generate a left endpoint anywhere in the box
        if tube_radius == 0:
            x_l = np.random.randint(1, self.size)  # nothing on boundaries
            y_l = np.random.randint(1, self.size)
        else:  # ensures no tube parts are outside grid
            x_l = np.random.randint(np.floor(tube_radius) + 1, self.size - np.floor(tube_radius))
            y_l = np.random.randint(np.floor(tube_radius) + 1, self.size - np.floor(tube_radius))
        good_angles = []
        if orientation == 'random':
            angle_range = list(range(0, 360))
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
            if (x_test > 0) & (x_test < self.size) & (y_test > 0) & (y_test < self.size):
                # angle and x_r choice inside box
                good_angles.append(i)
        if not good_angles:
            logging.error("Check box size and/or tube specs. No tubes can fit in the box.")
            raise SystemExit
        angle = np.random.choice(good_angles)  # in degrees
        x_r = int(round(radius * np.cos(np.deg2rad(angle)) + x_l))
        y_r = int(round(radius * np.sin(np.deg2rad(angle)) + y_l))
        if not ((x_l > 0) & (x_r < self.size + 1) & (y_l > 0) & (y_r < self.size + 1)):
            logging.error("Point found outside box.")
            raise SystemExit
        if x_l > x_r:  # this imposes left to right tube order WRT + x-axis
            x_l_temp = x_l
            x_r_temp = x_r
            y_l_temp = y_l
            y_r_temp = y_r
            x_l = x_r_temp
            x_r = x_l_temp
            y_l = y_r_temp
            y_r = y_l_temp
        x_c = round(radius * np.cos(np.deg2rad(angle)) + x_l) / 2
        y_c = round(radius * np.sin(np.deg2rad(angle)) + y_l) / 2
        # print x_l, y_l, x_r, y_r
        # print self.euc_dist(x_l, y_l, x_r, y_r)
        return x_l, y_l, x_r, y_r, x_c, y_c, angle

    def find_squares(self, start, end, tube_radius):
        """Bresenham's Line Algorithm
        Produces a list of tuples (bottom left corners of grid)
        All squares a tube passes through
        Also returns original start and end points in first and last positions respectively
        """
        # Setup initial conditions
        x1, y1 = start
        x2, y2 = end
        dx = x2 - x1
        dy = y2 - y1

        # Determine how steep the line is
        is_steep = abs(dy) > abs(dx)

        # Rotate line
        if is_steep:
            x1, y1 = y1, x1
            x2, y2 = y2, x2

        # Swap start and end points if necessary and store swap state
        swapped = False
        if x1 > x2:
            x1, x2 = x2, x1
            y1, y2 = y2, y1
            swapped = True

        # Recalculate differentials
        dx = x2 - x1
        dy = y2 - y1

        # Calculate error
        error = int(dx / 2.0)
        ystep = 1 if y1 < y2 else -1

        # Iterate over bounding box generating points between start and end
        y = y1
        points = []
        for x in range(x1, x2 + 1):
            coord = [y, x] if is_steep else [x, y]
            points.append(coord)
            error -= abs(dy)
            if error < 0:
                y += ystep
                error += dx

        # Reverse the list if the coordinates were swapped
        if swapped:
            points.reverse()

        # now check if tube radius > 1, if so add more volume points
        if tube_radius > 0.5:
            logging.info('Tube radius will be implemented here later if needed.')
        return points

    def generate_bd_and_vol(self, tube_squares, angle):
        """2D, This takes the bottom left squares from find_squares and creates volume and boundaries based
        on the tube radius, for one tube
        This function is not currently used for anything because boundaries are not implemented in this way
        """
        top = []
        bot = []
        left = []
        right = []
        l_end = tube_squares[0]
        r_end = tube_squares[-1]
        occupied_cubes = tube_squares  # this will be checked to ensure volume is exclusive too
        for i in range(1, self.tube_radius + 1):
            l_x_above = int(round(i * np.cos(np.deg2rad(angle + 90)) + l_end[0]))
            l_y_above = int(round(i * np.sin(np.deg2rad(angle + 90)) + l_end[1]))
            r_x_above = int(round(i * np.cos(np.deg2rad(angle + 90)) + r_end[0]))
            r_y_above = int(round(i * np.sin(np.deg2rad(angle + 90)) + r_end[1]))
            l_x_below = int(round(i * np.cos(np.deg2rad(angle - 90)) + l_end[0]))
            l_y_below = int(round(i * np.sin(np.deg2rad(angle - 90)) + l_end[1]))
            r_x_below = int(round(i * np.cos(np.deg2rad(angle - 90)) + r_end[0]))
            r_y_below = int(round(i * np.sin(np.deg2rad(angle - 90)) + r_end[1]))
            left.append([l_x_above, l_y_above])
            left.append([l_x_below, l_y_below])
            right.append([r_x_above, r_y_above])
            right.append([r_x_below, r_y_below])
        for i in range(len(tube_squares)):
            t_x = int(round(self.tube_radius * np.cos(np.deg2rad(angle + 90)) + tube_squares[i][0]))
            t_y = int(round(self.tube_radius * np.sin(np.deg2rad(angle + 90)) + tube_squares[i][1]))
            b_x = int(round(self.tube_radius * np.cos(np.deg2rad(angle - 90)) + tube_squares[i][0]))
            b_y = int(round(self.tube_radius * np.sin(np.deg2rad(angle - 90)) + tube_squares[i][1]))
            top.append([t_x, t_y])
            bot.append([b_x, b_y])
        total = top + bot + left + right
        return top, bot, left, right, total, occupied_cubes

    def generate_tube_check_array_2d(self):
        """To be used with no tube volume
        Generates a left and right lookup array that holds the index of the opposite endpoint"""
        tube_check_l = np.zeros((self.size + 1, self.size + 1), dtype=int)
        tube_check_r = np.zeros((self.size + 1, self.size + 1), dtype=int)
        bd = np.zeros((self.size + 1, self.size + 1), dtype=int)
        for i in range(0, len(self.tube_coords)):
            tube_check_l[self.tube_coords[i][0], self.tube_coords[i][1]] = i + 1  # THESE ARE OFFSET BY ONE
            tube_check_r[self.tube_coords[i][2], self.tube_coords[i][3]] = i + 1
            bd[self.tube_coords[i][0], self.tube_coords[i][1]] = 1  # endpoint
            bd[self.tube_coords[i][2], self.tube_coords[i][3]] = 1
            # holds index of tube_coords, if a walker on that position has a nonzero value in this array,
            # pull the right or left tube endpoint (array positions are at left and right endpoints respectively)
        # add boundary tags
        for i in range(self.size + 1):
            bd[0, i] = 10  # x = 0 is reflective
            bd[self.size, i] = 10  # x = grid.size is reflective
            bd[i, 0] = 20  # y = 0 is periodic
            bd[i, self.size] = 20  # y = grid.size is periodic
        # corners are special
        bd[0, 0] = 30
        bd[0, self.size] = 30
        bd[self.size, 0] = 30
        bd[self.size, self.size] = 30
        # np.set_printoptions(threshold=np.inf)
        # print tube_check_l
        return tube_check_l, tube_check_r, bd

    def generate_vol_check_array_2d(self, disable_func):
        """To be used with tube volume
        Generates a boundary/volume lookup array (0 nothing, 1 boundary, -1 volume)"""
        bd_vol = np.zeros((self.size + 1, self.size + 1), dtype=int)
        index = np.zeros((self.size + 1, self.size + 1), dtype=int)
        if disable_func:
            endpoint_val = -1  # treat endpoints as volume, changing the rules in the walk
        else:
            endpoint_val = 1  # leave it as endpoint
        for i in range(len(self.tube_coords)):
            bd_vol[self.tube_coords[i][0], self.tube_coords[i][1]] = endpoint_val  # left endpoints
            bd_vol[self.tube_coords[i][2], self.tube_coords[i][3]] = endpoint_val  # right endpoints
            index[self.tube_coords[i][0], self.tube_coords[i][1]] = i + 1  # THESE ARE OFFSET BY ONE
            index[self.tube_coords[i][2], self.tube_coords[i][3]] = i + 1  # THESE ARE OFFSET BY ONE
            for j in range(1, len(self.tube_squares[i]) - 1):
                bd_vol[self.tube_squares[i][j][0], self.tube_squares[i][j][1]] = -1  # volume points
                index[self.tube_squares[i][j][0], self.tube_squares[i][j][1]] = i + 1  # THESE ARE OFFSET BY ONE
        # add boundary tags
        for i in range(self.size + 1):
            bd_vol[0, i] = 10  # x = 0 is reflective
            bd_vol[self.size, i] = 10  # x = grid.size is reflective
            bd_vol[i, 0] = 20  # y = 0 is periodic
            bd_vol[i, self.size] = 20  # y = grid.size is periodic
        # corners are special
        bd_vol[0, 0] = 30
        bd_vol[0, self.size] = 30
        bd_vol[self.size, 0] = 30
        bd_vol[self.size, self.size] = 30
        # np.set_printoptions(threshold=np.inf)
        # print bd_vol
        return bd_vol, index

    def setup_tube_vol_check_array_2d(self):
        "Setup of tube check and index arrays, returns nothing"
        bd_vol = np.zeros((self.size + 1, self.size + 1), dtype=int)
        index = np.zeros((self.size + 1, self.size + 1), dtype=int)
        # add boundary tags
        for i in range(self.size + 1):
            bd_vol[0, i] = 10  # x = 0 is reflective
            bd_vol[self.size, i] = 10  # x = grid.size is reflective
            bd_vol[i, 0] = 20  # y = 0 is periodic
            bd_vol[i, self.size] = 20  # y = grid.size is periodic
        # corners are special
        bd_vol[0, 0] = 30
        bd_vol[0, self.size] = 30
        bd_vol[self.size, 0] = 30
        bd_vol[self.size, self.size] = 30
        self.tube_check_bd_vol = bd_vol
        self.tube_check_index = index

    def add_tube_vol_check_array_2d(self, new_tube_coords, new_tube_squares, disable_func):
        "Adds tube to the current check arrays"
        index_val = len(self.tube_coords) + 1  # THESE ARE OFFSET BY ONE
        if disable_func:
            endpoint_val = -1  # treat endpoints as volume, changing the rules in the walk
        else:
            endpoint_val = 1  # leave it as endpoint
        self.tube_check_bd_vol[new_tube_coords[0], new_tube_coords[1]] = endpoint_val  # left endpoints
        self.tube_check_bd_vol[new_tube_coords[2], new_tube_coords[3]] = endpoint_val  # right endpoints
        self.tube_check_index[new_tube_coords[0], new_tube_coords[1]] = index_val
        self.tube_check_index[new_tube_coords[2], new_tube_coords[3]] = index_val
        if new_tube_squares is not None:  # None used for tunneling only
            for j in range(1, len(new_tube_squares) - 1):
                self.tube_check_bd_vol[new_tube_squares[j][0], new_tube_squares[j][1]] = -1  # volume points
                self.tube_check_index[new_tube_squares[j][0], new_tube_squares[j][1]] = index_val

    def check_tube_unique_2d_arraymethod(self, new_tube_squares):
        "No volume"
        uni_flag = True
        check_l = self.tube_check_bd_vol[new_tube_squares[0], new_tube_squares[1]]
        check_r = self.tube_check_bd_vol[new_tube_squares[2], new_tube_squares[3]]
        if (check_l != 0) or (check_r != 0):
            uni_flag = False
        return uni_flag

    def check_tube_and_vol_unique_2d_arraymethod(self, new_tube_squares):
        "Volume"
        index_val = len(self.tube_coords) + 1  # current tube index
        uni_flag = True
        for l in range(len(new_tube_squares)):
            test_vol = self.tube_check_bd_vol[new_tube_squares[l][0], new_tube_squares[l][1]]
            if (test_vol == 1) or (test_vol == -1):  # new tube overlaps old volume or endpoint
                uni_flag = False
            # generate candidate check positions for tube crossing
            # this algorithm looks for interweaved diagonal clusters
            # check cur(tl,br) and old(tr,bl), pts wrt tl
            b_r = [new_tube_squares[l][0] + 1, new_tube_squares[l][1] - 1]
            old_t_r = self.tube_check_bd_vol[new_tube_squares[l][0] + 1, new_tube_squares[l][1]]
            old_b_l = self.tube_check_bd_vol[new_tube_squares[l][0], new_tube_squares[l][1] - 1]
            old_t_r_index = self.tube_check_index[new_tube_squares[l][0] + 1, new_tube_squares[l][1]]
            old_b_l_index = self.tube_check_index[new_tube_squares[l][0], new_tube_squares[l][1] - 1]
            if ((old_t_r == 1) or (old_t_r == -1)) and ((old_b_l == 1) or (old_b_l == -1)) \
                    and (b_r in new_tube_squares) and (old_t_r_index == old_b_l_index):  # we have a crossing
                uni_flag = False
            # check old(tl,br) and cur(tr,bl), pts wrt tr
            b_l = [new_tube_squares[l][0] - 1, new_tube_squares[l][1] - 1]
            old_t_l = self.tube_check_bd_vol[new_tube_squares[l][0] - 1, new_tube_squares[l][1]]
            old_b_r = self.tube_check_bd_vol[new_tube_squares[l][0], new_tube_squares[l][1] - 1]
            old_t_l_index = self.tube_check_index[new_tube_squares[l][0] - 1, new_tube_squares[l][1]]
            old_b_r_index = self.tube_check_index[new_tube_squares[l][0], new_tube_squares[l][1] - 1]
            if ((old_t_l == 1) or (old_t_l == -1)) and ((old_b_r == 1) or (old_b_r == -1)) \
                    and (b_l in new_tube_squares) and (old_t_l_index == old_b_r_index):  # we have a crossing
                uni_flag = False
        return uni_flag

    def check_tube_lengths(self):
        tube_lengths = np.zeros(len(self.tube_coords))
        for i in range(len(self.tube_coords)):
            dist = self.euc_dist(self.tube_coords[i][0], self.tube_coords[i][1], self.tube_coords[i][2],
                                 self.tube_coords[i][3])
            tube_lengths[i] = dist
        avg_tube_len = np.mean(tube_lengths)
        std_tube_len = np.std(tube_lengths, ddof=1)
        return avg_tube_len, std_tube_len, tube_lengths

    @staticmethod
    def check_tube_unique(coords_list, parallel, rank=None, size=None):
        uni_flag = None
        if not parallel:  # serial case
            temp = 1
        elif parallel:  # parallel case. numbers have true truth value.
            temp = size - rank  # SINCE INDEXING IS REVERSED
        current = coords_list[-temp]  # self.tube_coords[-1]
        cur_l = [current[0], current[1]]
        cur_r = [current[2], current[3]]
        # separate
        tube_l = []
        tube_r = []
        if not parallel:
            for i in range(len(coords_list) - 1):
                tube_l.append([coords_list[i][0], coords_list[i][1]])
                tube_r.append([coords_list[i][2], coords_list[i][3]])
        else:
            exclude_val = len(coords_list) - size + rank
            for i in range(len(coords_list)):
                if i != exclude_val:
                    tube_l.append([coords_list[i][0], coords_list[i][1]])
                    tube_r.append([coords_list[i][2], coords_list[i][3]])
        if (cur_l in tube_l) or (cur_l in tube_r) or (cur_r in tube_l) or (cur_r in tube_r):
            uni_flag = False
        else:
            uni_flag = True
        # ensures no endpoints, left or right, are in the same spot
        return uni_flag

    @staticmethod
    def check_tube_and_vol_unique(tube_squares, parallel, rank=None, size=None):
        """Checks all endpoints, boundaries, and volume for uniqueness and no overlap
        tube_squares holds the endpoints too, so just check the current one for uniqueness"""
        uni_flag = True
        if not parallel:  # serial case
            temp = 1
        elif parallel:  # parallel case. numbers have true truth value.
            temp = size - rank  # SINCE INDEXING IS REVERSED
        current_vol = tube_squares[-temp]  # self.tube_squares[-1], this holds all the points for one tube
        vol = []
        if not parallel:
            for i in range(len(tube_squares) - 1):
                for k in range(len(tube_squares[i])):
                    vol.append(tube_squares[i][k])
        else:
            exclude_val = len(tube_squares) - size + rank
            for i in range(len(tube_squares)):
                if i != exclude_val:
                    for k in range(len(tube_squares[i])):
                        vol.append(tube_squares[i][k])
        for l in range(len(current_vol)):
            if current_vol[l] in vol:
                uni_flag = False
            # generate candidate check positions for tube crossing
            ###########
            # this algorithm looks for interweaved diagonal clusters. The only 2 possibilities are checked
            t_l = current_vol[l]
            t_r = [current_vol[l][0] + 1, current_vol[l][1]]
            b_l = [current_vol[l][0], current_vol[l][1] - 1]
            b_r = [current_vol[l][0] + 1, current_vol[l][1] - 1]
            if (b_r in current_vol) and (t_r in vol) and (b_l in vol):
                uni_flag = False
            t_r = current_vol[l]
            t_l = [current_vol[l][0] - 1, current_vol[l][1]]
            b_l = [current_vol[l][0] - 1, current_vol[l][1] - 1]
            b_r = [current_vol[l][0], current_vol[l][1] - 1]
            if (t_l in vol) and (b_l in current_vol) and (b_r in vol):
                uni_flag = False
        return uni_flag

    @staticmethod
    def taxicab_dist(x0, y0, x1, y1):
        dist = np.abs(x1 - x0) + np.abs(y1 - y0)
        return dist

    @staticmethod
    def calc_radius(x, y):
        radius = np.sqrt(x ** 2 + y ** 2)
        return radius

    @staticmethod
    def euc_dist(x0, y0, x1, y1):
        dist = np.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2)
        return dist


class Walker2D_onlat(object):
    def __init__(self, grid_size, temp, rules_test):
        if not rules_test:
            if temp == 'hot':
                start_x = 0
            elif temp == 'cold':
                start_x = grid_size
            else:
                logging.error('Invalid walker temperature')
                raise SystemExit
        else:
            start_x = np.random.randint(0, grid_size + 1)
        start_y = np.random.randint(0, grid_size + 1)
        start = [start_x, start_y]
        self.pos = [start]

    def add_pos(self, newpos):  # add new position
        self.pos.append(list(newpos))

    def replace_pos(self, newpos):  # replace current position
        self.pos[-1] = list(newpos)

    def erase_prev_pos(self):
        '''Removes all but last position from memory. Important for constant flux simulation
        so that memory usage is kept low.'''
        self.pos = [self.pos[-1]]  # double nesting is important here for the way the item is parsed

