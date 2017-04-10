# //////////////////////////////////////////////////////////////////////////////////// #
# ////////////////////////////// ##  ##  ###  ## ### ### ///////////////////////////// #
# ////////////////////////////// # # # #  #  #   # #  #  ///////////////////////////// #
# ////////////////////////////// ##  ##   #  #   # #  #  ///////////////////////////// #
# ////////////////////////////// #   # #  #  #   # #  #  ///////////////////////////// #
# ////////////////////////////// #   # #  #   ## # #  #  ///////////////////////////// #
# ////////////////////////////// ###  #          ##           # ////////////////////// #
# //////////////////////////////  #      ###     # # # # ### ### ///////////////////// #
# //////////////////////////////  #   #  ###     ##  # # #    # ////////////////////// #
# //////////////////////////////  #   ## # #     # # ### #    ## ///////////////////// #
# //////////////////////////////  #              ## ////////////////////////////////// #
# //////////////////////////////////////////////////////////////////////////////////// #


from __future__ import division
import numpy as np
import logging
from conduction import *


class Grid3D_onlat(object):
    def __init__(self, grid_size, tube_length, num_tubes, orientation, tube_radius, parallel, plot_save_dir,
                 disable_func, rules_test, rank=None,
                 size=None):
        """Grid in first quadrant only for convenience"""
        self.size = grid_size
        self.tube_radius = tube_radius
        # serial implementation
        status_counter = 0
        counter = 0
        self.tube_coords = []
        self.tube_coords_l = []
        self.tube_coords_r = []
        self.tube_centers = []
        self.theta = []
        self.phi = []
        self.setup_tube_vol_check_array_3d()
        if rules_test:
            bound = [20, 20, 20]  # all periodic
        else:
            bound = [10, 20, 20]  # X reflective, YZ periodic
        self.bound = bound
        if tube_radius == 0:
            logging.info("Zero tube radius given. Tubes will have no volume.")
            fill_fract = 2.0 * float(num_tubes) / grid_size ** 3
            logging.info("Filling fraction is %.2f %%" % (fill_fract * 100.0))
            backend.save_fill_frac(plot_save_dir, fill_fract)
            if num_tubes > 0:  # tubes exist
                for i in range(num_tubes):  # currently no mean dist used, ADD LATER?
                    if (i % 50) == 0:
                        status_counter += 50
                        logging.info('Generating tube %d...' % (status_counter - 50))
                    x_l, y_l, z_l, x_r, y_r, z_r, x_c, y_c, z_c, theta, phi = self.generate_3d_tube(tube_length,
                                                                                                    orientation,
                                                                                                    tube_radius)
                    self.tube_centers.append([x_c, y_c, z_c])
                    self.tube_coords.append([x_l, y_l, z_l, x_r, y_r, z_r])
                    self.tube_coords_l.append([x_l, y_l, z_l])
                    self.tube_coords_r.append([x_r, y_r, z_r])
                    self.theta.append(theta)
                    self.phi.append(phi)
                    if i == 0:
                        self.add_tube_vol_check_array_3d([x_l, y_l, z_l, x_r, y_r, z_r], None, disable_func)
                    if i >= 1:
                        uni_flag = self.check_tube_unique_3d_arraymethod([x_l, y_l, z_l, x_r, y_r, z_r])
                        # uni_flag = self.check_tube_unique(self.tube_coords, False)
                        while not uni_flag:
                            counter += 1
                            self.tube_centers.pop()
                            self.tube_coords.pop()
                            self.tube_coords_l.pop()
                            self.tube_coords_r.pop()
                            self.theta.pop()
                            self.phi.pop()
                            x_l, y_l, z_l, x_r, y_r, z_r, x_c, y_c, z_c, theta, phi = self.generate_3d_tube(
                                tube_length,
                                orientation,
                                tube_radius)
                            self.tube_centers.append([x_c, y_c, z_c])
                            self.tube_coords.append([x_l, y_l, z_l, x_r, y_r, z_r])
                            self.tube_coords_l.append([x_l, y_l, z_l])
                            self.tube_coords_r.append([x_r, y_r, z_r])
                            self.theta.append(theta)
                            self.phi.append(phi)
                            uni_flag = self.check_tube_unique_3d_arraymethod([x_l, y_l, z_l, x_r, y_r, z_r])
                            # uni_flag = self.check_tube_unique(self.tube_coords, False)
                        self.add_tube_vol_check_array_3d([x_l, y_l, z_l, x_r, y_r, z_r], None, disable_func)
                logging.info("Tube generation complete")
                logging.info("Corrected %d overlapping tube endpoints" % counter)
            self.tube_check_l, self.tube_check_r, self.tube_check_bd = self.generate_tube_check_array_3d(rules_test)
        else:
            logging.info("Non-zero tube radius given. Tubes will have excluded volume.")
            l_d = tube_length / (2 * tube_radius)
            logging.info("L/D is %.4f." % l_d)
            self.tube_squares = []
            self.setup_tube_vol_check_array_3d()
            if num_tubes > 0:  # tubes exist
                for i in range(num_tubes):  # currently no mean dist used, ADD LATER?
                    if (i % 50) == 0:
                        status_counter += 50
                        logging.info('Generating tube %d...' % (status_counter - 50))
                    x_l, y_l, z_l, x_r, y_r, z_r, x_c, y_c, z_c, theta, phi = self.generate_3d_tube(tube_length,
                                                                                                    orientation,
                                                                                                    tube_radius)
                    tube_squares = self.find_cubes([x_l, y_l, z_l], [x_r, y_r, z_r])
                    self.tube_centers.append([x_c, y_c, z_c])
                    self.tube_coords.append([x_l, y_l, z_l, x_r, y_r, z_r])
                    self.tube_coords_l.append([x_l, y_l, z_l])
                    self.tube_coords_r.append([x_r, y_r, z_r])
                    self.theta.append(theta)
                    self.phi.append(phi)
                    self.tube_squares.append(tube_squares)
                    if i == 0:
                        self.add_tube_vol_check_array_3d([x_l, y_l, z_l, x_r, y_r, z_r], tube_squares, disable_func)
                    if i >= 1:
                        uni_flag = self.check_tube_and_vol_unique_3d_arraymethod(tube_squares)
                        while not uni_flag:
                            counter += 1
                            self.tube_centers.pop()
                            self.tube_coords.pop()
                            self.tube_coords_l.pop()
                            self.tube_coords_r.pop()
                            self.theta.pop()
                            self.phi.pop()
                            self.tube_squares.pop()
                            x_l, y_l, z_l, x_r, y_r, z_r, x_c, y_c, z_c, theta, phi = self.generate_3d_tube(
                                tube_length,
                                orientation,
                                tube_radius)
                            tube_squares = self.find_cubes([x_l, y_l, z_l], [x_r, y_r, z_r])
                            self.tube_centers.append([x_c, y_c, z_c])
                            self.tube_coords.append([x_l, y_l, z_l, x_r, y_r, z_r])
                            self.tube_coords_l.append([x_l, y_l, z_l])
                            self.tube_coords_r.append([x_r, y_r, z_r])
                            self.theta.append(theta)
                            self.phi.append(phi)
                            self.tube_squares.append(tube_squares)
                            uni_flag = self.check_tube_and_vol_unique_3d_arraymethod(tube_squares)
                        self.add_tube_vol_check_array_3d([x_l, y_l, z_l, x_r, y_r, z_r], tube_squares, disable_func)
            logging.info("Tube generation complete")
            logging.info("Corrected %d overlapping tube endpoints" % counter)
            # get number of squares filled
            cube_count = 0  # each cube has volume 1
            for i in range(len(self.tube_squares)):
                cube_count += len(self.tube_squares[i])
            fill_fract = float(cube_count) * 2.0 * tube_radius / grid_size ** 3
            # each cube has area 1, times the tube radius (important if not 1)
            logging.info("Filling fraction is %.2f %%" % (fill_fract * 100.0))
            backend.save_fill_frac(plot_save_dir, fill_fract)
            self.tube_check_l, self.tube_check_r, self.tube_check_bd = self.generate_tube_check_array_3d(rules_test)
            self.tube_check_bd_vol, self.tube_check_index = self.generate_vol_check_array_3d(disable_func)
        self.avg_tube_len, self.std_tube_len, self.tube_lengths = self.check_tube_lengths()
        logging.info("Actual tube length avg+std: %.4f +- %.4f" % (self.avg_tube_len, self.std_tube_len))

    def generate_3d_tube(self, radius, orientation, tube_radius):
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
                theta_angle_range = list(range(0, 360))  # y-z plane
                phi_angle_range = list(range(0, 360))  # x-y plane
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
            x_r = int(round(self.coord(radius, theta_angle, phi_angle)[0] + x_l))
            y_r = int(round(self.coord(radius, theta_angle, phi_angle)[1] + y_l))
            z_r = int(round(self.coord(radius, theta_angle, phi_angle)[2] + z_l))
            if (x_r > 0) and (x_r < self.size) and (y_r > 0) and (y_r < self.size) and (z_r > 0) \
                    and (z_r < self.size):
                inside_box = True
        if x_l > x_r:  # this imposes left to right tube order WRT + x-axis
            x_l_temp = x_l
            x_r_temp = x_r
            y_l_temp = y_l
            y_r_temp = y_r
            z_l_temp = z_l
            z_r_temp = z_r
            x_l = x_r_temp
            x_r = x_l_temp
            y_l = y_r_temp
            y_r = y_l_temp
            z_l = z_r_temp
            z_r = z_l_temp
        x_c = round(self.coord(radius, theta_angle, phi_angle)[0] + x_l) / 2
        y_c = round(self.coord(radius, theta_angle, phi_angle)[1] + y_l) / 2
        z_c = round(self.coord(radius, theta_angle, phi_angle)[2] + z_l) / 2
        return x_l, y_l, z_l, x_r, y_r, z_r, x_c, y_c, z_c, theta_angle, phi_angle

    def find_cubes(self, start, end):
        """Bresenham's Line Algorithm in 3D
        Produces a list of tuples (bottom left corners of grid)
        All cubes a tube passes through
        Also returns original start and end points in first and last positions respectively
        """
        p1 = np.asarray(start, dtype=float)
        p2 = np.asarray(end, dtype=float)
        p = p1.astype(float)
        d = p2 - p1
        N = int(max(abs(d)))
        s = d / float(N)
        points = []
        for i in range(N):
            p += s
            points.append(list(np.round(p).astype(int)))
        if list(p1) not in points:
            points.insert(0, list(p1.astype(int)))
        if list(p2) not in points:
            points.insert(-1, list(p2.astype(int)))
        # fig = plt.figure()
        # ax = fig.add_subplot(111, projection='3d')
        # for i in range(len(points)):
        #     ax.scatter(points[i][0],points[i][1],points[i][2])
        #     ax.set_xlabel('X')
        #     ax.set_ylabel('Y')
        #     ax.set_zlabel('Z')
        # plt.grid()
        # plt.show()
        return points

    def generate_tube_check_array_3d(self, rules_test):
        tube_check_l = np.zeros((self.size + 1, self.size + 1, self.size + 1), dtype=int)
        tube_check_r = np.zeros((self.size + 1, self.size + 1, self.size + 1), dtype=int)
        bd = np.zeros((self.size + 1, self.size + 1, self.size + 1), dtype=int)
        for i in range(len(self.tube_coords)):
            tube_check_l[self.tube_coords[i][0], self.tube_coords[i][1], self.tube_coords[i][2]] = i + 1
            # THESE ARE OFFSET BY ONE
            tube_check_r[self.tube_coords[i][3], self.tube_coords[i][4], self.tube_coords[i][5]] = i + 1
            # THESE ARE OFFSET BY ONE
            bd[self.tube_coords[i][0], self.tube_coords[i][1], self.tube_coords[i][2]] = 1
            bd[self.tube_coords[i][3], self.tube_coords[i][4], self.tube_coords[i][5]] = 1
            # holds index of tube_coords, if a walker on that position has a nonzero value in this array,
            # pull the right or left tube endpoint (array positions are at left and right endpoints respectively)
        # add boundary tags
        box_dims = [0, self.size]
        for i in range(self.size + 1):
            for j in range(self.size + 1):
                for k in range(self.size + 1):
                    if (i in box_dims) or (j in box_dims) or (k in box_dims):
                        bd[i, j, k] = -1000  # a boundary. no cnt volume or ends can be here. all 6 choices
                        # generated, running through bd function
        return tube_check_l, tube_check_r, bd

    def generate_vol_check_array_3d(self, disable_func):
        """To be used with tube volume
        Generates a boundary/volume lookup array (0 nothing, 1 boundary, -1 volume)"""
        bd_vol = np.zeros((self.size + 1, self.size + 1, self.size + 1), dtype=int)
        index = np.zeros((self.size + 1, self.size + 1, self.size + 1), dtype=int)
        if disable_func:
            endpoint_val = -1  # treat endpoints as volume, changing the rules in the walk
        else:
            endpoint_val = 1  # leave it as endpoint
        for i in range(len(self.tube_coords)):
            bd_vol[
                self.tube_coords[i][0], self.tube_coords[i][1], self.tube_coords[i][2]] = endpoint_val  # left endpoints
            bd_vol[self.tube_coords[i][3], self.tube_coords[i][4], self.tube_coords[i][
                5]] = endpoint_val  # right endpoints
            index[self.tube_coords[i][0], self.tube_coords[i][1], self.tube_coords[i][2]] = i + 1
            # THESE ARE OFFSET BY ONE
            index[self.tube_coords[i][3], self.tube_coords[i][4], self.tube_coords[i][5]] = i + 1
            # THESE ARE OFFSET BY ONE
            for j in range(1, len(self.tube_squares[i]) - 1):
                bd_vol[self.tube_squares[i][j][0], self.tube_squares[i][j][1], self.tube_squares[i][j][
                    2]] = -1  # volume points
                index[self.tube_squares[i][j][0], self.tube_squares[i][j][1], self.tube_squares[i][j][
                    2]] = i + 1  # THESE ARE OFFSET BY ONE
                # add boundary tags
                # entire planes hold bd conditions, so start at a corner and fill every plane
        box_dims = [0, self.size]
        for i in range(self.size + 1):
            for j in range(self.size + 1):
                for k in range(self.size + 1):
                    if (i in box_dims) or (j in box_dims) or (k in box_dims):
                        bd_vol[i, j, k] = -1000  # a boundary. no cnt volume or ends can be here. all 6 choices
                        # generated, running through bd function
        return bd_vol, index

    def setup_tube_vol_check_array_3d(self):
        "Setup of tube check arrays, returns nothing"
        bd_vol = np.zeros((self.size + 1, self.size + 1, self.size + 1), dtype=int)
        index = np.zeros((self.size + 1, self.size + 1, self.size + 1), dtype=int)
        self.tube_check_bd_vol = bd_vol
        self.tube_check_index = index

    def add_tube_vol_check_array_3d(self, new_tube_coords, new_tube_squares, disable_func):
        "Adds tube to the current check arrays"
        index_val = len(self.tube_coords) + 1  # THESE ARE OFFSET BY ONE
        if disable_func:
            endpoint_val = -1  # treat endpoints as volume, changing the rules in the walk
        else:
            endpoint_val = 1  # leave it as endpoint
        self.tube_check_bd_vol[
            new_tube_coords[0], new_tube_coords[1], new_tube_coords[2]] = endpoint_val  # left endpoints
        self.tube_check_bd_vol[
            new_tube_coords[3], new_tube_coords[4], new_tube_coords[5]] = endpoint_val  # right endpoints
        self.tube_check_index[new_tube_coords[0], new_tube_coords[1], new_tube_coords[2]] = index_val
        self.tube_check_index[new_tube_coords[3], new_tube_coords[4], new_tube_coords[5]] = index_val
        if new_tube_squares is not None:
            for j in range(1, len(new_tube_squares) - 1):
                self.tube_check_bd_vol[
                    new_tube_squares[j][0], new_tube_squares[j][1], new_tube_squares[j][2]] = -1  # volume points
                self.tube_check_index[
                    new_tube_squares[j][0], new_tube_squares[j][1], new_tube_squares[j][2]] = index_val

    def check_tube_unique_3d_arraymethod(self, new_tube_squares):
        "No volume"
        uni_flag = True
        check_l = self.tube_check_bd_vol[new_tube_squares[0], new_tube_squares[1], new_tube_squares[2]]
        check_r = self.tube_check_bd_vol[new_tube_squares[3], new_tube_squares[4], new_tube_squares[5]]
        if (check_l != 0) or (check_r != 0):
            uni_flag = False
        return uni_flag

    def check_tube_and_vol_unique_3d_arraymethod(self, new_tube_squares):
        "Volume"
        uni_flag = True
        index_val = len(self.tube_coords) + 1
        for l in range(len(new_tube_squares)):
            test_vol = self.tube_check_bd_vol[new_tube_squares[l][0], new_tube_squares[l][1], new_tube_squares[l][2]]
            if (test_vol == 1) or (test_vol == -1):  # new tube overlaps old volume or endpoint
                uni_flag = False
            # generate candidate check positions for tube crossing
            # this algorithm looks for interweaved diagonal clusters. 6 cases in 3D in array method.
            # x-y plane
            # cur(tl,br) old(tr,bl) pts wrt current tl
            b_r = [new_tube_squares[l][0] + 1, new_tube_squares[l][1] - 1, new_tube_squares[l][2]]
            old_t_r = self.tube_check_bd_vol[
                new_tube_squares[l][0] + 1, new_tube_squares[l][1], new_tube_squares[l][2]]
            old_b_l = self.tube_check_bd_vol[
                new_tube_squares[l][0], new_tube_squares[l][1] - 1, new_tube_squares[l][2]]
            old_t_r_index = self.tube_check_index[
                new_tube_squares[l][0] + 1, new_tube_squares[l][1], new_tube_squares[l][2]]
            old_b_l_index = self.tube_check_index[
                new_tube_squares[l][0], new_tube_squares[l][1] - 1, new_tube_squares[l][2]]
            if ((old_t_r == 1) or (old_t_r == -1)) and ((old_b_l == 1) or (old_b_l == -1)) \
                    and (b_r in new_tube_squares) and (old_t_r_index == old_b_l_index):  # we have a crossing
                uni_flag = False
            # old(tl,br) cur(tr,bl) pts wrt current tr
            b_l = [new_tube_squares[l][0] - 1, new_tube_squares[l][1] - 1, new_tube_squares[l][2]]
            old_t_l = self.tube_check_bd_vol[
                new_tube_squares[l][0] - 1, new_tube_squares[l][1], new_tube_squares[l][2]]
            old_b_r = self.tube_check_bd_vol[
                new_tube_squares[l][0], new_tube_squares[l][1] - 1, new_tube_squares[l][2]]
            old_t_l_index = self.tube_check_index[
                new_tube_squares[l][0] - 1, new_tube_squares[l][1], new_tube_squares[l][2]]
            old_b_r_index = self.tube_check_index[
                new_tube_squares[l][0], new_tube_squares[l][1] - 1, new_tube_squares[l][2]]
            if ((old_t_l == 1) or (old_t_l == -1)) and ((old_b_r == 1) or (old_b_r == -1)) \
                    and (b_l in new_tube_squares) and (old_t_l_index == old_b_r_index):  # we have a crossing
                uni_flag = False
            # y-z plane
            # cur(tl,br) old(tr,bl) pts wrt current tl
            b_r = [new_tube_squares[l][0], new_tube_squares[l][1] + 1, new_tube_squares[l][2] - 1]
            old_t_r = self.tube_check_bd_vol[
                new_tube_squares[l][0], new_tube_squares[l][1] + 1, new_tube_squares[l][2]]
            old_b_l = self.tube_check_bd_vol[
                new_tube_squares[l][0], new_tube_squares[l][1], new_tube_squares[l][2] - 1]
            old_t_r_index = self.tube_check_index[
                new_tube_squares[l][0], new_tube_squares[l][1] + 1, new_tube_squares[l][2]]
            old_b_l_index = self.tube_check_index[
                new_tube_squares[l][0], new_tube_squares[l][1], new_tube_squares[l][2] - 1]
            if ((old_t_r == 1) or (old_t_r == -1)) and ((old_b_l == 1) or (old_b_l == -1)) \
                    and (b_r in new_tube_squares) and (old_t_r_index == old_b_l_index):  # we have a crossing
                uni_flag = False
            # old(tl,br) cur(tr,bl) pts wrt current tr
            b_l = [new_tube_squares[l][0], new_tube_squares[l][1] - 1, new_tube_squares[l][2] - 1]
            old_t_l = self.tube_check_bd_vol[
                new_tube_squares[l][0], new_tube_squares[l][1] - 1, new_tube_squares[l][2]]
            old_b_r = self.tube_check_bd_vol[
                new_tube_squares[l][0], new_tube_squares[l][1], new_tube_squares[l][2] - 1]
            old_t_l_index = self.tube_check_index[
                new_tube_squares[l][0], new_tube_squares[l][1] - 1, new_tube_squares[l][2]]
            old_b_r_index = self.tube_check_index[
                new_tube_squares[l][0], new_tube_squares[l][1], new_tube_squares[l][2] - 1]
            if ((old_t_l == 1) or (old_t_l == -1)) and ((old_b_r == 1) or (old_b_r == -1)) \
                    and (b_l in new_tube_squares) and (old_t_l_index == old_b_r_index):  # we have a crossing
                uni_flag = False
            # x-z plane
            # cur(tl,br) old(tr,bl) pts wrt current tl
            b_r = [new_tube_squares[l][0] + 1, new_tube_squares[l][1], new_tube_squares[l][2] - 1]
            old_t_r = self.tube_check_bd_vol[
                new_tube_squares[l][0] + 1, new_tube_squares[l][1], new_tube_squares[l][2]]
            old_b_l = self.tube_check_bd_vol[
                new_tube_squares[l][0], new_tube_squares[l][1], new_tube_squares[l][2] - 1]
            old_t_r_index = self.tube_check_index[
                new_tube_squares[l][0] + 1, new_tube_squares[l][1], new_tube_squares[l][2]]
            old_b_l_index = self.tube_check_index[
                new_tube_squares[l][0], new_tube_squares[l][1], new_tube_squares[l][2] - 1]
            if ((old_t_r == 1) or (old_t_r == -1)) and ((old_b_l == 1) or (old_b_l == -1)) \
                    and (b_r in new_tube_squares) and (old_t_r_index == old_b_l_index):  # we have a crossing
                uni_flag = False
            # old(tl,br) cur(tr,bl) pts wrt current tr
            b_l = [new_tube_squares[l][0] - 1, new_tube_squares[l][1], new_tube_squares[l][2] - 1]
            old_t_l = self.tube_check_bd_vol[
                new_tube_squares[l][0] - 1, new_tube_squares[l][1], new_tube_squares[l][2]]
            old_b_r = self.tube_check_bd_vol[
                new_tube_squares[l][0], new_tube_squares[l][1], new_tube_squares[l][2] - 1]
            old_t_l_index = self.tube_check_index[
                new_tube_squares[l][0] - 1, new_tube_squares[l][1], new_tube_squares[l][2]]
            old_b_r_index = self.tube_check_index[
                new_tube_squares[l][0], new_tube_squares[l][1], new_tube_squares[l][2] - 1]
            if ((old_t_l == 1) or (old_t_l == -1)) and ((old_b_r == 1) or (old_b_r == -1)) \
                    and (b_l in new_tube_squares) and (old_t_l_index == old_b_r_index):  # we have a crossing
                uni_flag = False
        return uni_flag

    def check_tube_lengths(self):
        tube_lengths = np.zeros(len(self.tube_coords))
        for i in range(len(self.tube_coords)):
            dist = self.euc_dist(self.tube_coords[i][0], self.tube_coords[i][1], self.tube_coords[i][2],
                                 self.tube_coords[i][3], self.tube_coords[i][4], self.tube_coords[i][5])
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
        cur_l = [current[0], current[1], current[2]]
        cur_r = [current[3], current[4], current[5]]
        # separate
        tube_l = []
        tube_r = []
        if not parallel:
            for i in range(len(coords_list) - 1):  # -1 accounts for not including the current tube
                tube_l.append([coords_list[i][0], coords_list[i][1], coords_list[i][2]])
                tube_r.append([coords_list[i][3], coords_list[i][4], coords_list[i][5]])
        else:
            exclude_val = len(coords_list) - size + rank
            for i in range(len(coords_list)):
                if i != exclude_val:
                    tube_l.append([coords_list[i][0], coords_list[i][1], coords_list[i][2]])
                    tube_r.append([coords_list[i][3], coords_list[i][4], coords_list[i][5]])
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
            for i in range(len(tube_squares) - 1):  # -1 accounts for not including the current tube
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
                # this algorithm looks for interweaved diagonal clusters. The only 6 possibilities are checked
                # x-y plane
                t_l = current_vol[l]
                t_r = [current_vol[l][0] + 1, current_vol[l][1], current_vol[l][2]]
                b_l = [current_vol[l][0], current_vol[l][1] - 1, current_vol[l][2]]
                b_r = [current_vol[l][0] + 1, current_vol[l][1] - 1, current_vol[l][2]]
                if (b_r in current_vol) and (t_r in vol) and (b_l in vol):
                    uni_flag = False
                t_r = current_vol[l]
                t_l = [current_vol[l][0] - 1, current_vol[l][1], current_vol[l][2]]
                b_l = [current_vol[l][0] - 1, current_vol[l][1] - 1, current_vol[l][2]]
                b_r = [current_vol[l][0], current_vol[l][1] - 1, current_vol[l][2]]
                if (t_l in vol) and (b_l in current_vol) and (b_r in vol):
                    uni_flag = False
                # y-z plane
                t_l = current_vol[l]
                t_r = [current_vol[l][0], current_vol[l][1] + 1, current_vol[l][2]]
                b_l = [current_vol[l][0], current_vol[l][1], current_vol[l][2] - 1]
                b_r = [current_vol[l][0], current_vol[l][1] + 1, current_vol[l][2] - 1]
                if (b_r in current_vol) and (t_r in vol) and (b_l in vol):
                    uni_flag = False
                t_r = current_vol[l]
                t_l = [current_vol[l][0], current_vol[l][1] - 1, current_vol[l][2]]
                b_l = [current_vol[l][0], current_vol[l][1] - 1, current_vol[l][2] - 1]
                b_r = [current_vol[l][0], current_vol[l][1], current_vol[l][2] - 1]
                if (t_l in vol) and (b_l in current_vol) and (b_r in vol):
                    uni_flag = False
                # x-z plane
                t_l = current_vol[l]
                t_r = [current_vol[l][0] + 1, current_vol[l][1], current_vol[l][2]]
                b_l = [current_vol[l][0], current_vol[l][1], current_vol[l][2] - 1]
                b_r = [current_vol[l][0] + 1, current_vol[l][1], current_vol[l][2] - 1]
                if (b_r in current_vol) and (t_r in vol) and (b_l in vol):
                    uni_flag = False
                t_r = current_vol[l]
                t_l = [current_vol[l][0] - 1, current_vol[l][1], current_vol[l][2]]
                b_l = [current_vol[l][0] - 1, current_vol[l][1], current_vol[l][2] - 1]
                b_r = [current_vol[l][0], current_vol[l][1], current_vol[l][2] - 1]
                if (t_l in vol) and (b_l in current_vol) and (b_r in vol):
                    uni_flag = False
                    ##########
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
    def calc_radius(x, y, z):
        radius = np.sqrt(x ** 2 + y ** 2 + z ** 2)
        return radius

    @staticmethod
    def euc_dist(x0, y0, z0, x1, y1, z1):
        dist = np.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2 + (z1 - z0) ** 2)
        return dist


class Walker3D_onlat(object):
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
        start_z = np.random.randint(0, grid_size + 1)
        start = [start_x, start_y, start_z]
        self.pos = [start]

    def add_pos(self, newpos):  # add new position
        self.pos.append(list(newpos))

    def replace_pos(self, newpos):  # replace current position
        self.pos[-1] = list(newpos)

    def erase_prev_pos(self):
        '''Removes all but last position from memory. Important for constant flux simulation
        so that memory usage is kept low.'''
        self.pos = [self.pos[-1]]  # double nesting is important here for the way the item is parsed
