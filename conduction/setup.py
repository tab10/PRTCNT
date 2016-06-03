import numpy as np


class Grid2D(object):
    def __init__(self, grid_size, tube_length, num_tubes, mean_dist_tubes, exact_sol, orientation):
        if tube_length > grid_size:
            print 'Nanotube is too large for grid.'
            raise SystemExit
        self.tube_coords = []
        self.tube_centers = []
        if num_tubes > 0:  # tubes exist
            for i in range(num_tubes):  # currently no mean dist used, ADD LATER?
                self.tube_centers.append((np.random.randint(-grid_size + tube_length, grid_size - tube_length),
                                       np.random.randint(-grid_size + tube_length,
                                                         grid_size - tube_length)))
                # reduced to ensure tubes stay in grid
                self.tube_coords.append((self.generate_2d_tube(self.tube_centers[i][0], self.tube_centers[i][1],
                                                                   tube_length, grid_size, orientation, exact_sol)))

    @staticmethod
    def generate_2d_tube(x_c, y_c, radius, grid_size, orientation, int_sols=False):  # integer solutions finds exact solutions only
        # for endpoints with given radius, otherwise rounded, giving more possibilities
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
                if orientation == 'random':
                    angle = 2.0 * np.pi * np.random.random_sample()  # randomly chosen
                elif orientation == 'vertical':
                    angle = np.pi/2.0
                elif orientation == 'horizontal':
                    angle = np.pi
                else: # assuming input in DEGREES
                    angle = (float(orientation)/360.0) * 2 * np.pi
                x_vect = radius * np.cos(angle)
                y_vect = np.sqrt(radius ** 2 - x_vect ** 2)
                x_l = x_c - int(round(x_vect))
                x_r = x_c + int(round(x_vect))
                y_l = y_c - int(round(y_vect))
                y_r = y_c + int(round(y_vect))
        return x_l, y_l, x_r, y_r


def taxicab_dist(x0, y0, x1, y1):
    dist = np.abs(x1 - x0) + np.abs(y1 - y0)
    return dist


def radius(x, y):
    radius = np.sqrt(x ** 2 + y ** 2)
    return radius


def euc_dist(x0, y0, x1, y1):
    dist = np.sqrt((x1 - x0) ** 2 + (y1 - y0) ** 2)
    return dist


class Walker2D(object):
    def __init__(self, grid_size, walker_start_dir):
        if walker_start_dir == 'left':
            start_x = -grid_size
            start_y = np.random.randint(-grid_size, grid_size)
        elif walker_start_dir == 'right':
            start_x = grid_size
            start_y = np.random.randint(-grid_size, grid_size)
        elif walker_start_dir == 'top':
            start_x = np.random.randint(-grid_size, grid_size)
            start_y = grid_size
        elif walker_start_dir == 'bottom':
            start_x = np.random.randint(-grid_size, grid_size)
            start_y = -grid_size
        else:
            print('Not correct 2D start option for walkers')
            raise SystemExit
        start = [start_x, start_y]
        self.pos = start

    def updatepos(self,newpos):
        arr1 = np.array(self.pos)
        arr2 = np.array(newpos)
        self.pos = arr1 + arr2

def setupwalkers2d(num_walkers, grid_size, walker_start_dir):
    pos_all = [] # master list holding walker objects at every time step
    walkers_timestep = [] # this will be for each timestep
    for i in range(num_walkers):
        walkers_timestep.append(Walker2D(grid_size, walker_start_dir))
    pos_all.append(walkers_timestep)
    return pos_all


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