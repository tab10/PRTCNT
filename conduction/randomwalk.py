import logging

import numpy as np

import creation


def runrandomwalk_2d_onlat(grid, timesteps, temp):
    # Start the random walk, for one walker
    walker = creation.Walker2D_onlat(grid.size, temp)
    moves = [[0, 1], [1, 0], [0, -1], [-1, 0]]
    for i in range(1, timesteps + 1):
        d_pos = moves[np.random.randint(0, 4)]
        walker.add_dpos(d_pos)
        walker = applybdconditions2d(walker, grid)
        walker = applytubejumps2d(walker, grid)
    # logging.info('Done in %d steps' % i)
    return walker


def runrandomwalk_3d_onlat(grid, timesteps, temp):
    # Start the random walk, for one walker
    walker = creation.Walker3D_onlat(grid.size, temp)
    moves = [[0, 0, 1], [0, 1, 0], [1, 0, 0], [0, 0, -1], [0, -1, 0], [-1, 0, 0]]
    for i in range(1, timesteps + 1):
        d_pos = moves[np.random.randint(0, 6)]
        walker.add_dpos(d_pos)
        walker = applybdconditions3d(walker, grid)
        walker = applytubejumps3d(walker, grid)
    # logging.info('Done in %d steps' % i)
    return walker


def applytubejumps2d(walker, grid):
    jump_moves = [[0, 1], [1, 0], [0, -1], [-1, 0], [0, 1], [1, 0], [0, -1], [-1, 0]]
    # 8 possible jump positions for now if at tube end, 4 at left (first 4) and 4 at right (second 4)
    tube_l = []
    tube_r = []  # left and right coords of tubes (x,y)
    for i in range(len(grid.tube_coords)):
        tube_l.append(grid.tube_coords[i][0:2])
        tube_r.append(grid.tube_coords[i][2:4])
    cur_pos = list(walker.pos[-1])
    choice = np.random.randint(0, 8)
    d_pos = jump_moves[choice]
    # coord on left tube end jumps to right end
    tube_check_val_l = int(grid.tube_check_l[int(cur_pos[0]), int(cur_pos[1])])
    if tube_check_val_l > 0:
        # logging.debug("Jump found")
        if choice <= 3:  # stay at left end
            walker.replace_dpos(d_pos)
        else:  # jump across tube to right end
            newpos = np.array(tube_r[tube_check_val_l]) + np.array(d_pos)
            walker.replace_pos(newpos)
    # coord on right tube end jumps to left end
    tube_check_val_r = int(grid.tube_check_r[int(cur_pos[0]), int(cur_pos[1])])
    if tube_check_val_r > 0:
        # logging.debug("Jump found")
        if choice <= 3:  # stay at right end
            walker.replace_dpos(d_pos)
        else:  # jump across tube to left end
            newpos = np.array(tube_l[tube_check_val_r]) + np.array(d_pos)
            walker.replace_pos(newpos)
    return walker


def applytubejumps3d(walker, grid):
    jump_moves = [[0, 0, 1], [0, 1, 0], [1, 0, 0], [0, 0, -1], [0, -1, 0], [-1, 0, 0], [0, 0, 1], [0, 1, 0], [1, 0, 0],
                  [0, 0, -1], [0, -1, 0], [-1, 0, 0]]
    # 12 possible jump positions for now if at tube end, 6 at left (first 6) and 6 at right (second 6)
    tube_l = []
    tube_r = []  # left and right coords of tubes (x,y)
    for i in range(len(grid.tube_coords)):
        tube_l.append(grid.tube_coords[i][0:3])
        tube_r.append(grid.tube_coords[i][3:6])
    cur_pos = list(walker.pos[-1])
    choice = np.random.randint(0, 12)
    d_pos = jump_moves[choice]
    # coord on left tube end jumps to right end
    tube_check_val_l = int(grid.tube_check_l[int(cur_pos[0]), int(cur_pos[1]), int(cur_pos[2])])
    if tube_check_val_l > 0:
        # logging.debug("Jump found")
        if choice <= 5:  # stay at left end
            walker.replace_dpos(d_pos)
        else:  # jump across tube to right end
            newpos = np.array(tube_r[tube_check_val_l]) + np.array(d_pos)
            walker.replace_pos(newpos)
    # coord on right tube end jumps to left end
    tube_check_val_r = int(grid.tube_check_r[int(cur_pos[0]), int(cur_pos[1]), int(cur_pos[2])])
    if tube_check_val_r > 0:
        # logging.debug("Jump found")
        if choice <= 5:  # stay at right end
            walker.replace_dpos(d_pos)
        else:  # jump across tube to left end
            newpos = np.array(tube_l[tube_check_val_r]) + np.array(d_pos)
            walker.replace_pos(newpos)
    return walker


def applybdconditions2d(walker, grid):
    # periodic bd conditions for y=0 and y=grid.size
    # reflective bd conditions for x=0 and x=grid.size
    cur_pos = walker.pos[-1]
    if cur_pos[0] >= grid.size:  # x coordinate
        cur_pos[0] -= 1
    elif cur_pos[0] <= 0:
        cur_pos[0] += 1
    if cur_pos[1] > grid.size:  # y coordinate
        cur_pos[1] = 1
    elif cur_pos[1] < 0:
        cur_pos[1] = grid.size - 1
    walker.replace_pos(cur_pos)
    return walker


def applybdconditions3d(walker, grid):
    # periodic bd conditions for y,z=0 and y,z=grid.size
    # reflective bd conditions for x=0 and x=grid.size
    cur_pos = walker.pos[-1]
    if cur_pos[0] >= grid.size:  # x coordinate
        cur_pos[0] -= 1
    elif cur_pos[0] <= 0:
        cur_pos[0] += 1
    if cur_pos[1] > grid.size:  # y coordinate
        cur_pos[1] = 1
    elif cur_pos[1] < 0:
        cur_pos[1] = grid.size - 1
    if cur_pos[2] > grid.size:  # z coordinate
        cur_pos[2] = 1
    elif cur_pos[2] < 0:
        cur_pos[2] = grid.size - 1
    walker.replace_pos(cur_pos)
    return walker
