import logging

import numpy as np

import creation


def runrandomwalk_2d_onlat(grid, timesteps, temp):
    # Start the random walk, for one walker
    walker = creation.Walker2D_onlat(grid.size, temp)
    moves = [[0, 1], [1, 0], [0, -1], [-1, 0]]
    for i in range(1, timesteps + 1):
        d_pos = moves[np.random.randint(0, 4)]
        walker, jumped = applybdconditions2d(walker, grid)
        walker, jumped = applytubejumps2d(walker, grid, jumped)
        if not jumped:
            walker.add_dpos(d_pos)
    # logging.info('Done in %d steps' % i)
    return walker


def runrandomwalk_3d_onlat(grid, timesteps, temp):
    # Start the random walk, for one walker
    walker = creation.Walker3D_onlat(grid.size, temp)
    moves = [[0, 0, 1], [0, 1, 0], [1, 0, 0], [0, 0, -1], [0, -1, 0], [-1, 0, 0]]
    for i in range(1, timesteps + 1):
        d_pos = moves[np.random.randint(0, 6)]
        walker, jumped = applybdconditions3d(walker, grid)
        walker, jumped = applytubejumps3d(walker, grid, jumped)
        if not jumped:
            walker.add_dpos(d_pos)
    # logging.info('Done in %d steps' % i)
    return walker


def applytubejumps2d(walker, grid, jumped):
    if not jumped:
        jump_moves = [[0, 1], [1, 0], [0, -1], [-1, 0], [0, 1], [1, 0], [0, -1], [-1, 0]]
        # 8 possible jump positions for now if at tube end, 4 at left (first 4) and 4 at right (second 4)
        cur_pos = list(walker.pos[-1])
        # print cur_pos
        choice = np.random.randint(0, 8)
        d_pos = jump_moves[choice]
        # coord on left tube end jumps to right end
        tube_check_val_l = grid.tube_check_l[int(cur_pos[0]), int(cur_pos[1])]
        if tube_check_val_l > 0:
            check = True
            # logging.debug("Jump found")
            if choice <= 3:  # stay at left end
                walker.replace_dpos(d_pos)
            else:  # jump across tube to right end
                newpos = np.array(grid.tube_coords_r[tube_check_val_l]) + np.array(d_pos)
                walker.replace_pos(newpos)
        else:
            check = False
        # coord on right tube end jumps to left end
        tube_check_val_r = grid.tube_check_r[int(cur_pos[0]), int(cur_pos[1])]
        if tube_check_val_r > 0:
            check = True
            # logging.debug("Jump found")
            if choice <= 3:  # stay at right end
                walker.replace_dpos(d_pos)
            else:  # jump across tube to left end
                newpos = np.array(grid.tube_coords_l[tube_check_val_r]) + np.array(d_pos)
                walker.replace_pos(newpos)
        else:
            check = False
    check = jumped
    return walker, check


def applytubejumps3d(walker, grid, jumped):
    if not jumped:
        jump_moves = [[0, 0, 1], [0, 1, 0], [1, 0, 0], [0, 0, -1], [0, -1, 0], [-1, 0, 0], [0, 0, 1], [0, 1, 0],
                      [1, 0, 0],
                      [0, 0, -1], [0, -1, 0], [-1, 0, 0]]
        # 12 possible jump positions for now if at tube end, 6 at left (first 6) and 6 at right (second 6)
        cur_pos = list(walker.pos[-1])
        choice = np.random.randint(0, 12)
        d_pos = jump_moves[choice]
        # coord on left tube end jumps to right end
        tube_check_val_l = grid.tube_check_l[int(cur_pos[0]), int(cur_pos[1]), int(cur_pos[2])]
        if tube_check_val_l > 0:
            check = True
            # logging.debug("Jump found")
            if choice <= 5:  # stay at left end
                walker.replace_dpos(d_pos)
            else:  # jump across tube to right end
                newpos = np.array(grid.tube_coords_r[tube_check_val_l]) + np.array(d_pos)
                walker.replace_pos(newpos)
        else:
            check = False
        # coord on right tube end jumps to left end
        tube_check_val_r = grid.tube_check_r[int(cur_pos[0]), int(cur_pos[1]), int(cur_pos[2])]
        if tube_check_val_r > 0:
            check = True
            # logging.debug("Jump found")
            if choice <= 5:  # stay at right end
                walker.replace_dpos(d_pos)
            else:  # jump across tube to left end
                newpos = np.array(grid.tube_coords_l[tube_check_val_r]) + np.array(d_pos)
                walker.replace_pos(newpos)
        else:
            check = False
    check = jumped
    return walker, check


def applybdconditions2d(walker, grid):
    # periodic bd conditions for y=0 and y=grid.size
    # reflective bd conditions for x=0 and x=grid.size
    cur_pos = walker.pos[-1]
    #print cur_pos
    if cur_pos[0] >= grid.size:  # x coordinate
        cur_pos[0] = grid.size - 1
        jumped = True
    elif cur_pos[0] <= 0:
        cur_pos[0] = 1
        jumped = True
    else:
        jumped = False
    if cur_pos[1] > grid.size:  # y coordinate
        cur_pos[1] = 1
        jumped = True
    elif cur_pos[1] < 0:
        cur_pos[1] = grid.size - 1
        jumped = True
    else:
        jumped = False
    walker.replace_pos(cur_pos)
    return walker, jumped


def applybdconditions3d(walker, grid):
    # periodic bd conditions for y,z=0 and y,z=grid.size
    # reflective bd conditions for x=0 and x=grid.size
    cur_pos = walker.pos[-1]
    if cur_pos[0] >= grid.size:  # x coordinate
        cur_pos[0] = grid.size - 1
        jumped = True
    elif cur_pos[0] <= 0:
        cur_pos[0] = 1
        jumped = True
    else:
        jumped = False
    if cur_pos[1] > grid.size:  # y coordinate
        cur_pos[1] = 1
        jumped = True
    elif cur_pos[1] < 0:
        cur_pos[1] = grid.size - 1
        jumped = True
    else:
        jumped = False
    if cur_pos[2] > grid.size:  # z coordinate
        cur_pos[2] = 1
        jumped = True
    elif cur_pos[2] < 0:
        cur_pos[2] = grid.size - 1
        jumped = True
    else:
        jumped = False
    walker.replace_pos(cur_pos)
    return walker, jumped
