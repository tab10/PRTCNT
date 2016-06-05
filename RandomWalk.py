import numpy as np
import creation
import logging


def runrandomwalk_2d_onlat(grid, timesteps, temp):
    # Start the random walk, for one walker
    walker = creation.Walker2D_onlat(grid.size, temp)
    moves = [[0, 1], [1, 0], [0, -1], [-1, 0]]
    for i in range(1,timesteps+1):
        d_pos = moves[np.random.randint(0, 4)]
        walker.add_dpos(d_pos)
        walker = applytubejumps2d(walker, grid)
        walker = applybdconditions2d(walker, grid)
    logging.info('Done in %d steps' % i)
    return walker


def applytubejumps2d(walker, grid):
    # This code is slow, refactor to reduce load.
    jump_moves = [[0, 1], [1, 0], [0, -1], [-1, 0], [0, 1], [1, 0], [0, -1], [-1, 0]]
    # 8 possible jump positions for now if at tube end, 4 at left (first 4) and 4 at right (second 4)
    tube_l = []
    tube_r = []  # left and right coords of tubes (x,y)
    for i in range(len(grid.tube_coords)):
        tube_l.append(grid.tube_coords[i][0:2])
        tube_r.append(grid.tube_coords[i][2:4])
    cur_pos = list(walker.pos[-1])
    cur_pos_rev = list(np.flipud(walker.pos[-1]))
    for i in range(len(tube_l)):
        if tube_l[i] == (cur_pos or cur_pos_rev):  # coord on left tube end jumps to right end
            # logging.debug("Jump found")
            choice = np.random.randint(0, 8)
            d_pos = jump_moves[choice]
            if choice <= 3:  # stay at left end
                walker.replace_dpos(d_pos)
            else:  # jump across tube to right end
                newpos = np.array(tube_r[i]) + np.array(d_pos)
                walker.replace_pos(newpos)
        elif tube_r[i] == (cur_pos or cur_pos_rev):  # coord on left tube end jumps to right end
            # logging.debug("Jump found")
            choice = np.random.randint(0, 8)
            d_pos = jump_moves[choice]
            if choice <= 3:  # stay at left end
                walker.replace_dpos(d_pos)
            else:  # jump across tube to right end
                newpos = np.array(tube_l[i]) + np.array(d_pos)
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
    if cur_pos[1] > grid.size: # y coordinate
        cur_pos[1] = 1
    elif cur_pos[1] < 0:
        cur_pos[1] = grid.size - 1
    walker.replace_pos(cur_pos)
    return walker