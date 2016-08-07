import logging

import numpy as np

import creation


def runrandomwalk_2d_onlat(grid, timesteps, temp, kapitza, prob_m_cn):
    # Start the random walk, for one walker
    walker = creation.Walker2D_onlat(grid.size, temp)
    moves = [[0, 1], [1, 0], [0, -1], [-1, 0]]
    for i in range(1, timesteps + 1):
        # check where we are
        cur_pos = np.asarray(walker.pos[-1])
        # DEBUG
        # print i, cur_pos
        # Having tube radius doesn't matter if kapitza is off, apart from excluded volume
        if kapitza:
            cur_type = grid.tube_check_bd_vol[cur_pos[0], cur_pos[1]]  # type of square we're on
            random_num = np.random.random()  # [0,1)
            if cur_type == 1:  # endpoint
                # step any direction
                d_pos = np.asarray(moves[np.random.randint(0, 4)])
                final_pos = cur_pos + d_pos
                walker.add_pos(final_pos)
            elif cur_type == 0:  # matrix cell
                # generate candidate position
                d_pos = np.asarray(moves[np.random.randint(0, 4)])
                next_pos = cur_pos + d_pos
                next_type = grid.tube_check_bd_vol[next_pos[0], next_pos[1]]
                if next_type == -1:  # inside tube
                    if random_num > prob_m_cn:
                        final_pos = cur_pos
                        walker.add_pos(final_pos)
                    else:
                        final_pos = next_pos
                        walker.add_pos(final_pos)
                else:  # a normal step
                    final_pos = next_pos
                    walker.add_pos(final_pos)
            elif cur_type == -1:  # CNT cell
                # find index of current tube walker is in
                cur_tube_index = grid.tube_check_index[cur_pos[0], cur_pos[1]]
                cubes_in_tube = len(grid.tube_squares[cur_tube_index - 1])
                # -1 ABOVE AND BELOW BECAUSE OF +1 OFFSET IN CREATION TO AVOID ZERO INDEX
                next_pos = grid.tube_squares[cur_tube_index - 1][np.random.randint(0, cubes_in_tube)]
                # move to another random point in the tube
                d_pos = np.asarray(moves[np.random.randint(0, 4)])
                candidate_pos = np.asarray(next_pos) + d_pos
                # see where the candidate pos is
                candidate_loc = grid.tube_check_bd_vol[candidate_pos[0], candidate_pos[1]]
                # if candidate is in tube or on endpoint or random < kapitza move, else stay
                if (candidate_loc == -1) or (candidate_loc == 1) or (random_num < prob_m_cn):
                    final_pos = candidate_pos
                    walker.add_pos(final_pos)
                else:
                    final_pos = next_pos
                    walker.add_pos(final_pos)
        else:  # standard hopping through tubes
            cur_type = grid.tube_check_bd[cur_pos[0], cur_pos[1]]  # type of square we're on
            if cur_type == 0:  # matrix cell
                d_pos = np.asarray(moves[np.random.randint(0, 4)])
                final_pos = cur_pos + d_pos
                walker.add_pos(final_pos)
            elif cur_type == 1:  # endpoint
                jump_moves = [[0, 1], [1, 0], [0, -1], [-1, 0], [0, 1], [1, 0], [0, -1], [-1, 0]]
                # 8 possible jump positions for now if at tube end, 4 at left (first 4) and 4 at right (second 4)
                choice = np.random.randint(0, 8)
                d_pos = np.asarray(jump_moves[choice])
                # coord on left tube end jumps to right end
                tube_check_val_l = grid.tube_check_l[cur_pos[0], cur_pos[1]]
                if tube_check_val_l > 0:
                    if choice <= 3:  # stay at left end
                        final_pos = cur_pos + d_pos
                        walker.add_pos(final_pos)
                    else:  # jump across tube to right end
                        # -1 BELOW BECAUSE OF +1 OFFSET IN CREATION TO AVOID ZERO INDEX
                        final_pos = np.asarray(grid.tube_coords_r[tube_check_val_l - 1]) + np.asarray(d_pos)
                        walker.add_pos(final_pos)
                # coord on right tube end jumps to left end
                tube_check_val_r = grid.tube_check_r[cur_pos[0], cur_pos[1]]
                if tube_check_val_r > 0:
                    if choice <= 3:  # stay at right end
                        final_pos = cur_pos + d_pos
                        walker.add_pos(final_pos)
                    else:  # jump across tube to left end
                        # -1 BELOW BECAUSE OF +1 OFFSET IN CREATION TO AVOID ZERO INDEX
                        final_pos = np.asarray(grid.tube_coords_l[tube_check_val_r - 1]) + np.asarray(d_pos)
                        walker.add_pos(final_pos)
        if cur_type == 10:  # reflective boundary
            d_pos = np.asarray(moves[np.random.randint(0, 4)])
            next_pos = cur_pos + d_pos
            if (next_pos[0] > grid.size) or (next_pos[0] < 0):
                final_pos = cur_pos
                walker.add_pos(final_pos)
            else:
                final_pos = next_pos
                walker.add_pos(final_pos)
        elif cur_type == 20:  # periodic boundary
            d_pos = np.asarray(moves[np.random.randint(0, 4)])
            next_pos = cur_pos + d_pos
            if (next_pos[1] > grid.size) or (next_pos[1] < 0):
                final_pos = [next_pos[0], 1]
                walker.add_pos(final_pos)
            else:
                final_pos = next_pos
                walker.add_pos(final_pos)
        elif cur_type == 30:  # corner
            d_pos = np.asarray(moves[np.random.randint(0, 4)])
            next_pos = cur_pos + d_pos
            if (next_pos[0] > grid.size) or (next_pos[0] < 0) or (next_pos[1] > grid.size) or (next_pos[1] < 0):
                final_pos = cur_pos
                walker.add_pos(final_pos)
            else:
                final_pos = next_pos
                walker.add_pos(final_pos)
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
    # NOT CURRENTLY USED
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
    # NOT CURRENTLY USED
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
