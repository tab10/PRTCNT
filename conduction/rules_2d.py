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


from builtins import range
import logging
from conduction import *
import numpy as np


def runrandomwalk_2d_onlat(grid, timesteps, temp, kapitza, prob_m_cn, object):
    # Start the random walk, for one walker
    walker = creation.Walker2D_onlat(grid.size, temp)
    for i in range(1, timesteps + 1):
        walker = apply_moves_2d(walker, kapitza, grid, prob_m_cn, object)
    return walker


def apply_moves_2d(walker, kapitza, grid, prob_m_cn, object):
    '''Maybe given walker object or coordinates directly, check object'''
    moves = [[0, 1], [1, 0], [0, -1], [-1, 0]]
    # check where we are
    if object:
        cur_pos = np.asarray(walker.pos[-1])
    else:
        cur_pos = np.asarray(walker[-1])
    # DEBUG
    # Having tube radius doesn't matter if kapitza is off, apart from excluded volume
    if kapitza:
        cur_type = grid.tube_check_bd_vol[cur_pos[0], cur_pos[1]]  # type of square we're on
        random_num = np.random.random()  # [0,1)
        if cur_type == 1:  # endpoint
            # step any direction
            d_pos = np.asarray(moves[np.random.randint(0, 4)])
            final_pos = cur_pos + d_pos
            # walker.add_pos(final_pos)
        elif cur_type == 0:  # matrix cell
            # generate candidate position
            d_pos = np.asarray(moves[np.random.randint(0, 4)])
            candidate_pos = cur_pos + d_pos
            candidate_type = grid.tube_check_bd_vol[candidate_pos[0], candidate_pos[1]]
            if candidate_type == -1:  # inside tube
                if random_num > prob_m_cn:
                    final_pos = cur_pos
                    # walker.add_pos(final_pos)
                else:
                    final_pos = candidate_pos
                    # walker.add_pos(final_pos)
            else:  # a normal step
                final_pos = candidate_pos
                # walker.add_pos(final_pos)
        elif cur_type == -1:  # CNT cell
            # find index of current tube walker is in
            cur_index = grid.tube_check_index[cur_pos[0], cur_pos[1]]
            cubes_in_tube = len(grid.tube_squares[cur_index - 1])
            # -1 ABOVE AND BELOW BECAUSE OF +1 OFFSET IN CREATION TO AVOID ZERO INDEX
            candidate_pos = grid.tube_squares[cur_index - 1][np.random.randint(0, cubes_in_tube)]
            # move to another random point in the tube
            d_pos = np.asarray(moves[np.random.randint(0, 4)])
            candidate_pos = np.asarray(candidate_pos) + d_pos
            # see where the candidate pos is
            candidate_loc = grid.tube_check_bd_vol[candidate_pos[0], candidate_pos[1]]
            # if candidate is in tube or on endpoint or random < kapitza move, else stay
            if (candidate_loc == -1) or (candidate_loc == 1) or (random_num < prob_m_cn):
                final_pos = candidate_pos
                # walker.add_pos(final_pos)
            else:
                final_pos = candidate_pos
                # walker.add_pos(final_pos)
    else:  # standard hopping through tubes
        cur_type = grid.tube_check_bd[cur_pos[0], cur_pos[1]]  # type of square we're on
        if cur_type == 0:  # matrix cell
            d_pos = np.asarray(moves[np.random.randint(0, 4)])
            final_pos = cur_pos + d_pos
            # walker.add_pos(final_pos)
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
                    # walker.add_pos(final_pos)
                else:  # jump across tube to right end
                    # -1 BELOW BECAUSE OF +1 OFFSET IN CREATION TO AVOID ZERO INDEX
                    final_pos = np.asarray(grid.tube_coords_r[tube_check_val_l - 1]) + np.asarray(d_pos)
                    # walker.add_pos(final_pos)
            # coord on right tube end jumps to left end
            tube_check_val_r = grid.tube_check_r[cur_pos[0], cur_pos[1]]
            if tube_check_val_r > 0:
                if choice <= 3:  # stay at right end
                    final_pos = cur_pos + d_pos
                    # walker.add_pos(final_pos)
                else:  # jump across tube to left end
                    # -1 BELOW BECAUSE OF +1 OFFSET IN CREATION TO AVOID ZERO INDEX
                    final_pos = np.asarray(grid.tube_coords_l[tube_check_val_r - 1]) + np.asarray(d_pos)
                    # walker.add_pos(final_pos)
    if cur_type == 10:  # reflective boundary
        d_pos = np.asarray(moves[np.random.randint(0, 4)])
        candidate_pos = cur_pos + d_pos
        if (candidate_pos[0] > grid.size) or (candidate_pos[0] < 0):
            final_pos = cur_pos
            # walker.add_pos(final_pos)
        else:
            final_pos = candidate_pos
            # walker.add_pos(final_pos)
    elif cur_type == 20:  # periodic boundary
        d_pos = np.asarray(moves[np.random.randint(0, 4)])
        candidate_pos = cur_pos + d_pos
        if candidate_pos[1] > grid.size:
            final_pos = [candidate_pos[0], 1]
            # walker.add_pos(final_pos)
        elif candidate_pos[1] < 0:
            final_pos = [candidate_pos[0], grid.size - 1]
            # walker.add_pos(final_pos)
        else:
            final_pos = candidate_pos
            # walker.add_pos(final_pos)
    elif cur_type == 30:  # corner
        d_pos = np.asarray(moves[np.random.randint(0, 4)])
        candidate_pos = cur_pos + d_pos
        if (candidate_pos[0] > grid.size) or (candidate_pos[0] < 0) or (candidate_pos[1] > grid.size) or (
            candidate_pos[1] < 0):
            final_pos = cur_pos
            # walker.add_pos(final_pos)
        else:
            final_pos = candidate_pos
    if object:
        walker.add_pos(final_pos)
    else:
        walker.append(list(final_pos))
    return walker
