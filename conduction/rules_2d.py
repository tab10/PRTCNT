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


import logging
from conduction import *
import numpy as np


def kill(message="Invalid random walk rule. Check rules."):
    """Stop the program"""
    logging.error(message)
    raise SystemExit


def runrandomwalk_2d_onlat(grid, timesteps, temp, kapitza, prob_m_cn, bound, rules_test):
    # Start the random walk, for one walker
    walker = creation_2d.Walker2D_onlat(grid.size, temp, rules_test)
    inside_cnt = False
    for i in range(1, timesteps + 1):
        walker, inside_cnt = apply_moves_2d(walker, kapitza, grid, prob_m_cn, inside_cnt, bound)
    return walker


def generate_novol_choices_2d(grid, moves_2d, cur_pos, tube_index, kapitza, return_pos=True):
    """Check all directions (4) and remove those that are still part of the SAME CNT volume.
    If not kapitza given, will remove ALL volumes. tube_index not used."""
    possible_moves = []
    possible_locs = []
    for i in range(len(moves_2d)):
        candidate_temp = cur_pos + np.asarray(moves_2d[i])
        candidate_temp_index = grid.tube_check_index[candidate_temp[0], candidate_temp[1]] - 1
        candidate_temp_type = grid.tube_check_bd_vol[candidate_temp[0], candidate_temp[1]]  # type of square we're on
        if kapitza:
            if not ((candidate_temp_type == -1) and (tube_index == candidate_temp_index)):
                # candidate must not be CNT volume from the SAME tube
                # cannot move back into same tube
                possible_moves.append(list(moves_2d[i]))
                possible_locs.append(list(candidate_temp))
        elif not kapitza:  # candidate must not be ANY CNT volume (tunneling)
            if not candidate_temp_index == -1:
                # cannot move to any volume points
                possible_moves.append(list(moves_2d[i]))
                possible_locs.append(list(candidate_temp))
        else:
            kill('Check kapitza value.')
    num_possible_moves = len(possible_moves)
    num_possible_locs = len(possible_locs)
    if not return_pos:
        return possible_moves, num_possible_moves
    else:
        return possible_locs, num_possible_locs


def apply_bd_cond_2d(grid, moves_2d, cur_pos, bound):
    choices = generate_bd_choices_2d(grid, cur_pos, moves_2d, bound)
    # pick random choice
    new_pos = np.asarray(choices[np.random.randint(0, len(choices))])
    final_pos = new_pos
    return final_pos


def generate_bd_choices_2d(grid, cur_pos, moves_2d, bound):
    choices = cur_pos + moves_2d
    # get edge locations
    idx_remove = []
    min_val = 0
    max_val = grid.size  # not + 1 as in setup as we can walk on 0 or 100
    for i in range(len(choices)):
        temp = choices[i]
        for j in range(len(temp)):
            type_bound = bound[j]
            coord_temp = temp[j]
            if type_bound == 10:  # reflective
                if coord_temp < min_val:
                    # stay
                    coord_temp = min_val
                elif coord_temp > max_val:
                    coord_temp = max_val
                    # coord_temp = None
                    # idx_remove.append(i)
            elif type_bound == 20:  # periodic
                coord_temp = (coord_temp % max_val)
            temp[j] = coord_temp
        choices[i] = list(temp)
        # print(choices[i])
    # for index in sorted(idx_remove, reverse=True):
    #    del choices[index]
    return choices


def apply_moves_2d(walker, kapitza, grid, prob_m_cn, inside_cnt, bound):
    # having the moves as lists is OK since numpy arrays + list is the standard + we want
    inert_vol = grid.inert_vol
    moves_2d = [[0, 1], [1, 0], [0, -1], [-1, 0]]
    jump_moves_2d = [[0, 1], [1, 0], [0, -1], [-1, 0], [0, 1], [1, 0], [0, -1], [-1, 0]]
    moves_2d_diag = [[0, 1], [1, 0], [0, -1], [-1, 0], [1, 1], [1, -1], [-1, 1], [-1, -1]]
    jump_moves_2d_diag = [[0, 1], [1, 0], [0, -1], [-1, 0], [1, 1], [1, -1], [-1, 1], [-1, -1],
                          [0, 1], [1, 0], [0, -1], [-1, 0], [1, 1], [1, -1], [-1, 1], [-1, -1]]
    cur_pos = np.asarray(walker.pos[-1])
    if kapitza:
        cur_type = grid.tube_check_bd_vol[cur_pos[0], cur_pos[1]]  # type of square we're on
        cur_index = grid.tube_check_index[cur_pos[0], cur_pos[1]] - 1  # index>0 of CNT (or 0 for not one)
        if cur_type == 1:  # CNT end
            final_pos, inside_cnt = kapitza_cntend(grid, moves_2d, kapitza, cur_pos, cur_index, prob_m_cn,
                                                   inside_cnt)
            walker.add_pos(final_pos)
        elif cur_type == 0:  # matrix cell
            final_pos, inside_cnt = kapitza_matrix(grid, moves_2d, cur_pos, cur_index, prob_m_cn)
            walker.add_pos(final_pos)
        elif cur_type == -1:  # CNT volume
            final_pos, inside_cnt = kapitza_cntvol(grid, moves_2d, kapitza, cur_pos, cur_index, prob_m_cn,
                                                   inside_cnt)
            walker.add_pos(final_pos)
        elif cur_type == -1000:  # boundary
            final_pos = apply_bd_cond_2d(grid, moves_2d, cur_pos, bound)
            walker.add_pos(final_pos)
        else:
            exit()
    elif not kapitza:  # tunneling with/without inert volume models
        if inert_vol:
            cur_type = grid.tube_check_bd_vol[cur_pos[0], cur_pos[1]]  # type of square we're on
            cur_index = grid.tube_check_index[cur_pos[0], cur_pos[1]] - 1  # index>0 of CNT (or 0 for not one)
        else:
            cur_type = grid.tube_check_bd[cur_pos[0], cur_pos[1]]  # type of square we're on
            cur_index = None
        if cur_type == 0:  # matrix cell
            final_pos = tunneling_matrix(grid, moves_2d, cur_pos, cur_index, inert_vol)
            walker.add_pos(final_pos)
        elif cur_type == 1:  # endpoint
            final_pos = tunneling_cntend(grid, jump_moves_2d, cur_pos, inert_vol)
            walker.add_pos(final_pos)
        elif cur_type == -1000:  # boundary
            final_pos = apply_bd_cond_2d(grid, moves_2d, cur_pos, bound)
            walker.add_pos(final_pos)
        else:
            exit()
    return walker, inside_cnt


def kapitza_cntend(grid, moves_2d, kapitza, cur_pos, cur_index, prob_m_cn, inside_cnt):
    choice = np.random.randint(0, len(moves_2d))
    d_pos = np.asarray(moves_2d[choice])
    candidate_pos = cur_pos + d_pos
    candidate_type = grid.tube_check_bd_vol[candidate_pos[0], candidate_pos[1]]
    candidate_index = grid.tube_check_index[candidate_pos[0], candidate_pos[1]] - 1
    if candidate_type == -1:  # CNT volume
        if candidate_index == cur_index:  # Same tube, send it back in
            # move to random volume/endpoint within same CNT
            final_pos = np.asarray(
                grid.tube_squares[cur_index][np.random.randint(0, len(grid.tube_squares[cur_index]))])
            inside_cnt = True
        else:  # NOT same tube, check Kapitza
            random_num = np.random.random()  # [0.0, 1.0)
            enter = (random_num < prob_m_cn)
            if enter:  # move to random volume/endpoint within new CNT
                final_pos = np.asarray(
                    grid.tube_squares[candidate_index][np.random.randint(0, len(grid.tube_squares[candidate_index]))])
                inside_cnt = True
            else:  # stay put
                final_pos = cur_pos
    else:  # matrix, CNT end, or boundary (walk off) NO MORE TUNNELING AS OF 5_9_17 TAB
        final_pos = candidate_pos
        inside_cnt = False
    return final_pos, inside_cnt


def kapitza_matrix(grid, moves_2d, cur_pos, cur_index, prob_m_cn):
    # generate candidate position
    d_pos = np.asarray(moves_2d[np.random.randint(0, len(moves_2d))])
    candidate_pos = cur_pos + d_pos
    candidate_type = grid.tube_check_bd_vol[candidate_pos[0], candidate_pos[1]]
    candidate_idx = grid.tube_check_index[candidate_pos[0], candidate_pos[1]] - 1
    if candidate_type == -1:  # CNT volume
        random_num = np.random.random()  # [0.0, 1.0)
        kap_enter = (random_num < prob_m_cn)
        if kap_enter:
            # move to random volume/endpoint within same CNT
            final_pos = np.asarray(
                grid.tube_squares[candidate_idx][np.random.randint(0, len(grid.tube_squares[candidate_idx]))])
            inside_cnt = True
        else:
            ### SIT
            final_pos = cur_pos
            inside_cnt = False
    elif candidate_type == 1:  # CNT end
        # move to random point within tube
        final_pos = np.asarray(
            grid.tube_squares[cur_index][np.random.randint(0, len(grid.tube_squares[cur_index]))])
        inside_cnt = True
    else:  # CNT boundary or matrix
        # move there
        final_pos = candidate_pos
        inside_cnt = False
    return final_pos, inside_cnt


def kapitza_cntvol(grid, moves_2d, kapitza, cur_pos, cur_index, prob_m_cn, inside_cnt):
    d_pos = np.asarray(moves_2d[np.random.randint(0, len(moves_2d))])
    candidate_pos = cur_pos + d_pos
    candidate_type = grid.tube_check_bd_vol[candidate_pos[0], candidate_pos[1]]
    candidate_idx = grid.tube_check_index[candidate_pos[0], candidate_pos[1]] - 1
    if (candidate_type == 0) or (candidate_type == -1000):  # matrix or boundary
        # check if the walker is inside or outside of a CNT
        if inside_cnt:  # wants to leave
            random_num = np.random.random()  # [0.0, 1.0)
            stay = (random_num > prob_m_cn)
            if stay:  # move to random volume/endpoint within same CNT
                final_pos = np.asarray(
                    grid.tube_squares[cur_index][np.random.randint(0, len(grid.tube_squares[cur_index]))])
                inside_cnt = True
            else:  # walk outside tube
                final_pos = np.asarray(candidate_pos)
                inside_cnt = False
        else:  # on CNT volume NOT in tube, can happen if walker is spawned on a CNT
            # move to a safe area
            final_pos = np.asarray(candidate_pos)
            inside_cnt = False
    elif (candidate_type == -1) or (candidate_type == 1):  # CNT volume or end
        if candidate_idx == cur_index:  # want to go to CNT volume in same tube
            final_pos = np.asarray(
                grid.tube_squares[cur_index][np.random.randint(0, len(grid.tube_squares[cur_index]))])
            inside_cnt = True
        else:  # wants to enter a new tube
            random_num = np.random.random()  # [0.0, 1.0)
            stay = (random_num > prob_m_cn)
            if stay:  # move to random volume/endpoint within same CNT
                final_pos = np.asarray(
                    grid.tube_squares[cur_index][np.random.randint(0, len(grid.tube_squares[cur_index]))])
                inside_cnt = True
            else:  # exit to new
                final_pos = np.asarray(
                    grid.tube_squares[candidate_idx][np.random.randint(0, len(grid.tube_squares[candidate_idx]))])
                inside_cnt = True
    else:
        exit()
    return final_pos, inside_cnt


def tunneling_matrix(grid, moves_2d, cur_pos, cur_index, inert_vol):
    if inert_vol:  # checking that no CNT volume is a possibility
        possible_locs, num_possible_locs = generate_novol_choices_2d(grid, moves_2d, cur_pos, cur_index, False,
                                                                 return_pos=True)  # turning off Kapitza for tunneling
        final_pos = np.asarray(possible_locs[np.random.randint(0, num_possible_locs)])
    else:
        final_pos = np.asarray(moves_2d[np.random.randint(0, len(moves_2d))])
    return final_pos


def tunneling_cntend(grid, jump_moves_2d, cur_pos, inert_vol):
    # no CNT volume, so this rule remains unchanged from the originals (3/30/17 TB)
    # walk off either end, 12 3D choices always
    choice = np.random.randint(0, 8)
    d_pos = np.asarray(jump_moves_2d[choice])
    # coord on left tube end jumps to right end
    tube_check_val_l = grid.tube_check_l[cur_pos[0], cur_pos[1]]
    # coord on right tube end jumps to left end
    tube_check_val_r = grid.tube_check_r[cur_pos[0], cur_pos[1]]
    if (tube_check_val_l > 0) and (tube_check_val_r == 0):
        # check that pixel cannot be a left and right endpoint
        if choice <= 3:  # stay at left end
            final_pos = cur_pos + d_pos
        else:  # jump across tube to right end
            # -1 BELOW BECAUSE OF +1 OFFSET IN CREATION TO AVOID ZERO INDEX
            final_pos = np.asarray(grid.tube_coords_r[tube_check_val_l - 1]) + np.asarray(d_pos)
    elif (tube_check_val_r > 0) and (tube_check_val_l == 0):
        # check that pixel cannot be a left and right endpoint
        if choice <= 3:  # stay at right end
            final_pos = cur_pos + d_pos
        else:  # jump across tube to left end
            # -1 BELOW BECAUSE OF +1 OFFSET IN CREATION TO AVOID ZERO INDEX
            final_pos = np.asarray(grid.tube_coords_l[tube_check_val_r - 1]) + np.asarray(d_pos)
    else:
        kill()
    return final_pos
