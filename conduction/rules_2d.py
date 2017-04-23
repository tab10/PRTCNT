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
            final_pos, inside_cnt = kapitza_matrix(grid, moves_2d, kapitza, cur_pos, cur_index, prob_m_cn,
                                                   inside_cnt)
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


def kapitza_cntend(grid, jump_moves_2d_diag, kapitza, cur_pos, cur_index, prob_m_cn, inside_cnt):
    l_moves = jump_moves_2d_diag
    r_moves = jump_moves_2d_diag
    # we need to check if the endpoint has the possibility to jump back in, or we will have to add it
    p0 = np.asarray(grid.tube_squares[cur_index][0])
    p1 = np.asarray(grid.tube_squares[cur_index][1])
    p2 = np.asarray(grid.tube_squares[cur_index][-2])
    p3 = np.asarray(grid.tube_squares[cur_index][-1])
    checkl = np.abs(p1 - p0)
    checkr = np.abs(p3 - p2)
    if np.all(checkl == np.asarray([1, 1])):
        l_moves.append(p1 - p0)
    if np.all(checkr == np.asarray([1, 1])):
        r_moves.append(p3 - p2)
    # generate candidate position
    ends = ['l', 'r']
    c_end = np.random.choice(ends)
    if c_end == 'l':
        choice = np.random.randint(0, len(l_moves))
        d_pos = np.asarray(l_moves[choice])
    elif c_end == 'r':
        choice = np.random.randint(0, len(r_moves))
        d_pos = np.asarray(r_moves[choice])
    candidate_pos = cur_pos + d_pos
    # d_pos = np.asarray(jump_moves_2d_diag[choice])
    # # coord on left tube end jumps to right end
    # tube_check_val_l = grid.tube_check_l[cur_pos[0], cur_pos[1]]
    # # coord on right tube end jumps to left end
    # tube_check_val_r = grid.tube_check_r[cur_pos[0], cur_pos[1]]
    # if (tube_check_val_l > 0) and (tube_check_val_r == 0):
    #     # check that pixel cannot be a left and right endpoint
    #     if choice <= split_pt:  # stay at left end
    #         candidate_pos = cur_pos + d_pos
    #     else:  # jump across tube to right end
    #         # -1 BELOW BECAUSE OF +1 OFFSET IN CREATION TO AVOID ZERO INDEX
    #         candidate_pos = np.asarray(grid.tube_coords_r[tube_check_val_l - 1]) + np.asarray(d_pos)
    # elif (tube_check_val_r > 0) and (tube_check_val_l == 0):
    #     # check that pixel cannot be a left and right endpoint
    #     if choice <= split_pt:  # stay at right end
    #         candidate_pos = cur_pos + d_pos
    #     else:  # jump across tube to left end
    #         # -1 BELOW BECAUSE OF +1 OFFSET IN CREATION TO AVOID ZERO INDEX
    #         candidate_pos = np.asarray(grid.tube_coords_l[tube_check_val_r - 1]) + np.asarray(d_pos)
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
            not_enter = (random_num > prob_m_cn)
            if enter:  # move to random volume/endpoint within new CNT
                final_pos = np.asarray(
                    grid.tube_squares[candidate_index][np.random.randint(0, len(grid.tube_squares[candidate_index]))])
                inside_cnt = True
            else:  # stay put
                final_pos = cur_pos
    else:  # matrix, CNT end, or boundary (walk off)
        final_pos = candidate_pos
        inside_cnt = False
    return final_pos, inside_cnt


    # if inside_cnt:
    #     random_num = np.random.random()  # [0.0, 1.0)
    #     stay = (random_num > prob_m_cn)
    #     leave = (random_num < prob_m_cn)
    #     if stay:  # move to random volume/endpoint within same CNT
    #         final_pos = np.asarray(
    #             grid.tube_squares[cur_index][np.random.randint(0, len(grid.tube_squares[cur_index]))])
    #         inside_cnt = True
    #     else:  # sit there or tunnel
    #         if not jump:
    #             final_pos = np.asarray(cur_end)
    #         else:
    #             final_pos = np.asarray(jump_end)
    #         inside_cnt = True
    #
    # if list(cur_pos) == grid.tube_squares[cur_index][0]:
    #     print('L')
    #     cur_end = grid.tube_squares[cur_index][0]
    #     jump_end = grid.tube_squares[cur_index][-1]
    # elif list(cur_pos) == grid.tube_squares[cur_index][-1]:
    #     print('R')
    #     cur_end = grid.tube_squares[cur_index][-1]
    #     jump_end = grid.tube_squares[cur_index][0]
    # else:
    #     exit()
    # if choice <= 7:
    #     candidate_pos = np.asarray(cur_end) + d_pos
    # else:
    #     candidate_pos = np.asarray(jump_end) + d_pos


    #     choice = np.random.randint(0, 8)
    #     d_pos = np.asarray(jump_moves_2d[choice])
    #     # coord on left tube end jumps to right end
    #     tube_check_val_l = grid.tube_check_l[cur_pos[0], cur_pos[1]]
    #     # coord on right tube end jumps to left end
    #     tube_check_val_r = grid.tube_check_r[cur_pos[0], cur_pos[1]]
    #     if (tube_check_val_l > 0) and (tube_check_val_r == 0):
    #         # check that pixel cannot be a left and right endpoint
    #         if choice <= 3:  # stay at left end
    #             final_pos = cur_pos + d_pos
    #         else:  # jump across tube to right end
    #             # -1 BELOW BECAUSE OF +1 OFFSET IN CREATION TO AVOID ZERO INDEX
    #             final_pos = np.asarray(grid.tube_coords_r[tube_check_val_l - 1]) + np.asarray(d_pos)
    #     elif (tube_check_val_r > 0) and (tube_check_val_l == 0):
    #         # check that pixel cannot be a left and right endpoint
    #         if choice <= 3:  # stay at right end
    #             final_pos = cur_pos + d_pos
    #         else:  # jump across tube to left end
    #             # -1 BELOW BECAUSE OF +1 OFFSET IN CREATION TO AVOID ZERO INDEX
    #             final_pos = np.asarray(grid.tube_coords_l[tube_check_val_r - 1]) + np.asarray(d_pos)
    #
    #     #  remove current spot from choices
    #     new_choices = []
    #     coord_del = [list(cur_pos)]
    #     for x in grid.tube_squares[cur_index]:
    #         if x not in coord_del:
    #             new_choices.append(x)
    #     # -1 because of above statement
    #     num_new_choices = len(new_choices)
    #     final_pos = np.asarray(new_choices[np.random.randint(0, num_new_choices)])
    #     inside_cnt = True
    # else:  # exit CNT (just go there easy!)
    #     # coord on left tube end jumps to right end
    #     tube_check_val_l = grid.tube_check_l[cur_pos[0], cur_pos[1]]
    #     # coord on right tube end jumps to left end
    #     tube_check_val_r = grid.tube_check_r[cur_pos[0], cur_pos[1]]
    #     if (tube_check_val_l > 0) and (tube_check_val_r == 0):
    #         # check that pixel cannot be a left and right endpoint
    #         if choice <= 8:  # stay at left end
    #             final_pos = cur_pos + d_pos
    #         else:  # jump across tube to right end
    #             # -1 BELOW BECAUSE OF +1 OFFSET IN CREATION TO AVOID ZERO INDEX
    #             final_pos = np.asarray(grid.tube_coords_r[tube_check_val_l - 1]) + np.asarray(d_pos)
    #     elif (tube_check_val_r > 0) and (tube_check_val_l == 0):
    #         # check that pixel cannot be a left and right endpoint
    #         if choice <= 8:  # stay at right end
    #             final_pos = cur_pos + d_pos
    #         else:  # jump across tube to left end
    #             # -1 BELOW BECAUSE OF +1 OFFSET IN CREATION TO AVOID ZERO INDEX
    #             final_pos = np.asarray(grid.tube_coords_l[tube_check_val_r - 1]) + np.asarray(d_pos)
    #     else:
    #         kill()
    #     inside_cnt = False
    # return final_pos, inside_cnt

    # # generate candidate position
    # d_pos = np.asarray(moves_2d[np.random.randint(0, 4)])
    # candidate_pos = cur_pos + d_pos
    # candidate_type = grid.tube_check_bd_vol[candidate_pos[0], candidate_pos[1]]
    #
    # choices = cur_pos + moves_2d
    # probs = [2.0 / 8.0, 6.0 / 8.0]  # detailed balance
    # d_b = ['stay_enter', 'leave_notenter']  # two possibilities
    # d_b_choice = np.random.choice(d_b, p=probs)
    # if d_b_choice == 'stay_enter':  # move to random volume/endpoint within same CNT
    #     #  remove current spot from choices
    #     new_choices = []
    #     coord_del = [list(cur_pos)]
    #     for x in grid.tube_squares[cur_index - 1]:
    #         if x not in coord_del:
    #             new_choices.append(x)
    #     # -1 because of above statement
    #     num_new_choices = len(new_choices)
    #     final_pos = np.asarray(new_choices[np.random.randint(0, num_new_choices)])
    #     inside_cnt = True
    # elif d_b_choice == 'leave_notenter':  # exit CNT, on EITHER side randomly
    #     # it will walk off either of the two ends, checking that current CNT volume is not a possibility
    #     # collect position of both ends on current CNT
    #     end_1 = grid.tube_squares[cur_index - 1][0]
    #     end_2 = grid.tube_squares[cur_index - 1][-1]
    #     possible_locs_end1, num_possible_locs_end1 = generate_novol_choices_2d(grid, moves_2d, end_1, cur_index,
    #                                                                            kapitza,
    #                                                                            return_pos=True)
    #     possible_locs_end2, num_possible_locs_end2 = generate_novol_choices_2d(grid, moves_2d, end_2, cur_index,
    #                                                                            kapitza,
    #                                                                            return_pos=True)
    #     possible_locs = list(possible_locs_end1) + list(possible_locs_end2)
    #     num_possible_locs = num_possible_locs_end1 + num_possible_locs_end2
    #     final_pos = np.asarray(possible_locs[np.random.randint(0, num_possible_locs)])
    #     inside_cnt = False
    # else:
    #     kill()
    # return final_pos, inside_cnt


def kapitza_matrix(grid, moves_2d, kapitza, cur_pos, cur_index, prob_m_cn, inside_cnt):
    # generate candidate position
    d_pos = np.asarray(moves_2d[np.random.randint(0, len(moves_2d))])
    candidate_pos = cur_pos + d_pos
    candidate_type = grid.tube_check_bd_vol[candidate_pos[0], candidate_pos[1]]
    candidate_idx = grid.tube_check_index[candidate_pos[0], candidate_pos[1]] - 1
    if candidate_type == -1:  # CNT volume
        random_num = np.random.random()  # [0.0, 1.0)
        kap_enter = (random_num < prob_m_cn)
        if kap_enter:
            # move to random volume/endpoint within same CNT,
            # remove current spot from choices (NOT NOW)
            # coord_del = [list(cur_pos)]
            # new_choices = []
            # for x in grid.tube_squares[cur_index]:
            #     if x not in coord_del:
            #         new_choices.append(x)
            # num_new_choices = len(new_choices)
            # final_pos = np.asarray(new_choices[np.random.randint(0, num_new_choices)])
            final_pos = np.asarray(
                grid.tube_squares[candidate_idx][np.random.randint(0, len(grid.tube_squares[candidate_idx]))])
            inside_cnt = True
        else:
            ### SIT
            final_pos = cur_pos
            # walk away, checking that current CNT volume is not a possibility
            # possible_locs, num_possible_locs = generate_novol_choices_2d(grid, moves_2d_diag, cur_pos, cur_index,
            # kapitza,
            # return_pos=True)
            # final_pos = np.asarray(possible_locs[np.random.randint(0, num_possible_locs)])
            inside_cnt = False
            # random_num = np.random.random()  # [0.0, 1.0)
            # if random_num > prob_m_cn:
            #     # walk away, checking that current CNT volume is not a possibility
            #     # inside_cnt SHOULD be false here always
            #     possible_locs, num_possible_locs = generate_novol_choices_2d(grid, moves_2d, cur_pos, cur_index, kapitza,
            #                                                                  return_pos=True)
            #     final_pos = np.asarray(possible_locs[np.random.randint(0, num_possible_locs)])
            #     inside_cnt = False
            # elif random_num < prob_m_cn:  # move to random volume/endpoint within CNT
            #     # get candidate CNT index
            #     candidate_index = grid.tube_check_index[candidate_pos[0], candidate_pos[1]] - 1
            #     #  DON'T remove candidate_pos from choices (since outside CNT here)
            #     new_choices = grid.tube_squares[candidate_index - 1]
            #     num_new_choices = len(new_choices)
            #     final_pos = np.asarray(new_choices[np.random.randint(0, num_new_choices)])
            #     inside_cnt = True
            # else:
            #    kill()
    else:  # CNT end, boundary, or matrix
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
            leave = (random_num < prob_m_cn)
            if stay:  # move to random volume/endpoint within same CNT
                final_pos = np.asarray(
                    grid.tube_squares[cur_index][np.random.randint(0, len(grid.tube_squares[cur_index]))])
                inside_cnt = True
            else:  # walk outside tube
                final_pos = np.asarray(candidate_pos)
                inside_cnt = False
        else:  # on CNT volume NOT in tube, can happen if walker is spawned on a CNT
            # use current pos. instead of candidate in this case
            final_pos = np.asarray(
                grid.tube_squares[cur_index][np.random.randint(0, len(grid.tube_squares[cur_index]))])
            inside_cnt = True
            # random_num = np.random.random()  # [0.0, 1.0)
            # enter = (random_num < prob_m_cn)
            # not_enter = (random_num > prob_m_cn)
            # if enter:  # move to random volume/endpoint within new CNT it was sitting on
            #     final_pos = np.asarray(
            #         grid.tube_squares[cur_index][np.random.randint(0, len(grid.tube_squares[cur_index]))])
            #     inside_cnt = True
            # else:  # move off of tube
            #     final_pos = candidate_pos
            #     inside_cnt = False
            # final_pos = candidate_pos
            # enter = (random_num < prob_m_cn)
            # not_enter = (random_num > prob_m_cn)
            # if enter:  # move to random volume/endpoint within same CNT
            #     final_pos = np.asarray(
            #         grid.tube_squares[candidate_idx][np.random.randint(0, len(grid.tube_squares[candidate_idx]))])
            #     inside_cnt = True
            # else:
            #     final_pos = np.asarray(candidate_pos)
            #     inside_cnt = False
    elif candidate_type == -1:  # CNT volume
        if candidate_idx == cur_index:  # want to go to CNT volume in same tube
            final_pos = np.asarray(
                grid.tube_squares[cur_index][np.random.randint(0, len(grid.tube_squares[cur_index]))])
            inside_cnt = True
        else:  # wants to enter a new tube
            random_num = np.random.random()  # [0.0, 1.0)
            stay = (random_num > prob_m_cn)
            leave = (random_num < prob_m_cn)
            if stay:  # move to random volume/endpoint within same CNT
                final_pos = np.asarray(
                    grid.tube_squares[cur_index][np.random.randint(0, len(grid.tube_squares[cur_index]))])
                inside_cnt = True
            else:  # exit to new
                final_pos = np.asarray(
                    grid.tube_squares[candidate_idx][np.random.randint(0, len(grid.tube_squares[candidate_idx]))])
                inside_cnt = True
    elif candidate_type == 1:  # CNT end
        if candidate_idx == cur_index:  # want to go to CNT end in same tube, go to random Vol or End in same tube
            final_pos = np.asarray(
                grid.tube_squares[cur_index][np.random.randint(0, len(grid.tube_squares[cur_index]))])
            inside_cnt = True
        else:  # wants to enter a new tube
            random_num = np.random.random()  # [0.0, 1.0)
            stay = (random_num > prob_m_cn)
            leave = (random_num < prob_m_cn)
            if stay:  # move to random volume/endpoint within same CNT
                final_pos = np.asarray(
                    grid.tube_squares[cur_index][np.random.randint(0, len(grid.tube_squares[cur_index]))])
                inside_cnt = True
            else:  # exit, moving to new CNT end
                final_pos = candidate_pos
                inside_cnt = False
    else:
        exit()
    return final_pos, inside_cnt
    #
    # if kap_stay_enter:
    #     # move to random volume/endpoint within same CNT, remove current spot from choices
    #     new_choices = []
    #     coord_del = [list(cur_pos)]
    #     for x in grid.tube_squares[cur_index - 1]:
    #         if x not in coord_del:
    #             new_choices.append(x)
    #     num_new_choices = len(new_choices)
    #     final_pos = np.asarray(new_choices[np.random.randint(0, num_new_choices)])
    #     final_pos = np.asarray(grid.tube_bds[cur_index][np.random.randint(0, len(grid.tube_bds[cur_index]))])
    # elif kap_leave_notenter:
    #     # walk away, checking that current CNT volume is not a possibility
    #     possible_locs, num_possible_locs = generate_novol_choices_2d(grid, moves_2d, cur_pos, cur_index, kapitza,
    #                                                                  return_pos=True)
    #     final_pos = np.asarray(possible_locs[np.random.randint(0, num_possible_locs)])
    # else:
    #     kill()
    # #
    # inside_cnt = not inside_cnt
    #
    #return final_pos, inside_cnt


# def kapitza_cntvol(grid, moves_2d_diag, kapitza, cur_pos, cur_index, prob_m_cn, inside_cnt):
#     p_cn_m = grid.p_cn_m[cur_pos[0], cur_pos[1]]
#     # generate candidate position
#     d_pos = np.asarray(moves_2d_diag[np.random.randint(0, len(moves_2d_diag))])
#     candidate_pos = cur_pos + d_pos
#     candidate_type = grid.tube_check_bd_vol[candidate_pos[0], candidate_pos[1]]
#     candidate_index = grid.tube_check_index[candidate_pos[0], candidate_pos[1]] - 1
#     #if (candidate_type == -1) and (candidate_index == cur_index):
#     ###
#     if inside_cnt == False:
#         random_num = np.random.random()  # [0.0, 1.0)
#         kap_enter = (random_num < prob_m_cn)
#         if kap_enter:
#             # move to random volume/endpoint within same CNT,
#             # remove current spot from choices (NOT NOW)
#             # coord_del = [list(cur_pos)]
#             # new_choices = []
#             # for x in grid.tube_squares[cur_index]:
#             #     if x not in coord_del:
#             #         new_choices.append(x)
#             # num_new_choices = len(new_choices)
#             #final_pos = np.asarray(new_choices[np.random.randint(0, num_new_choices)])
#             final_pos = np.asarray(grid.tube_squares[cur_index][np.random.randint(0, len(grid.tube_squares[cur_index]))])
#             inside_cnt = True
#         else:
#             ### SIT
#             final_pos = cur_pos
#             # walk away, checking that current CNT volume is not a possibility
#             #possible_locs, num_possible_locs = generate_novol_choices_2d(grid, moves_2d_diag, cur_pos, cur_index,
#                                                                          #kapitza,
#                                                                          #return_pos=True)
#             #final_pos = np.asarray(possible_locs[np.random.randint(0, num_possible_locs)])
#             inside_cnt = False
#     else:
#         random_num = np.random.random()  # [0.0, 1.0)
#         kap_exit = (random_num < prob_m_cn)
#         if kap_exit:
#             #cur_pos = final_pos
#             #cur_index = grid.tube_check_index[cur_pos[0], cur_pos[1]] - 1
#             #bd_choice = np.asarray(grid.tube_bds[cur_index][np.random.randint(0, len(grid.tube_bds[cur_index]))])
#             #final_pos = bd_choice
#             #inside_cnt = False
#             ###
#             # walk away, checking that current CNT volume is not a possibility
#             #possible_locs, num_possible_locs = generate_novol_choices_2d(grid, moves_2d_diag, cur_pos, cur_index,
#             #                                                             kapitza,
#             #                                                             return_pos=True)
#             #final_pos = np.asarray(possible_locs[np.random.randint(0, num_possible_locs)])
#             final_pos = np.asarray(
#                 grid.tube_bds[cur_index][np.random.randint(0, len(grid.tube_bds[cur_index]))])
#             inside_cnt = False
#         random_num = np.random.random()  # [0.0, 1.0)
#         kap_redistrib = (random_num > p_cn_m)
#         if kap_redistrib:
#             # move to random volume/endpoint within same CNT, remove current spot from choices
#             # coord_del = [list(cur_pos)]
#             # new_choices = []
#             # for x in grid.tube_squares[cur_index]:
#             #     if x not in coord_del:
#             #         new_choices.append(x)
#             # num_new_choices = len(new_choices)
#             # final_pos = np.asarray(new_choices[np.random.randint(0, num_new_choices)])
#             final_pos = np.asarray(grid.tube_squares[cur_index][np.random.randint(0, len(grid.tube_squares[cur_index]))])
#             inside_cnt = True
#         else:
#             ### SIT
#             #final_pos = cur_pos
#             #inside_cnt = True
#             final_pos = np.asarray(
#                 grid.tube_bds[cur_index][np.random.randint(0, len(grid.tube_bds[cur_index]))])
#             inside_cnt = False
#     return final_pos, inside_cnt
#     ###
#     # else:  # on a CNT volume DIFFERENT from current position
#     #     # check if the walker is inside or outside of a CNT
#     #     random_num = np.random.random()  # [0.0, 1.0)
#     #     if inside_cnt:
#     #         kap_stay_enter = (random_num > prob_m_cn)
#     #         kap_leave_notenter = (random_num < prob_m_cn)
#     #     else:
#     #         kap_stay_enter = (random_num < prob_m_cn)
#     #         kap_leave_notenter = (random_num > prob_m_cn)
#     #     if kap_stay_enter:
#     #         # move to random volume/endpoint within same CNT, remove current spot from choices
#     #         coord_del = [list(cur_pos)]
#     #         new_choices = []
#     #         for x in grid.tube_squares[cur_index]:
#     #             if x not in coord_del:
#     #                 new_choices.append(x)
#     #         num_new_choices = len(new_choices)
#     #         final_pos = np.asarray(new_choices[np.random.randint(0, num_new_choices)])
#     #         inside_cnt = True
#     #     elif kap_leave_notenter:
#     #         # walk away, checking that current CNT volume is not a possibility
#     #         possible_locs, num_possible_locs = generate_novol_choices_2d(grid, moves_2d_diag, cur_pos, cur_index,
#     #                                                                      kapitza,
#     #                                                                      return_pos=True)
#     #         final_pos = np.asarray(possible_locs[np.random.randint(0, num_possible_locs)])
#     #     else:
#     #         kill()
#     #return final_pos, inside_cnt


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
