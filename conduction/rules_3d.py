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


# from builtins import range
import logging
from conduction import *
import numpy as np


def kill(message="Invalid random walk rule. Check rules."):
    """Stop the program"""
    logging.error(message)
    raise SystemExit


def runrandomwalk_3d_onlat(grid, timesteps, temp, kapitza, prob_m_cn, bound, rules_test):
    # Start the random walk, for one walker
    walker = creation.Walker3D_onlat(grid.size, temp, rules_test)
    inside_cnt = False
    for i in range(1, timesteps + 1):
        walker, inside_cnt = apply_moves_3d(walker, kapitza, grid, prob_m_cn, inside_cnt, bound)
    return walker


def generate_novol_choices_3d(grid, moves_3d, cur_pos, tube_index, kapitza, return_pos=True):
    """Check all directions (6) and remove those that are still part of the SAME CNT volume.
    If not kapitza given, will remove ALL volumes. tube_index not used."""
    possible_moves = []
    possible_locs = []
    for i in range(len(moves_3d)):
        candidate_temp = cur_pos + np.asarray(moves_3d[i])
        candidate_temp_index = grid.tube_check_index[candidate_temp[0], candidate_temp[1], candidate_temp[2]]
        if kapitza:
            if (candidate_temp_index != -1) and (tube_index != candidate_temp_index):
                # candidate must not be CNT volume from the SAME tube
                # cannot move back into same tube
                possible_moves.append(moves_3d[i])
                possible_locs.append(candidate_temp)
        elif not kapitza:  # candidate must not be ANY CNT volume (tunneling)
            if candidate_temp_index != -1:
                # cannot move to any volume points
                possible_moves.append(moves_3d[i])
                possible_locs.append(candidate_temp)
        else:
            kill('Check kapitza value.')
    num_possible_moves = len(possible_moves)
    num_possible_locs = len(possible_locs)
    if not return_pos:
        return possible_moves, num_possible_moves
    else:
        return possible_locs, num_possible_locs


def apply_bd_cond_3d(grid, moves_3d, cur_pos, cur_type, bound):
    choices = generate_bc_choices_3d(grid, cur_pos, moves_3d, bound)
    # pick random choice
    new_pos = np.asarray(choices[np.random.randint(0, 6)])
    final_pos = new_pos
    # if rules_test:
    #     choices = generate_periodic_bc_choices_3d(grid, cur_pos, moves_3d)
    #     # pick random choice
    #     new_pos = np.asarray(choices[np.random.randint(0, 6)])
    #     final_pos = new_pos
    # if cur_type == 10:  # reflective boundary
    #     choices = generate_periodic_bc_choices_3d(grid, cur_pos, moves_3d)
    #     # pick random choice
    #     d_pos = np.asarray(moves_3d[np.random.randint(0, 6)])
    #     candidate_pos = cur_pos + d_pos
    #     if (candidate_pos[0] > grid.size) or (candidate_pos[0] < 0):
    #         final_pos = cur_pos
    #     else:
    #         final_pos = candidate_pos
    # elif cur_type == 20:  # periodic boundary
    #     d_pos = np.asarray(moves_3d[np.random.randint(0, 6)])
    #     candidate_pos = cur_pos + d_pos
    #     if candidate_pos[1] > grid.size:
    #         final_pos = [candidate_pos[0], 1, candidate_pos[2]]
    #     elif candidate_pos[1] < 0:
    #         final_pos = [candidate_pos[0], grid.size - 1, candidate_pos[2]]
    #     elif candidate_pos[2] > grid.size:
    #         final_pos = [candidate_pos[0], candidate_pos[1], 1]
    #     elif candidate_pos[2] < 0:
    #         final_pos = [candidate_pos[0], candidate_pos[1], grid.size - 1]
    #     else:
    #         final_pos = candidate_pos
    # elif cur_type == 30:  # corner
    #     d_pos = np.asarray(moves_3d[np.random.randint(0, 6)])
    #     candidate_pos = cur_pos + d_pos
    #     if (candidate_pos[0] > grid.size) or (candidate_pos[0] < 0) or (candidate_pos[1] > grid.size) or (candidate_pos[1] < 0) \
    #             or (candidate_pos[2] > grid.size) or (candidate_pos[2] < 0):
    #         final_pos = cur_pos
    #     else:
    #         final_pos = candidate_pos
    # else:
    #     kill()
    return final_pos


def generate_bc_choices_3d(grid, cur_pos, moves_3d, bound):
    choices = cur_pos + moves_3d
    # get edge locations
    min_val = 0
    max_val = grid.size  # not + 1 as in setup as we can walk on 0 or 100
    for i in range(len(choices)):
        temp = choices[i]
        for j in range(len(bound)):
            type_bound = bound[j]
            coord_temp = temp[j]
            if type_bound == 10:  # reflective
                if (coord_temp < min_val) or (coord_temp > max_val):
                    coord_temp = None
            elif type_bound == 20:  # periodic
                coord_temp = (coord_temp % max_val)
            temp[j] = coord_temp
        if None in temp:
            del choices[i]
        else:
            choices[i] = temp
    return choices


# def generate_periodic_bc_choices_3d(grid, cur_pos, moves_3d):
#     """Returns a list of move choices for a walker on a periodic boundary"""
#     choices = cur_pos + moves_3d
#     # get edge locations
#     min_val = 0
#     max_val = grid.size  # not + 1 as in setup as we can walk on 0 or 100
#     # check choices for crossover
#     for i in range(len(choices)):
#         temp = choices[i]
#         for j in range(len(temp)):
#             if temp[j] < min_val:
#                 temp[j] = max_val
#             elif temp[j] > max_val:
#                 temp[j] = min_val
#         choices[i] = temp
#     return choices


# def generate_reflective_bc_choices_3d(grid, cur_pos, moves_3d):
#     """Returns a list of move choices for a walker on a reflective boundary"""
#     choices = cur_pos + moves_3d
#     # get edge locations
#     min_val = 0
#     max_val = grid.size  # not + 1 as in setup as we can walk on 0 or 100
#     # we are guaranteed that that coordinate outside bd is reflective because we're told that :P
#     for i in range(len(choices)):
#         temp = choices[i]
#         for j in range(len(temp)):
#             if (temp[j] < min_val) or (temp[j] > max_val):
#                 temp[j] = None
#                 break
#         if None in temp:
#             del choices[i]
#         else:
#             choices[i] = temp
#     return choices


# def generate_corner_bc_choices_3d(grid, cur_pos, moves_3d):
#     """Returns a list of move choices for a walker on a corner boundary"""
#     choices = cur_pos + moves_3d
#     # get edge locations
#     min_val = 0
#     max_val = grid.size  # not + 1 as in setup as we can walk on 0 or 100
#     # one will be at max and the other at min, that defines the type
#     for i in range(len(choices)):
#         temp = choices[i]
#         # 3d here, X reflective Y,Z periodic
#         if (temp[0] < min_val) or (temp[0] > max_val):  # X reflective, delete choice
#             del choices[i]
#             continue
#         #
#         if (temp[1] < min_val) and (temp[2] < min_val):
#             temp[1] = max_val
#             temp[2] = max_val
#         elif (temp[1] > max_val) and (temp[2] > max_val):
#             temp[1] = min_val
#             temp[2] = min_val
#         elif (temp[1] > max_val) and (temp[2] < min_val):
#             temp[1] = min_val
#             temp[2] = max_val
#         elif (temp[1] < min_val) and (temp[2] > max_val):
#             temp[1] = max_val
#             temp[2] = min_val
#         choices[i] = temp
#     return choices


def apply_moves_3d(walker, kapitza, grid, prob_m_cn, inside_cnt, bound):
    # having the moves as lists is OK since numpy arrays + list is the standard + we want
    moves_3d = [[0, 0, 1], [0, 1, 0], [1, 0, 0], [0, 0, -1], [0, -1, 0], [-1, 0, 0]]
    jump_moves_3d = [[0, 0, 1], [0, 1, 0], [1, 0, 0], [0, 0, -1], [0, -1, 0], [-1, 0, 0],
                     [0, 0, 1], [0, 1, 0], [1, 0, 0], [0, 0, -1], [0, -1, 0], [-1, 0, 0]]
    cur_pos = np.asarray(walker.pos[-1])
    if kapitza:
        cur_type = grid.tube_check_bd_vol[cur_pos[0], cur_pos[1], cur_pos[2]]  # type of square we're on
        cur_index = grid.tube_check_index[cur_pos[0], cur_pos[1], cur_pos[2]]  # index>0 of CNT (or 0 for not one)
        if cur_type == 1:  # CNT end
            final_pos = kapitza_cntend(grid, moves_3d, kapitza, cur_pos, cur_index)
            walker.add_pos(final_pos)
        elif cur_type == 0:  # matrix cell
            final_pos, inside_cnt = kapitza_matrix(grid, moves_3d, kapitza, cur_pos, cur_index, prob_m_cn, inside_cnt)
            walker.add_pos(final_pos)
        elif cur_type == -1:  # CNT volume
            final_pos, inside_cnt = kapitza_cntvol(grid, moves_3d, kapitza, cur_pos, cur_index, prob_m_cn, inside_cnt)
            walker.add_pos(final_pos)
        else:  # boundary rules
            final_pos = apply_bd_cond_3d(grid, moves_3d, cur_pos, cur_type, bound)
            walker.add_pos(final_pos)
    elif not kapitza:  # tunneling with/without inert volume models
        # check if we have inert volume, boolean
        inert_vol = ('grid.tube_check_bd_vol' in locals()) or ('grid.tube_check_bd_vol' in globals())
        if inert_vol:
            cur_type = grid.tube_check_bd_vol[cur_pos[0], cur_pos[1], cur_pos[2]]  # type of square we're on
            cur_index = grid.tube_check_index[cur_pos[0], cur_pos[1], cur_pos[2]]  # index>0 of CNT (or 0 for not one)
        else:
            cur_type = grid.tube_check_bd[cur_pos[0], cur_pos[1], cur_pos[2]]  # type of square we're on
            cur_index = None
        if cur_type == 0:  # matrix cell
            final_pos = tunneling_matrix(grid, moves_3d, cur_pos, cur_index)
            walker.add_pos(final_pos)
        elif cur_type == 1:  # endpoint
            final_pos = tunneling_cntend(grid, jump_moves_3d, cur_pos)
            walker.add_pos(final_pos)
        else:  # boundary rules
            final_pos = apply_bd_cond_3d(grid, moves_3d, cur_pos, cur_type, bound)
            walker.add_pos(final_pos)
    return walker, inside_cnt


def kapitza_cntend(grid, moves_3d, kapitza, cur_pos, cur_index):
    probs = [2.0 / 12.0, 10.0 / 12.0]  # detailed balance
    d_b = ['stay_enter', 'leave_notenter']  # two possibilities
    d_b_choice = np.random.choice(d_b, p=probs)
    if d_b_choice == 'stay_enter':  # move to random volume/endpoint within same CNT
        #  remove current spot from choices
        new_choices = grid.tube_squares[cur_index - 1].pop(cur_pos)  # -1 because of above statement
        num_new_choices = len(new_choices)
        final_pos = np.asarray(new_choices[np.random.randint(0, num_new_choices)])
    elif d_b_choice == 'leave_notenter':  # exit CNT, on EITHER side randomly
        # it will walk off either of the two ends, checking that current CNT volume is not a possibility
        # collect position of both ends on current CNT
        end_1 = grid.tube_squares[cur_index - 1][0]
        end_2 = grid.tube_squares[cur_index - 1][-1]
        possible_locs_end1, num_possible_locs_end1 = generate_novol_choices_3d(grid, moves_3d, end_1, cur_index,
                                                                               kapitza,
                                                                               return_pos=True)
        possible_locs_end2, num_possible_locs_end2 = generate_novol_choices_3d(grid, moves_3d, end_2, cur_index,
                                                                               kapitza,
                                                                               return_pos=True)
        possible_locs = list(possible_locs_end1) + list(possible_locs_end2)
        num_possible_locs = num_possible_locs_end1 + num_possible_locs_end2
        final_pos = np.asarray(possible_locs[np.random.randint(0, num_possible_locs)])
    else:
        kill()
    return final_pos


def kapitza_matrix(grid, moves_3d, kapitza, cur_pos, cur_index, prob_m_cn, inside_cnt):
    # generate candidate position
    d_pos = np.asarray(moves_3d[np.random.randint(0, 6)])
    candidate_pos = cur_pos + d_pos
    candidate_type = grid.tube_check_bd_vol[candidate_pos[0], candidate_pos[1], candidate_pos[2]]
    if candidate_type == -1:  # CNT volume
        random_num = np.random.random()  # [0.0, 1.0)
        if random_num > prob_m_cn:
            # walk away, checking that current CNT volume is not a possibility
            # inside_cnt SHOULD be false here always
            possible_locs, num_possible_locs = generate_novol_choices_3d(grid, moves_3d, cur_pos, cur_index, kapitza,
                                                                         return_pos=True)
            final_pos = np.asarray(possible_locs[np.random.randint(0, num_possible_locs)])
        elif random_num < prob_m_cn:  # move to random volume/endpoint within CNT
            # get candidate CNT index
            candidate_index = grid.tube_check_index[
                candidate_pos[0], candidate_pos[1], candidate_pos[2]]
            #  DON'T remove candidate_pos from choices (since outside CNT here)
            new_choices = grid.tube_squares[candidate_index - 1]
            num_new_choices = len(new_choices)
            final_pos = np.asarray(new_choices[np.random.randint(0, num_new_choices)])
            inside_cnt = True
        else:
            kill()
    else:  # CNT end or boundary
        # move there
        final_pos = candidate_pos
    return final_pos, inside_cnt


def kapitza_cntvol(grid, moves_3d, kapitza, cur_pos, cur_index, prob_m_cn, inside_cnt):
    # this has been written to follow coding best practices, for in or out of CNT in the same IF statement
    random_num = np.random.random()  # [0.0, 1.0)
    # check if the walker is inside or outside of a CNT
    if inside_cnt:
        probs = [2.0 / 6.0, 4.0 / 6.0]  # detailed balance
        kap_stay_enter = (random_num > prob_m_cn)
        kap_leave_notenter = (random_num < prob_m_cn)
    else:
        probs = [1.0 / 6.0, 5.0 / 6.0]  # detailed balance
        kap_stay_enter = (random_num < prob_m_cn)
        kap_leave_notenter = (random_num > prob_m_cn)
    d_b = ['stay_enter', 'leave_notenter']  # two possibilities
    d_b_choice = np.random.choice(d_b, p=probs)
    if ((d_b_choice == 'stay_enter') and kap_stay_enter) \
            or ((d_b_choice == 'leave_notenter') and kap_stay_enter):
        # (Detailed balance stay) AND (Kapitza stay)
        # OR
        # (Detailed balance leave) AND (Kapitza stay)
        # move to random volume/endpoint within same CNT, remove current spot from choices
        new_choices = grid.tube_squares[cur_index - 1].pop(cur_pos)
        num_new_choices = len(new_choices)
        final_pos = np.asarray(new_choices[np.random.randint(0, num_new_choices)])
    elif ((d_b_choice == 'stay_enter') and kap_leave_notenter) \
            or ((d_b_choice == 'leave_notenter') and kap_leave_notenter):
        # (Detailed balance stay) AND (Kapitza leave)
        # OR
        # (Detailed balance leave) AND (Kapitza leave)
        # walk away, checking that current CNT volume is not a possibility
        possible_locs, num_possible_locs = generate_novol_choices_3d(grid, moves_3d, cur_pos, cur_index, kapitza,
                                                                     return_pos=True)
        final_pos = np.asarray(possible_locs[np.random.randint(0, num_possible_locs)])
    else:
        kill()
    #
    inside_cnt = not inside_cnt
    #
    return final_pos, inside_cnt


def tunneling_matrix(grid, moves_3d, cur_pos, cur_index):
    # checking that no CNT volume is a possibility
    possible_locs, num_possible_locs = generate_novol_choices_3d(grid, moves_3d, cur_pos, cur_index, False,
                                                                 return_pos=True)  # turning off Kapitza for tunneling
    final_pos = np.asarray(possible_locs[np.random.randint(0, num_possible_locs)])
    return final_pos


def tunneling_cntend(grid, jump_moves_3d, cur_pos):
    # no CNT volume, so this rule remains unchanged from the originals (3/30/17 TB)
    # walk off either end, 12 3D choices always
    choice = np.random.randint(0, 12)
    d_pos = np.asarray(jump_moves_3d[choice])
    # coord on left tube end jumps to right end
    tube_check_val_l = grid.tube_check_l[cur_pos[0], cur_pos[1], cur_pos[2]]
    # coord on right tube end jumps to left end
    tube_check_val_r = grid.tube_check_r[cur_pos[0], cur_pos[1], cur_pos[2]]
    if (tube_check_val_l > 0) and (tube_check_val_r == 0):
        # check that pixel cannot be a left and right endpoint
        if choice <= 5:  # stay at left end
            final_pos = cur_pos + d_pos
        else:  # jump across tube to right end
            # -1 BELOW BECAUSE OF +1 OFFSET IN CREATION TO AVOID ZERO INDEX
            final_pos = np.asarray(grid.tube_coords_r[tube_check_val_l - 1]) + np.asarray(d_pos)
    elif (tube_check_val_r > 0) and (tube_check_val_l == 0):
        # check that pixel cannot be a left and right endpoint
        if choice <= 5:  # stay at right end
            final_pos = cur_pos + d_pos
        else:  # jump across tube to left end
            # -1 BELOW BECAUSE OF +1 OFFSET IN CREATION TO AVOID ZERO INDEX
            final_pos = np.asarray(grid.tube_coords_l[tube_check_val_r - 1]) + np.asarray(d_pos)
    else:
        kill()
    return final_pos
