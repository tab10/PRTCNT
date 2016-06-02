import numpy as np
from mpi4py import MPI
import setup

comm = MPI.COMM_WORLD
rank = comm.Get_rank()


def runrandomwalk2d(walkers_pos, grid_size, timesteps, bd_condition, walker_start_dir):
    # Start the random walk
    step = 0
    moves = [[0, 1], [1, 0], [0, -1], [-1, 0]]
    for i in range(1,timesteps+1):
        step += 1
        print "On step {0}".format(step)
        d_pos = []
        for j in range(len(walkers_pos)):
            d_pos.append(moves[np.random.randint(0, 3)])
            # array of all new positions for all walkers
        updatewalkers2d(walkers_pos, d_pos)
        applybdconditions(walkers_pos, bd_condition, grid_size, walker_start_dir)
    print 'Done in %d steps' % step
    return walkers_pos


def updatewalkers2d(pos_all, newpos):
    pos_all.append(pos_all[-1])  # copies last timestep to current one
    for i in range(len(pos_all[-1])):
        pos_all[-1][i].updatepos(newpos)
    return pos_all


def applybdconditions(pos_all, bd_condition, grid_size, walker_start_dir):
    for i in range(len(pos_all[-1])):
        cur_pos = pos_all[-1][i].pos
        if bd_condition == 'reflect':
            if cur_pos[0] > grid_size:  # x coordinate
                cur_pos[0] = grid_size
            elif cur_pos[0] < -grid_size:
                cur_pos[0] = -grid_size
            if cur_pos[1] > grid_size: # y coordinate
                cur_pos[1] = grid_size
            elif cur_pos[1] < -grid_size:
                cur_pos[1] = -grid_size
            pos_all[-1][i].updatepos(cur_pos)
        elif bd_condition == 'exit':
            pos_all[-1][i] = setup.Walker2D(grid_size, walker_start_dir) # starts back at random origin point
        else:
            print('Check bd condition argument.')
            raise SystemExit
    return pos_all
            