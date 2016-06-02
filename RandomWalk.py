

def two_d_random_walk(grid_size, tube_coords):
    ### Start the random walk
    step = 0
    moves = [(0, 1), (1, 0), (0, -1), (-1, 0)]
    x = []
    y = []
    x.append(x_i)
    y.append(y_i)
    while (x[-1] != x_f) or (y[-1] != y_f):
        step += 1
        dx, dy = moves[np.random.randint(0, 3)]
        x.append(x[-1] + dx)
        y.append(y[-1] + dy)
        # bd conditions
        if x[-1] > grid_size:
            x[-1] = -grid_size
        elif x[-1] < -grid_size:
            x[-1] = np.abs(grid_size)
        if y[-1] > grid_size:
            y[-1] = -grid_size
        elif y[-1] < -grid_size:
            y[-1] = np.abs(grid_size)
    print 'Done in %d steps' % step
    return x, y