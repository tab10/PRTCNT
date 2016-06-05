import numpy as np
import logging
import numericalunits


def calc_conductivity_2d_onlat(num_walkers, temp_gradient_x, grid_size, timesteps):
    """Plots temperature gradient for all walkers"""
    # heat_flux - [# walkers]/([time][length]**2)
    # dT(x)/dx - [# walkers]/[length]
    # k - 1/([time][length])
    gradient_avg = float(np.mean(temp_gradient_x[3:]))
    # disregard first few x= slices as close to the wall and values have large errors
    heat_flux = (2 * float(num_walkers)) / (float(grid_size) ** 2 * float(timesteps))
    k = - heat_flux / gradient_avg
    logging.info("Average dT(x)/dx: %.4E" % gradient_avg)
    logging.info("Heat flux: %.4E" % heat_flux)
    logging.info("Conductivity: %.4E" % k)
    return k
