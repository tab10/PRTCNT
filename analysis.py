import numpy as np
import logging


def calc_conductivity_2d_onlat(num_walkers, temp_gradient_x, grid_size, timesteps):
    """Plots temperature gradient for all walkers"""
    gradient_avg = float(np.mean(temp_gradient_x))
    heat_flux = (2 * float(num_walkers)) / (float(grid_size) ** 2 * float(timesteps))
    k = - heat_flux / gradient_avg
    logging.info("Average dT(x)/dx: %.4E" % gradient_avg)
    logging.info("Heat flux: %.4E" % heat_flux)
    logging.info("Conductivity: %.4E" % k)
