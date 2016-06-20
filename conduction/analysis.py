import numpy as np
import logging
import numericalunits
from scipy import stats


# def sim_setup_3d(grid_size, tube_diameter, l_d, )
#        tube_length = l_d * tube_diameter


def filling_fraction(tube_coords, grid_size):
    """Calculates pixels covered by the tubes"""
    x = []
    y = []
    fill_fract_list = []  # index corresponds to tube_coords index
    for i in range(len(tube_coords)):
        x.append(np.abs(tube_coords[i][2] - tube_coords[i][0]))
        y.append(np.abs(tube_coords[i][3] - tube_coords[i][1]))
        fill_fract_list.append(x[i] * y[i])
    fill_fract_sum = sum(fill_fract_list)
    fill_fract = fill_fract_sum / grid_size
    return fill_fract


def final_conductivity_2d_onlat(num_walkers, grid_size, timesteps, slope, gradient_err, k_err, num_tubes, cur_dir,
                                k_convergence_val, gradient_cutoff):
    """Final conductivity calculation"""
    # heat_flux - [# walkers]/([time][length]**2)
    # dT(x)/dx - [# walkers]/[length]
    # k - 1/([time][length])
    logging.info("Using gradient cutoff at x=%d" % gradient_cutoff)
    # gradient_avg = float(np.mean(temp_gradient_x[gradient_cutoff:]))
    # gradient_std = float(np.std(temp_gradient_x[gradient_cutoff:], ddof = 1))
    gradient_avg = slope / float(timesteps)
    gradient_std = gradient_err
    # disregard first few x= slices as close to the wall and values have large errors
    heat_flux = float(num_walkers) / (float(grid_size) ** 2 * float(timesteps))
    # k = - heat_flux / gradient_avg
    # k_err = (heat_flux/gradient_avg**2)*gradient_std
    k = k_convergence_val
    logging.info("Average dT(x)/dx: %.4E +/- %.4E" % (gradient_avg, gradient_std))
    logging.info("Heat flux: %.4E" % heat_flux)
    logging.info("Conductivity: %.4E +/- %.4E" % (k, k_err))
    f = open("%s/k.txt" % cur_dir, 'a')
    f.write("%d %.4E %.4E\n" % (num_tubes, k, k_err))
    f.close()
    return k


def final_conductivity_3d_onlat(num_walkers, grid_size, timesteps, slope, gradient_err, k_err, num_tubes, cur_dir,
                                k_convergence_val, gradient_cutoff):
    """Final conductivity calculation"""
    # heat_flux - [# walkers]/([time][length]**2)
    # dT(x)/dx - [# walkers]/[length]
    # k - 1/([time][length])
    logging.info("Using gradient cutoff at x=%d" % gradient_cutoff)
    # gradient_avg = float(np.mean(temp_gradient_x[gradient_cutoff:]))
    # gradient_std = float(np.std(temp_gradient_x[gradient_cutoff:], ddof = 1))
    gradient_avg = slope / float(timesteps)
    gradient_std = gradient_err
    # disregard first few x= slices as close to the wall and values have large errors
    heat_flux = float(num_walkers) / (float(grid_size) ** 3 * float(timesteps))
    # k = - heat_flux / gradient_avg
    # k_err = (heat_flux/gradient_avg**2)*gradient_std
    k = k_convergence_val
    logging.info("Average dT(x)/dx: %.4E +/- %.4E" % (gradient_avg, gradient_std))
    logging.info("Heat flux: %.4E" % heat_flux)
    logging.info("Conductivity: %.4E +/- %.4E" % (k, k_err))
    f = open("%s/k.txt" % cur_dir, 'a')
    f.write("%d %.4E %.4E\n" % (num_tubes, k, k_err))
    f.close()
    return k

def check_convergence_2d_onlat(H_tot, cur_num_walkers, grid_size, timesteps):
    temp_profile = H_tot
    test_mean = np.mean(temp_profile[1:], axis=1)
    test_std = np.std(temp_profile[1:], axis=1, ddof=1)
    heat_flux = float(cur_num_walkers) / (float(grid_size) ** 2 * float(timesteps))
    gradient_err = np.mean(test_std)
    x = np.arange(1, len(test_mean) + 1)
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, test_mean)  # slope is dT(x)/dx
    slope /= float(timesteps)
    k = - heat_flux / slope
    k_err = (heat_flux / slope ** 2) * gradient_err
    return slope, heat_flux, gradient_err, k, k_err, r_value**2


def check_convergence_3d_onlat(H_tot, cur_num_walkers, grid_size, timesteps):
    temp_profile = H_tot  # temp_profile is 3D. Collapse y and z dimensions (periodic)
    temp_profile_sum = np.zeros((len(temp_profile), len(temp_profile)))
    for i in range(len(temp_profile)):
        for j in range(len(temp_profile)):
            temp_profile_sum[i][j] = np.sum(temp_profile[i][j])
    test_mean = np.mean(temp_profile_sum[1:], axis=1)
    test_std = np.std(temp_profile_sum[1:], axis=1, ddof=1)
    heat_flux = float(cur_num_walkers) / (float(grid_size) ** 3 * float(timesteps))
    gradient_err = np.mean(test_std)
    x = np.arange(1, len(test_mean) + 1)
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, test_mean)  # slope is dT(x)/dx
    slope /= float(timesteps)
    k = - heat_flux / slope
    k_err = (heat_flux / slope ** 2) * gradient_err
    return slope, heat_flux, gradient_err, k, k_err, r_value ** 2, temp_profile_sum