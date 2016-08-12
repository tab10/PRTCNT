import numpy as np
import logging
# import numericalunits
from scipy import stats


def final_conductivity_2d_onlat(num_walkers, grid_size, timesteps, slope, gradient_err, k_err, num_tubes, cur_dir,
                                k_convergence_val, prob_m_cn, gradient_cutoff):
    """Final conductivity calculation"""
    # heat_flux - [# walkers]/([time][length]**2)
    # dT(x)/dx - [# walkers]/[length]
    # k - 1/([time][length])
    logging.info("Using gradient cutoff at x=%d" % gradient_cutoff)
    # gradient_avg = float(np.mean(temp_gradient_x[gradient_cutoff:]))
    # gradient_std = float(np.std(temp_gradient_x[gradient_cutoff:], ddof = 1))
    gradient_avg = slope  # / float(timesteps)
    gradient_std = gradient_err
    # disregard first few x= slices as close to the wall and values have large errors
    heat_flux = float(num_walkers) / (float(grid_size + 1) ** 2)  # * float(timesteps))
    # k = - heat_flux / gradient_avg
    # k_err = (heat_flux/gradient_avg**2)*gradient_std
    k = k_convergence_val
    logging.info("Average dT(x)/dx: %.4E +/- %.4E" % (gradient_avg, gradient_std))
    logging.info("Heat flux: %.4E" % heat_flux)
    logging.info("Conductivity: %.4E +/- %.4E" % (k, k_err))

    f = open("%s/k.txt" % cur_dir, 'w')
    f.write("%.4E\n" % k)
    f.close()

    f = open("%s/prob_m_cn.txt" % cur_dir, 'w')
    f.write("%.4E\n" % prob_m_cn)
    f.close()
    return k


def final_conductivity_3d_onlat(num_walkers, grid_size, timesteps, slope, gradient_err, k_err, num_tubes, cur_dir,
                                k_convergence_val, prob_m_cn, gradient_cutoff):
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
    heat_flux = float(num_walkers) / (float(grid_size + 1) ** 2 * float(timesteps))
    # k = - heat_flux / gradient_avg
    # k_err = (heat_flux/gradient_avg**2)*gradient_std
    k = k_convergence_val
    logging.info("Average dT(x)/dx: %.4E +/- %.4E" % (gradient_avg, gradient_std))
    logging.info("Heat flux: %.4E" % heat_flux)
    logging.info("Conductivity: %.4E +/- %.4E" % (k, k_err))

    f = open("%s/k.txt" % cur_dir, 'w')
    f.write("%.4E\n" % k)
    f.close()

    f = open("%s/prob_m_cn.txt" % cur_dir, 'w')
    f.write("%.4E\n" % prob_m_cn)
    f.close()
    return k


def check_convergence_2d_onlat(H_tot, cur_num_walkers, grid_size, timesteps):
    temp_profile = H_tot
    cutoff_dist = int(0.25 * grid_size)
    test_mean = np.mean(temp_profile[cutoff_dist:grid_size - cutoff_dist], axis=1)
    test_std = np.std(temp_profile[cutoff_dist:grid_size - cutoff_dist], axis=1, ddof=1)
    # test_mean /= (float(timesteps)*float(cur_num_walkers))
    test_mean /= float(cur_num_walkers)
    # test_mean /= float(timesteps)
    heat_flux = float(cur_num_walkers) / ((float(grid_size + 1)) * float(timesteps))
    gradient_err = np.mean(test_std)
    x = np.arange(cutoff_dist, grid_size - cutoff_dist)
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, test_mean)  # slope is dT(x)/dx
    k = - heat_flux / slope
    k_err = (heat_flux / slope ** 2) * gradient_err
    return slope, heat_flux, gradient_err, k, k_err, r_value**2


def check_convergence_3d_onlat(H_tot, cur_num_walkers, grid_size, timesteps):
    temp_profile = H_tot  # temp_profile is 3D. Collapse y and z dimensions (periodic)
    temp_profile_sum = np.zeros((len(temp_profile), len(temp_profile)))
    for i in range(len(temp_profile)):
        for j in range(len(temp_profile)):
            temp_profile_sum[i][j] = np.sum(temp_profile[i][j])
    cutoff_dist = int(0.25 * grid_size)
    test_mean = np.mean(temp_profile_sum[cutoff_dist:grid_size - cutoff_dist], axis=1)
    test_std = np.std(temp_profile_sum[cutoff_dist:grid_size - cutoff_dist], axis=1, ddof=1)
    heat_flux = float(cur_num_walkers) / (float(grid_size + 1) ** 2 * float(timesteps))
    gradient_err = np.mean(test_std)
    x = np.arange(cutoff_dist, grid_size - cutoff_dist)
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, test_mean)  # slope is dT(x)/dx
    slope /= float(timesteps)
    k = - heat_flux / slope
    k_err = (heat_flux / slope ** 2) * gradient_err
    return slope, heat_flux, gradient_err, k, k_err, r_value ** 2, temp_profile_sum
