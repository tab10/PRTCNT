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


from __future__ import division  # this ALWAYS gives float division for integers
import numpy as np
import logging
from scipy import stats


def final_conductivity_onlat(cur_dir, prob_m_cn, dt_dx_list, k_list, k_conv_error_buffer):
    """Final conductivity calculation, the best way to do this is averaging the last so many k values"""
    # heat_flux - [# walkers]/([time][length]**2)
    # dT(x)/dx - [# walkers]/[length]
    # k - 1/([time][length])
    k_conv_error_buffer = int(k_conv_error_buffer)
    k_mean = np.mean(k_list[-k_conv_error_buffer:])
    k_std = np.std(k_list[-k_conv_error_buffer:], ddof=1)
    dt_dx_mean = np.mean(dt_dx_list[-k_conv_error_buffer:])
    dt_dx_std = np.std(dt_dx_list[-k_conv_error_buffer:], ddof=1)
    logging.info("Average dT(x)/dx: %.4E +/- %.4E" % (dt_dx_mean, dt_dx_std))
    logging.info("Conductivity: %.4E +/- %.4E" % (k_mean, k_std))

    f = open("%s/k.txt" % cur_dir, 'w')
    f.write("%.4E\n" % k_mean)
    f.close()

    f = open("%s/prob_m_cn.txt" % cur_dir, 'w')
    f.write("%.4E\n" % prob_m_cn)
    f.close()


def check_convergence_2d_onlat(H_tot, cur_num_walkers, grid_size, timesteps):
    temp_profile = H_tot
    cutoff_dist = int(0.05 * grid_size)  # k is very sensitive to this, give it the most possible data to fit
    test_mean = np.mean(temp_profile[cutoff_dist:grid_size - cutoff_dist], axis=1)
    test_std = np.std(temp_profile[cutoff_dist:grid_size - cutoff_dist], axis=1, ddof=1)
    heat_flux = float(cur_num_walkers) / (float(grid_size + 2) * float(timesteps))
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
    cutoff_dist = int(0.05 * grid_size)  # k is very sensitive to this, give it the most possible data to fit
    test_mean = np.mean(temp_profile_sum[cutoff_dist:grid_size - cutoff_dist], axis=1)
    test_std = np.std(temp_profile_sum[cutoff_dist:grid_size - cutoff_dist], axis=1, ddof=1)
    heat_flux = float(cur_num_walkers) / (float(grid_size + 2) ** 2 * float(timesteps))
    gradient_err = np.mean(test_std)
    x = np.arange(cutoff_dist, grid_size - cutoff_dist)
    slope, intercept, r_value, p_value, std_err = stats.linregress(x, test_mean)  # slope is dT(x)/dx
    k = - heat_flux / slope
    k_err = (heat_flux / slope ** 2) * gradient_err
    return slope, heat_flux, gradient_err, k, k_err, r_value ** 2, temp_profile_sum


def rules_test_analysis(H_tot, cur_num_walkers, tot_time):
    """As this is the rules test, temperature is times visited by walkers. Works for 2D and 3D."""
    temp_profile = H_tot
    # generate histogram in terms of probability of visiting a spot (normalize)
    temp_profile_norm = np.divide(H_tot, (float(cur_num_walkers) * float(tot_time)))
    # averaging
    mean_temp = np.mean(temp_profile)
    std_temp = np.std(temp_profile, ddof=1)
    mean_temp_norm = np.mean(temp_profile_norm)
    std_temp_norm = np.std(temp_profile_norm, ddof=1)
    return mean_temp, mean_temp_norm, std_temp, std_temp_norm, temp_profile_norm
