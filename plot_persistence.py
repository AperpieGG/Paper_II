"""
Description: This script is used to take image data and plot the results.
The image data is recorded with the Marana.

It measures the mean dark signal from the FFI for each image exposure at different temperatures.
It then finds the gradient/slope for each different temperature and estimates the dark current.
It plots the dark current as a function of temperature.
The errors for each slope for the dark current vs temperature are estimated using the covariance matrix.
"""

from astropy.io import fits
import numpy as np
import glob
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def plot_images():
    plt.rcParams['figure.dpi'] = 300
    plt.rcParams['xtick.top'] = True
    plt.rcParams['xtick.labeltop'] = False
    plt.rcParams['xtick.labelbottom'] = True
    plt.rcParams['xtick.bottom'] = True
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['xtick.minor.visible'] = True
    plt.rcParams['xtick.major.top'] = True
    plt.rcParams['xtick.minor.top'] = True
    plt.rcParams['xtick.minor.bottom'] = True
    plt.rcParams['xtick.alignment'] = 'center'

    plt.rcParams['ytick.left'] = True
    plt.rcParams['ytick.labelleft'] = True
    plt.rcParams['ytick.right'] = True
    plt.rcParams['ytick.minor.visible'] = True
    plt.rcParams['ytick.major.right'] = True
    plt.rcParams['ytick.major.left'] = True
    plt.rcParams['ytick.minor.right'] = True
    plt.rcParams['ytick.minor.left'] = True

    plt.rcParams['font.family'] = 'Times New Roman'
    plt.rcParams['font.size'] = 14

    plt.rcParams['legend.frameon'] = True
    plt.rcParams['legend.framealpha'] = 0.8
    plt.rcParams['legend.loc'] = 'best'
    plt.rcParams['legend.fancybox'] = True
    plt.rcParams['legend.fontsize'] = 14


def get_dark_data(path, suffix):
    dark_list = glob.glob(path + f'/DATA_{suffix}/dark*.fits')
    print("The shape of the dark_list is: ", np.shape(dark_list))
    dark_data = [fits.getdata(dark)[0] for dark in dark_list]
    dark_data = dark_data + np.array(99.26)
    return np.array(dark_data)


def get_exposure_values(path, suffix):
    dark_list = glob.glob(path + f'/DATA_{suffix}/dark*.fits')
    exposure = [fits.getheader(dark)['EXPOSURE'] for dark in dark_list]
    exposure = np.array(exposure)
    return exposure


def get_mean_dark_MID(dark_data):
    return [np.mean(dark) for dark in dark_data]


def get_std_dark_MID(dark_data):
    return [np.std(dark) for dark in dark_data]


def linear_fit(x, a, b):
    return a * x + b


def fit_in_data(exposure_times, mean_dark_current):
    popt, pcov = curve_fit(linear_fit, exposure_times, mean_dark_current)
    return popt, pcov


def plot_ds_exp_and_residuals(exposure, mean_dark_MID_list, std_dark_MID_list, temperature, slope_errors):
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, gridspec_kw={'height_ratios': [3, 1]})
    cmap = plt.get_cmap('cividis')
    colors = cmap(np.linspace(0, 1, len(temperature)))
    selected_temps = [-25, -45, 15]
    for temp, mean_dark_MID, std_dark_MID, color, slope_err in zip(temperature, mean_dark_MID_list, std_dark_MID_list,
                                                                   colors, slope_errors):
        if temp in selected_temps:
            popt, pcov = fit_in_data(exposure, mean_dark_MID)
            fitted_line = linear_fit(exposure, *popt)
            residuals = mean_dark_MID - fitted_line
            ax1.plot(exposure, fitted_line, color=color, linestyle='-')
            ax2.errorbar(exposure, residuals, yerr=std_dark_MID, fmt='o', color=color,
                         label=f'{temp} \N{DEGREE SIGN}C Residuals')
            ax1.errorbar(exposure, mean_dark_MID, yerr=std_dark_MID, fmt='o', color=color,
                         label=f'{temp} \N{DEGREE SIGN}C')
            ax1.errorbar(exposure, fitted_line, yerr=slope_err, fmt='o', color=color)

    ax1.set_ylabel('Dark signal (ADU)')
    ax2.set_xlabel('Exposure Time (s)')
    ax2.set_ylabel('Residuals (ADU)')

    # Set x-ticks for both axes
    ax1.set_xticks(np.arange(1, 11, 1))
    ax2.set_xticks(np.arange(1, 11, 1))

    # Create a custom colorbar for ax2
    norm = plt.Normalize(vmin=-55, vmax=15)
    sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
    sm.set_array([])  # dummy variable required for ScalarMappable
    cbar = fig.colorbar(sm, ax=ax2, ticks=selected_temps, aspect=6.5, pad=0.01, fraction=0.055, shrink=1)
    temp_ticks = [-50, -40, -30, -20, -10, 0, 10]
    cbar_1 = fig.colorbar(plt.cm.ScalarMappable(cmap=cmap, norm=plt.Normalize(vmin=-55, vmax=15)
                                                ), ax=ax1, ticks=temp_ticks, pad=0.01, fraction=0.057, shrink=1,
                          aspect=20,
                          extend='both', extendrect=True, extendfrac='auto')
    cbar_1.set_label('Temperature (\N{DEGREE SIGN}C)')
    plt.tight_layout()
    plt.savefig('DC_vs_Exposure.pdf', bbox_inches='tight', dpi=300)
    plt.show()


def main():
    plot_images()
    path = '/Users/u5500483/Documents/GitHub/Paper_I/Results/Images/DARK_DATA/'
    temperature = [-55, -50, -45, -40, -35, -30, -25, -20, -15, -10, -5, 0, 5, 10, 15]
    number_of_pixels = 4194304
    dark_data_list = [get_dark_data(path, suffix) for suffix in temperature]
    exposure_list = [get_exposure_values(path, suffix) for suffix in temperature]
    mean_dark_MID_list = [get_mean_dark_MID(dark_data) for dark_data in dark_data_list]
    std_dark_MID_list = [get_std_dark_MID(dark_data) for dark_data in dark_data_list] / np.sqrt(number_of_pixels)

    # Curve fitting
    popts, pcovs = zip(*[fit_in_data(exposure, mean_dark_MID) for exposure, mean_dark_MID
                         in zip(exposure_list, mean_dark_MID_list)])

    # Extract slope errors from the covariance matrix
    slope_errors = [np.sqrt(np.diag(pcov))[0] for pcov in
                    pcovs]  # Only take the slope errors (first element of diagonal)
    print('The slopes are:', popts)
    print('The errors are:', slope_errors)
    plot_ds_exp_and_residuals(exposure_list[0], mean_dark_MID_list, std_dark_MID_list, temperature, slope_errors)


if __name__ == '__main__':
    main()
