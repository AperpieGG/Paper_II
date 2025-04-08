#!/usr/bin/env python

import batman
import numpy as np
import matplotlib.pyplot as plt
import json
from plot_images import plot_images
import argparse

# Argument parser for input
parser = argparse.ArgumentParser(description='Plot the transit model for a given TIC ID and camera numbers.')
parser.add_argument('cam1', type=str, help='First camera number (CCD or CMOS)')
parser.add_argument('cam2', type=str, help='Second camera number (CCD or CMOS)')
parser.add_argument('target', type=str, help='Target name')
args = parser.parse_args()
cam1 = args.cam1
cam2 = args.cam2
target = args.target

plot_images()


# Function to process data and generate plots
def process_camera(cam, target):
    # Set up the transit parameters
    params = batman.TransitParams()
    params.t0 = 2454823.58592  # time of inferior conjunction (BJD)
    params.per = 1.338231602  # orbital period (days)
    params.rp = 0.155241747  # planet radius (in units of stellar radii)
    params.a = 5.451  # semi-major axis (in units of stellar radii)
    params.inc = 88.39  # orbital inclination (degrees)
    params.ecc = 0.0053  # eccentricity
    params.w = 288  # longitude of periastron (degrees)
    params.u = [0.4451, 0.2297]  # limb darkening coefficients [u1, u2]
    params.limb_dark = "quadratic"  # limb darkening model

    # Load data from JSON file
    with open(f'data_transit/target_light_curve_{target}_{cam}_0802.json', 'r') as json_file:
        data = json.load(json_file)

    # Extract data
    time = np.array(data['Time_BJD'])
    flux = np.array(data['Relative_Flux'])
    flux_err = np.array(data['Relative_Flux_err'])

    # # Exclude last 300 points only for the first camera
    # if cam == cam1:
    #     time = time[:-500]
    #     flux = flux[:-500]
    #     flux_err = flux_err[:-500]

    # Bin the data
    time_binned, flux_binned, fluxerr_binned = bin_time_flux_error(time, flux, flux_err, 30)

    print(f'The maximum value is: {np.max(flux_binned)} ')
    print(f'The minimum value is: {np.min(flux_binned)}')
    print(f'The transit depth is : {(np.max(flux_binned) - np.min(flux_binned)) * 1000000}')

    # Calculate mean of normalised data
    mean_dt_flux = np.mean(flux_binned)
    adjustment = mean_dt_flux - 0.9948

    # Adjust flux
    dt_flux_adjusted = flux_binned - adjustment
    flux_adjusted = flux - adjustment

    # Normalize time array to be centered around the transit
    time_centered = time - params.t0

    # Initialize the transit model with the centered time array
    m = batman.TransitModel(params, time_binned)
    model_flux = m.light_curve(params)

    return time, flux_adjusted, flux_err, time_binned, dt_flux_adjusted, fluxerr_binned, model_flux


def bin_time_flux_error(time, flux, error, bin_fact):
    """
    Use reshape to bin light curve data, clip under filled bins
    Works with 2D arrays of flux and errors

    Note: under filled bins are clipped off the end of the series

    Parameters
    ----------
    time : array         of times to bin
    flux : array         of flux values to bin
    error : array         of error values to bin
    bin_fact : int
        Number of measurements to combine

    Returns
    -------
    times_b : array
        Binned times
    flux_b : array
        Binned fluxes
    error_b : array
        Binned errors

    Raises
    ------
    None
    """
    n_binned = int(len(time) / bin_fact)
    clip = n_binned * bin_fact
    time_b = np.average(time[:clip].reshape(n_binned, bin_fact), axis=1)
    # determine if 1 or 2d flux/err inputs
    if len(flux.shape) == 1:
        flux_b = np.average(flux[:clip].reshape(n_binned, bin_fact), axis=1)
        error_b = np.sqrt(np.sum(error[:clip].reshape(n_binned, bin_fact) ** 2, axis=1)) / bin_fact
    else:
        # assumed 2d with 1 row per star
        n_stars = len(flux)
        flux_b = np.average(flux[:clip].reshape((n_stars, n_binned, bin_fact)), axis=2)
        error_b = np.sqrt(np.sum(error[:clip].reshape((n_stars, n_binned, bin_fact)) ** 2, axis=2)) / bin_fact
    return time_b, flux_b, error_b


# Process both cameras
(time1, flux_adjusted1, flux_err1, time_binned1, dt_flux_adjusted1,
 fluxerr_binned1, model_flux1) = process_camera(cam1, target)
(time2, flux_adjusted2, flux_err2, time_binned2, dt_flux_adjusted2,
 fluxerr_binned2, model_flux2) = process_camera(cam2, target)

# Plotting
fig, axes = plt.subplots(1, 2, figsize=(12, 4), dpi=100, sharey=True, constrained_layout=False)
plt.subplots_adjust(wspace=0)  # Remove space between plots

# First subplot
axes[0].plot(time1, flux_adjusted1, '.', label=f"{cam1} Unbinned", color="grey", alpha=0.5)
axes[0].plot(time_binned1, dt_flux_adjusted1, 'o', label=f"{cam1} 5 min bin", color="red")
axes[0].plot(time_binned1, model_flux1, label=f"{cam1} Transit Model", color="black", linestyle='-')
axes[0].set_xlabel("Time (BJD)")
axes[0].set_ylabel("Relative flux")
axes[0].set_ylim(0.96, 1.03)
# axes[0].set_xlim(2460551 .591, 2460551.84)
axes[0].text(
    0.03, 0.05, "CMOS", transform=axes[0].transAxes,
    fontsize=12, verticalalignment='bottom', horizontalalignment='left',
    bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3')
)

# Second subplot
axes[1].plot(time2, flux_adjusted2, '.', label=f"{cam2} Unbinned", color="grey", alpha=0.5)
axes[1].plot(time_binned2, dt_flux_adjusted2, 'o', label=f"{cam2} 5 min bin", color="red")
axes[1].plot(time_binned2, model_flux2, label=f"{cam2} Transit Model", color="black", linestyle='-')
axes[1].set_xlabel("Time (BJD)")
axes[1].set_ylim(0.96, 1.03)
# axes[1].set_xlim(2460551.591, 2460551.84)
axes[1].text(
    0.03, 0.05, "CCD", transform=axes[1].transAxes,
    fontsize=12, verticalalignment='bottom', horizontalalignment='left',
    bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3')
)

# Save and show
plt.savefig(f'comparison_lc_{target}.pdf', dpi=300)
plt.show()