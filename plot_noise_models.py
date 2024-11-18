#!/usr/bin/env python
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from plot_images import plot_images


def load_rms_mags_data():
    """
    Load RMS and magnitude data from JSON file.
    """
    filename_1 = 'rms_mags_rel_phot_NG2320-1302_4_1_0705.json'
    filename_2 = 'rms_mags_rel_phot_NG2320-1302_5_1_0705.json'
    with open(filename_1, 'r') as file:
        data_1 = json.load(file)

    with open(filename_2, 'r') as file:
        data_2 = json.load(file)

    return data_1, data_2


def plot_noise_model(ax, data):
    """
    Plot RMS noise model on a given axis.
    """
    RMS_list = data['RMS_list']
    Tmag_list = data['Tmag_list']
    color_list = data['COLOR']
    synthetic_mag = data['synthetic_mag']
    RNS = data['RNS']
    photon_shot_noise = data['photon_shot_noise']
    read_noise = data['read_noise']
    dc_noise = data['dc_noise']
    sky_noise = data['sky_noise']
    N = data['N']
    print(f'The average scintillation noise is: {np.mean(N)}')

    # Filter out stars with missing color information
    total_mags, total_RMS, total_colors = [], [], []
    for i in range(len(Tmag_list)):
        if color_list[i] is not None:  # Include only stars with color information
            total_mags.append(Tmag_list[i])
            total_RMS.append(RMS_list[i])
            total_colors.append(color_list[i])

    # Scatter plot with remaining stars
    scatter = ax.scatter(total_mags, total_RMS, c=total_colors, cmap='coolwarm', vmin=0.5, vmax=1.5)

    # Plot various noise sources
    ax.plot(synthetic_mag, RNS, color='black', label='total noise')
    ax.plot(synthetic_mag, photon_shot_noise, color='green', label='photon shot', linestyle='--')
    ax.plot(synthetic_mag, read_noise, color='red', label='read noise', linestyle='--')
    ax.plot(synthetic_mag, dc_noise, color='purple', label='dark noise', linestyle='--')
    ax.plot(synthetic_mag, sky_noise, color='blue', label='sky bkg', linestyle='--')
    ax.plot(synthetic_mag, np.ones(len(synthetic_mag)) * N, color='orange', label='scintillation noise', linestyle='--')

    ax.set_xlabel('TESS Magnitude')
    ax.set_yscale('log')
    ax.set_xlim(7.5, 14)
    ax.set_ylim(1000, 100000)
    ax.invert_xaxis()
    ax.xaxis.set_major_locator(MultipleLocator(2))  # Set x-axis ticks to step by 1
    return scatter


def main():
    # Set plot parameters
    plot_images()
    # Load RMS and magnitude data from JSON files
    data_1, data_2 = load_rms_mags_data()

    # Create a figure with two subplots (1 row, 2 columns) with shared y-axis
    fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, constrained_layout=True)

    # Plot RMS noise models for both datasets
    scatter1 = plot_noise_model(ax1, data_1)
    scatter2 = plot_noise_model(ax2, data_2)

    # Set the shared y-axis label only on the first plot
    ax1.set_ylabel('RMS (ppm)')

    # Create a single colorbar for both plots
    cbar = fig.colorbar(scatter1, ax=[ax1, ax2], orientation='vertical')
    cbar.set_label(label='$\mathdefault{G_{BP}-G_{RP}}$')
    plt.savefig('Noise_Dark.pdf', dpi=300)
    plt.show()


if __name__ == "__main__":
    main()