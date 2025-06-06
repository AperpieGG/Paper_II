#!/usr/bin/env python
import json
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

from plot_images import plot_images


def load_rms_mags_data():
    """
    Load RMS and magnitude data from JSON files.
    """
    filenames = [
        'files/rms_mags_rel_phot_NG2320-1302_4_1_0705.json',
        'files/rms_mags_rel_phot_NG2320-1302_5_1_0705.json',
        'files/rms_mags_rel_phot_NG2320-1302_4_1_0622.json',
        'files/rms_mags_rel_phot_NG2320-1302_5_1_0622.json'
    ]
    data_list = []
    for filename in filenames:
        with open(filename, 'r') as file:
            data_list.append(json.load(file))
    return data_list


def plot_noise_model(ax, data, label, no_moon_label=None, full_moon_label=None):
    """
    Plot RMS noise model on a given axis, with main and moon phase labels.
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
        if color_list[i] is not None:  # is not None and color_list[i] > 1.4
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

    ax.set_yscale('log')
    ax.set_xlim(7.5, 14)
    ax.set_ylim(1000, 100000)
    ax.invert_xaxis()

    # Add the main label in the bottom-left corner
    ax.text(
        0.02, 0.03, label, transform=ax.transAxes,
        fontsize=15, fontweight='bold', ha='left', va='bottom',
        bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3')
    )

    # Add the moon phase label, depending on the input
    if no_moon_label:
        ax.text(
            0.815, 0.08, no_moon_label, transform=ax.transAxes,
            fontsize=15, fontweight='bold', ha='left', va='top',
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3')
        )
    if full_moon_label:
        ax.text(
            0.8, 0.08, full_moon_label, transform=ax.transAxes,
            fontsize=15, fontweight='bold', ha='left', va='top',
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3')
        )

    return scatter


def main():
    # Set plot parameters
    plot_images()
    # Load RMS and magnitude data from JSON files
    data_list = load_rms_mags_data()

    # Create a figure with 2 rows and 2 columns of subplots
    fig, axes = plt.subplots(2, 2, sharey=True, constrained_layout=True, figsize=(12, 8))

    # Labels for the subplots
    labels = ["CCD", "CMOS", "CCD", "CMOS"]

    # Plot RMS noise models for all datasets
    scatters = []
    for i, ax in enumerate(axes.flat):
        # Determine moon labels based on the subplot's row
        no_moon_label = 'No-moon' if i < 2 else None
        full_moon_label = 'Full-moon' if i >= 2 else None

        # Plot the noise model with appropriate labels
        scatter = plot_noise_model(ax, data_list[i], labels[i], no_moon_label, full_moon_label)

        if i in [0, 2]:  # Set the shared y-axis label only on the first column
            ax.set_ylabel('RMS per 10 seconds (ppm)')
        if i in [2, 3]:  # Set the shared x-axis label only on the last row
            ax.set_xlabel('TESS Magnitude')
        scatters.append(scatter)

    # Create a single colorbar for the entire figure
    cbar = fig.colorbar(scatters[0], ax=axes, orientation='vertical', fraction=0.1, pad=0.05)
    cbar.set_label(label='$\mathdefault{G_{BP}-G_{RP}}$')

    # Save and show the plot
    plt.savefig('Noise_NEW_FULL.pdf', dpi=300)
    plt.show()


if __name__ == "__main__":
    main()