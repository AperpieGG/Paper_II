#!/usr/bin/env python
import json
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import sigmaclip  # Import for sigma clipping
from plot_images import plot_images

plot_images()


def load_rms_mags_data():
    """
    Load RMS and magnitude data from a single JSON file.
    """
    filename = 'files/rms_mags_rel_phot_NG2320-1302_4_1_0705.json'
    with open(filename, 'r') as file:
        data = json.load(file)
    return data


def apply_sigma_clip(deviations, colors, mags, sigma_threshold=3.0):
    """
    Apply sigma clipping to deviation data and filter corresponding colors and magnitudes.
    """
    deviations = np.array(deviations)
    colors = np.array(colors)
    mags = np.array(mags)

    # Remove NaNs before sigma clipping
    valid_indices = ~np.isnan(deviations)
    deviations = deviations[valid_indices]
    colors = colors[valid_indices]
    mags = mags[valid_indices]

    # Select only stars with colors in the range [0.7, 1.2] for sigma clipping
    clip_mask = (colors >= 0.3) & (colors <= 2)

    deviations_clip = deviations[clip_mask]
    colors_clip = colors[clip_mask]
    mags_clip = mags[clip_mask]

    if len(deviations_clip) > 0:  # Only apply clipping if there are enough data points
        clipped_devs, lower, upper = sigmaclip(deviations_clip, low=sigma_threshold, high=sigma_threshold)

        # Keep only data within sigma-clipped range
        clip_valid = (deviations_clip >= lower) & (deviations_clip <= upper)

        deviations_clip = deviations_clip[clip_valid]
        colors_clip = colors_clip[clip_valid]
        mags_clip = mags_clip[clip_valid]

    # Keep all original data outside the 0.7 - 1.2 color range
    non_clip_mask = ~clip_mask  # Select data outside the color range
    deviations_final = np.concatenate([deviations_clip, deviations[non_clip_mask]])
    colors_final = np.concatenate([colors_clip, colors[non_clip_mask]])
    mags_final = np.concatenate([mags_clip, mags[non_clip_mask]])

    return colors_final, deviations_final, mags_final


def compute_deviation(mags, rms, RNS, synthetic_mag):
    """
    Compute the deviation of RMS values from RNS values correctly for a log-scale y-axis.

    Deviation is computed as:
        deviation = | log10(RMS_measured) - log10(RNS_expected) |
    """
    deviations = []
    for i in range(len(mags)):
        closest_rns_index = (np.abs(synthetic_mag - mags[i])).argmin()
        expected_rns = RNS[closest_rns_index]

        # Ensure values are positive before taking log10
        if rms[i] > 0 and expected_rns > 0:
            deviation = np.abs(np.log10(rms[i]) - np.log10(expected_rns))  # Log-scale absolute difference
        else:
            deviation = np.nan  # Assign NaN if invalid values occur
        deviations.append(deviation)

    return deviations


def plot_noise_model(ax, data, label):
    """
    Plot RMS noise model on a given axis, with main and moon phase labels.
    """
    RMS_list = np.array(data['RMS_list'])
    Tmag_list = np.array(data['Tmag_list'])
    color_list = np.array(data['COLOR'])
    synthetic_mag = np.array(data['synthetic_mag'])
    RNS = np.array(data['RNS'])

    # Print average scintillation noise
    print(f'The average scintillation noise is: {np.mean(data["N"])}')

    # Filter out stars with missing color information
    valid_indices = [i for i in range(len(Tmag_list)) if color_list[i] is not None and Tmag_list[i] > 9]

    # Extract valid data
    mags_all = Tmag_list[valid_indices]
    rms_all = RMS_list[valid_indices]
    colors_all = color_list[valid_indices]

    # Compute deviations for all stars
    deviations_all = compute_deviation(mags_all, rms_all, RNS, synthetic_mag)

    # Scatter plot of stars
    scatter = ax.scatter(mags_all, rms_all, c=colors_all, cmap='coolwarm', vmin=0.5, vmax=1.5)

    # Add a colorbar
    cbar = plt.colorbar(scatter, ax=ax, label='Color index (G_BP - G_RP)')

    # Plot various noise sources
    ax.plot(synthetic_mag, RNS, color='black', label='Total noise')
    ax.set_yscale('log')
    ax.set_xlim(7.5, 14)
    ax.set_ylim(1000, 100000)
    ax.invert_xaxis()
    ax.legend()

    # Add labels
    ax.text(
        0.02, 0.03, label, transform=ax.transAxes,
        fontsize=15, fontweight='bold', ha='left', va='bottom',
        bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3')
    )

    return colors_all, deviations_all, mags_all


def main():
    # Load a single dataset
    data = load_rms_mags_data()

    # Create figure for noise model
    fig, ax = plt.subplots(figsize=(8, 6))
    label = "CMOS"  # Adjust the label as needed

    # Store data for deviation plot
    colors_all, deviations_all, mags_all = plot_noise_model(ax, data, label)
    colors_filtered, deviations_filtered, mags_filtered = apply_sigma_clip(deviations_all, colors_all, mags_all)

    # Create a new figure for deviation plot
    fig_dev, ax_dev = plt.subplots(figsize=(8, 6))
    sc = ax_dev.scatter(colors_filtered, deviations_filtered, c=mags_filtered, cmap='viridis')

    # Add a colorbar
    cbar = plt.colorbar(sc, ax=ax_dev, label='TESS Magnitude')

    ax_dev.set_xlabel('$\mathdefault{G_{BP}-G_{RP}}$')
    ax_dev.set_ylabel('Deviation from Model (ppm)')
    # set y limits
    ax_dev.set_ylim(-0.2, 2)

    plt.show()


if __name__ == "__main__":
    main()