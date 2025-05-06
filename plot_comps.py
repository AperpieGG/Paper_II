#!/usr/bin/env python
import json
import matplotlib.pyplot as plt
import numpy as np

from plot_images import plot_images

plot_images()


def load_rms_mags_data():
    """
    Load RMS and magnitude data from JSON file.
    """
    filename_1 = 'files/rms_vs_mag_comps_0802.json'
    with open(filename_1, 'r') as file:
        data = json.load(file)

    return data


def remove_lowest_point(data, lower_bound, upper_bound):
    """
    Find and remove the lowest RMS point between specified magnitude bounds.
    """
    # Convert RMS to numpy arrays
    good_rms = np.array(data['good_rms'])
    good_mags = np.array(data['good_mags'])

    # Find points within the specified Tmag range
    mask = (good_mags >= lower_bound) & (good_mags <= upper_bound)
    rms_in_range = good_rms[mask]
    mags_in_range = good_mags[mask]

    if len(rms_in_range) > 0:
        # Identify the index of the lowest RMS within the range
        lowest_index = np.argmin(rms_in_range)
        lowest_rms = rms_in_range[lowest_index]
        lowest_mag = mags_in_range[lowest_index]
        print(f"Removing lowest point: Tmag={lowest_mag}, RMS={lowest_rms}")

        # Remove the lowest point from good_rms and good_mags
        keep_mask = ~((good_mags == lowest_mag) & (good_rms == lowest_rms))
        data['good_rms'] = good_rms[keep_mask].tolist()
        data['good_mags'] = good_mags[keep_mask].tolist()


def plot_noise_model(data):
    """
    Plot RMS noise model on a given axis.
    """
    good_rms = np.array(data['good_rms']) * 1e6
    good_mags = data['good_mags']
    bad_rms = np.array(data['bad_rms']) * 1e6
    bad_mags = data['bad_mags']

    plt.figure()
    plt.scatter(good_mags, good_rms, color='black', label='Good Stars')
    plt.scatter(bad_mags, bad_rms, marker='s', color='red', label='Bad Stars')
    plt.xlabel('TESS Magnitude')
    plt.ylabel('RMS (ppm)')
    plt.yscale('log')
    plt.ylim(5000, 200000)  # Set y-axis limit based on dimmest good star RMS
    plt.xlim(11, 13.5)
    plt.gca().invert_xaxis()
    # plt.legend()
    plt.tight_layout()
    plt.savefig('rms_vs_mag_comps.pdf', dpi=300)
    plt.show()


if __name__ == "__main__":
    # Load the data
    data = load_rms_mags_data()

    # Remove the lowest RMS point in the specified range
    remove_lowest_point(data, lower_bound=12.6, upper_bound=12.8)

    # Plot the updated data
    plot_noise_model(data)