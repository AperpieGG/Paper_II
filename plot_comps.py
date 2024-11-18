#!/usr/bin/env python
import json
import matplotlib.pyplot as plt
from plot_images import plot_images

plot_images()


def load_rms_mags_data():
    """
    Load RMS and magnitude data from JSON file.
    """
    filename_1 = 'rms_vs_mag_comps.json'
    with open(filename_1, 'r') as file:
        data = json.load(file)

    return data


def plot_noise_model(data):
    """
    Plot RMS noise model on a given axis.
    """
    good_rms = data['good_rms']
    good_mags = data['good_mags']
    bad_rms = data['bad_rms']
    bad_mags = data['bad_mags']

    dimmest_good_star_rms = min(good_rms)
    y_limit_high = 3 * dimmest_good_star_rms  # Adjust multiplier as needed for clarity
    y_limit_low = 0.01 * dimmest_good_star_rms  # Adjust multiplier as needed for clarity

    plt.figure()
    plt.scatter(good_mags, good_rms, color='black')
    plt.scatter(bad_mags, bad_rms, color='red')
    plt.xlabel('TESS Magnitude')
    plt.ylabel('RMS')
    plt.ylim(0.01, 0.06)  # Set y-axis limit based on dimmest good star RMS
    plt.xlim(10.6, 13)
    plt.ylim()
    # plt.title('RMS vs. Magnitude of Comparison Stars')
    plt.tight_layout()
    plt.savefig('rms_vs_mag_comps.pdf', dpi=300)
    plt.show()


if __name__ == "__main__":
    data = load_rms_mags_data()
    plot_noise_model(data)
