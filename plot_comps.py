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
    filename_1 = 'files/rms_vs_mag_comps.json'
    with open(filename_1, 'r') as file:
        data = json.load(file)

    return data


def plot_noise_model(data):
    """
    Plot RMS noise model on a given axis.
    """
    good_rms = np.array(data['good_rms']) * 1e6
    good_mags = data['good_mags']
    bad_rms = np.array(data['bad_rms']) * 1e6
    bad_mags = data['bad_mags']

    plt.figure()
    plt.scatter(good_mags, good_rms, color='black')
    plt.scatter(bad_mags, bad_rms, color='red')
    plt.xlabel('TESS Magnitude')
    plt.ylabel('RMS (ppm)')
    plt.yscale('log')
    plt.ylim(5000, 200000)  # Set y-axis limit based on dimmest good star RMS
    plt.xlim(11, 13.5)
    plt.gca().invert_xaxis()
    plt.ylim()
    # plt.title('RMS vs. Magnitude of Comparison Stars')
    plt.tight_layout()
    plt.savefig('rms_vs_mag_comps.pdf', dpi=300)
    plt.show()


if __name__ == "__main__":
    data = load_rms_mags_data()
    plot_noise_model(data)
