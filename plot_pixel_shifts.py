#!/usr/bin/env python
import os
import fnmatch
import numpy as np
from astropy.io import fits
import matplotlib.pyplot as plt
import argparse
from plot_images import plot_images

plot_images()


def get_phot_files(directory, phot_file_pattern):
    phot_files = []
    for root, _, files in os.walk(directory):
        for filename in fnmatch.filter(files, phot_file_pattern):
            phot_files.append(os.path.join(root, filename))
    return sorted(phot_files)


def read_phot_file(phot_file):
    try:
        with fits.open(phot_file) as ff:
            return ff[1].data
    except Exception as e:
        print(f"Error reading {phot_file}: {e}")
        return None


def get_pixel_positions_for_tic(phot_files, tic_id):
    """
    Collect x and y pixel positions for a specific TIC ID across all files.
    Returns two lists: x_positions, y_positions
    """
    x_positions = []
    y_positions = []

    for phot_file in phot_files:
        tab = read_phot_file(phot_file)
        if tab is None:
            continue

        # Find all rows where tic_id matches (elementwise)
        matches = tab['tic_id'] == tic_id
        if np.any(matches):
            x_positions.extend(tab['x'][matches])
            y_positions.extend(tab['y'][matches])

    return x_positions, y_positions


def plot_absolute_pixel_positions(x_positions, y_positions, tic_id):
    """
    Plot absolute pixel positions in two subplots (X and Y) for a given TIC ID.

    Parameters
    ----------
    x_positions : list or array-like
        List of X pixel positions over time or image index.
    y_positions : list or array-like
        List of Y pixel positions over time or image index.
    tic_id : str or int
        TIC ID of the star.
    """
    fig, axs = plt.subplots(2, 1, figsize=(8, 6), sharex=True, dpi=100)

    axs[0].plot(x_positions, 'b.')
    axs[0].set_ylabel("X Pixel")
    axs[0].set_title(f"X and Y Pixel Positions for TIC ID {tic_id}")
    axs[0].grid(True)
    axs[0].set_ylim(np.mean(x_positions) - 1.5, np.mean(x_positions) + 1.5)  # Keep the same Y limits for both plots

    axs[1].plot(y_positions, 'r.')
    axs[1].set_xlabel("Image Index")
    axs[1].set_ylabel("Y Pixel")
    axs[1].grid(True)
    axs[1].set_ylim(np.mean(y_positions) - 1.5, np.mean(y_positions) + 1.5)  # Keep the same Y limits for both plots

    plt.tight_layout()
    plt.show()


def main():
    parser = argparse.ArgumentParser(description="Plot absolute pixel positions for a given TIC ID.")
    parser.add_argument("phot_file", type=str, help="Photometry file pattern (e.g. 'phot_*.fits')")
    parser.add_argument("tic_id", type=int, help="TIC ID of the target star")
    args = parser.parse_args()

    directory = '/Users/u5500483/Downloads/'
    phot_files = get_phot_files(directory, args.phot_file)
    if not phot_files:
        print("No photometry files found.")
        return

    x_positions, y_positions = get_pixel_positions_for_tic(phot_files, args.tic_id)
    if not x_positions:
        print(f"TIC ID {args.tic_id} not found in any photometry file.")
        return

    plot_absolute_pixel_positions(x_positions, y_positions, args.tic_id)


if __name__ == "__main__":
    main()