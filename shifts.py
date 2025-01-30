#!/usr/bin/env python3

"""
This script runs Donuts on a set of images, computes the shifts,
extracts 'DATE-OBS' from the headers, converts it to Julian Date,
and plots the x-y shifts as a function of time.
"""

import os
import sys
import glob
import warnings
from datetime import datetime
from matplotlib import pyplot as plt
from astropy.io import fits
from astropy.time import Time
from donuts import Donuts
from plot_images import plot_images

plot_images()

warnings.filterwarnings("ignore", category=UserWarning,
                        module="numpy.core.fromnumeric")
warnings.filterwarnings("ignore", category=UserWarning,
                        module="donuts.image")


def utc_to_jd(utc_time_str):
    """
    Convert UTC time string to Julian Date.

    Parameters
    ----------
    utc_time_str : str
        UTC time string in the format 'YYYY-MM-DDTHH:MM:SS.SSS'

    Returns
    -------
    jd : float
        Julian Date
    """
    t = Time(utc_time_str, format='isot', scale='utc')
    jd = t.jd
    return jd


def find_first_image_of_each_prefix(directory, save_path):
    # List all items in the given directory
    items = os.listdir(directory)
    # Filter out unwanted files
    filtered_items = [item for item in items
                      if "flat" not in item.lower()
                      and "bias" not in item.lower()
                      and "dark" not in item.lower()]
    # Dictionary to store the first image of each prefix
    first_image_of_each_prefix = {}
    # Words to exclude in the prefix
    exclude_words = ["evening", "morning", "flat", "bias", "dark"]
    # Iterate through filtered items
    for item in filtered_items:
        # Check if any exclude word is in the item
        if any(word in item.lower() for word in exclude_words):
            continue
        # Extract the first 11 characters as prefix
        prefix = item[:11]
        # Check if the prefix is already a key
        if prefix not in first_image_of_each_prefix:
            # Get matching files
            pattern = os.path.join(directory, f'{prefix}*.fits')
            matching_files = glob.glob(pattern)
            matching_files = sorted(matching_files)
            if matching_files:
                first_image_of_each_prefix[prefix] = matching_files[0]
    # Process each prefix
    for prefix, first_image in first_image_of_each_prefix.items():
        print(f"Processing prefix: {prefix}, First Image: {first_image}")
        run_donuts(directory, prefix, save_path)
    if not first_image_of_each_prefix:
        print(f"No images found in {directory} with the specified prefix.")


def run_donuts(directory, prefix, save_path):
    path = os.path.join(directory, '')
    image_names = glob.glob(os.path.join(path, f'{prefix}*.fits'))
    image_names = sorted(image_names)
    if not image_names:
        print(f"No images found for prefix: {prefix}")
        return
    reference_image_name = image_names[0]
    print(f"Using {reference_image_name} as reference for prefix: {prefix}\n")
    science_image_names = image_names[1:]  # Exclude the reference image
    d = Donuts(
        refimage=reference_image_name,
        image_ext=0,
        overscan_width=20,
        prescan_width=20,
        border=64,
        normalise=True,
        exposure='EXPTIME',
        subtract_bkg=True,
        ntiles=32)
    x_shifts = []
    y_shifts = []
    times = []
    for image in science_image_names:
        shift_result = d.measure_shift(image)
        x = shift_result.x.value
        y = shift_result.y.value
        # Get JD time
        with fits.open(image) as hdulist:
            header = hdulist[0].header
            utc_time_str = header.get('DATE-OBS')
            if utc_time_str:
                jd = utc_to_jd(utc_time_str)
                times.append(jd)
            else:
                print(f"DATE-OBS not found in {image}")
                continue
        if abs(x) < 0.5 and abs(y) < 0.5:
            print(f"Image {os.path.basename(image)} shifts (x, y): {x}, {y}")
        else:
            print(f"WARNING: Image {os.path.basename(image)} is not aligned "
                  f"with shifts (x, y): {x}, {y}")
        x_shifts.append(x)
        y_shifts.append(y)
    num_large_shifts = sum(1 for x, y in zip(x_shifts, y_shifts)
                           if abs(x) >= 0.5 or abs(y) >= 0.5)
    print(f"The number of images with shifts greater than 0.5 pixels is: "
          f"{num_large_shifts}\n")
    plot_shifts(x_shifts, y_shifts, times, save_path, prefix)


def plot_shifts(x_shifts, y_shifts, times, save_path, prefix):
    """
    Plot X and Y shifts as a function of time (JD).

    Parameters
    ----------
    x_shifts :
        list of X-shift values for each image.
    y_shifts :
        List of Y-shift values for each image.
    times :
        List of time values (JD) from image headers.
    save_path : str
        Directory where the plot should be saved.
    prefix : str
        Prefix to identify the dataset.
    """
    # Plot X and Y shifts as a function of time
    fig, ax = plt.subplots(figsize=(8, 6))

    # Remove a constant from the time values
    time = [t - 2460551 for t in times]
    scatter = ax.scatter(x_shifts, y_shifts, c=time, cmap='viridis',
                         label='Shifts for field: {}'.format(prefix), marker='o')
    plt.xlabel('X Shift (pixels)')
    plt.ylabel('Y Shift (pixels)')
    plt.title('Shifts to ref image')
    plt.axhline(0, color='black', linestyle='-', linewidth=1)  # Add horizontal line at y=0
    plt.axvline(0, color='black', linestyle='-', linewidth=1)  # Add vertical line at x=0
    plt.legend()

    # Set the axes limits to center (0, 0)
    plt.xlim(-1, 1)
    plt.ylim(-1, 1)

    # Add colorbar
    cbar = plt.colorbar(scatter, label='Time')

    # Save the plot
    pdf_file_path = os.path.join(save_path, f"{prefix}_shifts.pdf")
    fig.savefig(pdf_file_path, bbox_inches='tight')
    print(f"Plot saved to: {pdf_file_path}\n")
    plt.show()


def main():
    # Define base paths
    base_path_1 = '/Users/u5500483/Downloads/DATA_MAC/CMOS/'
    base_path_2 = '/home/ops/data/'
    # Check which base path exists
    if os.path.exists(base_path_1):
        base_path = base_path_1
    else:
        base_path = base_path_2
    save_path = os.path.join(base_path, 'shifts_plots')
    if not os.path.exists(save_path):
        os.makedirs(save_path)
    # Use current directory as the data directory
    current_directory = '/Users/u5500483/Downloads/DATA_MAC/CMOS/20240828/'
    if os.path.exists(current_directory):
        print(f"Processing images in: {current_directory}")
        find_first_image_of_each_prefix(current_directory, save_path)
    else:
        print("Data directory not found.")
        sys.exit(1)


if __name__ == "__main__":
    main()