#!/usr/bin/env python
import argparse
import json
import matplotlib.pyplot as plt
from plot_images import plot_images


def load_json_files(ccd_file, cmos_file):
    """
    Load CCD and CMOS JSON files for analysis.

    Parameters
    ----------
    ccd_file : str
        Path to the CCD JSON file.
    cmos_file : str
        Path to the CMOS JSON file.

    Returns
    -------
    tuple
        A tuple containing data from the two JSON files (CCD and CMOS).
    """
    with open(ccd_file, 'r') as f1, open(cmos_file, 'r') as f2:
        print(f"Loading {ccd_file} and {cmos_file}")
        data1 = json.load(f1)
        data2 = json.load(f2)

    return data1, data2


def compute_rms_ratios(data1, data2):
    # Extract fields as lists
    tic_ids1 = data1["TIC_IDs"]
    rms1 = data1["RMS_list"]
    tmag1 = data1["Tmag_list"]
    color = data1["COLOR"]

    tic_ids2 = data2["TIC_IDs"]
    rms2 = data2["RMS_list"]

    if len(tic_ids1) != len(tic_ids2) or len(rms1) != len(rms2):
        raise ValueError("Mismatched data lengths between JSON files. Ensure both files have the same TIC_IDs.")

    # Compute RMS ratio and collect Tmag and color values
    rms_ratio = []
    tmag_values = []
    color_values = []

    for i in range(len(tic_ids1)):
        if rms2[i] != 0:  # Avoid division by zero
            ratio = rms1[i] / rms2[i]
            rms_ratio.append(ratio)
            tmag_values.append(tmag1[i])
            color_values.append(color[i])

    return tmag_values, rms_ratio, color_values


def plot_rms_ratio(tmag_values, rms_ratio, color_values):
    plt.figure()
    scatter = plt.scatter(tmag_values, rms_ratio, c=color_values, cmap='coolwarm', vmin=0.5, vmax=1.5)
    plt.colorbar(scatter, label='$\mathdefault{G_{BP}-G_{RP}}$')  # Add colorbar for the COLOR field
    plt.axhline(y=1, color='black', linestyle='--')
    plt.xlabel('Tmag')
    plt.ylabel('CCD / CMOS RMS Ratio')
    plt.grid(True)
    plt.ylim(-0.3, 3)
    plt.gca().invert_xaxis()
    plt.tight_layout()
    plt.show()


def main():
    parser = argparse.ArgumentParser(description="Plot RMS Ratio of CCD and CMOS for a set of stars.")
    parser.add_argument('ccd_file', type=str, help='Path to the CCD JSON file')
    parser.add_argument('cmos_file', type=str, help='Path to the CMOS JSON file')
    args = parser.parse_args()

    data1, data2 = load_json_files(args.ccd_file, args.cmos_file)
    tmag_values, rms_ratio, color_values = compute_rms_ratios(data1, data2)
    plot_rms_ratio(tmag_values, rms_ratio, color_values)


if __name__ == "__main__":
    plot_images()
    main()