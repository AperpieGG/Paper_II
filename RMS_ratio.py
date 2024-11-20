#!/usr/bin/env python
import argparse
import json
import matplotlib.pyplot as plt
import numpy as np
from plot_images import plot_images


def load_json_files():
    """
    Load CCD and CMOS JSON files for analysis.
    Returns
    -------
    tuple
        A tuple containing data from the two JSON files (CCD and CMOS).
    """
    ccd_file_bright = 'files/rms_mags_rel_phot_NG2320-1302_4_1_0622.json'
    cmos_file_bright = 'files/rms_mags_rel_phot_NG2320-1302_5_1_0622.json'
    ccd_file_dark = 'files/rms_mags_rel_phot_NG2320-1302_4_1_0705.json'
    cmos_file_dark = 'files/rms_mags_rel_phot_NG2320-1302_5_1_0705.json'
    with open(ccd_file_bright, 'r') as f1, open(ccd_file_dark, 'r') as f2, open(cmos_file_bright, 'r') as f3, open(cmos_file_dark, 'r') as f4:
        data_ccd_bright = json.load(f1)
        data_ccd_dark = json.load(f2)
        data_cmos_bright = json.load(f3)
        data_cmos_dark = json.load(f4)

    return data_ccd_bright, data_ccd_dark, data_cmos_bright, data_cmos_dark


def compute_rms_ratios(data_ccd, data_cmos):
    # Extract fields as lists
    tic_ids_ccd = data_ccd["TIC_IDs"]
    rms_ccd = data_ccd["RMS_list"]
    tmag_ccd = data_ccd["Tmag_list"]
    color = data_ccd["COLOR"]

    tic_ids_cmos = data_cmos["TIC_IDs"]
    rms_cmos = data_cmos["RMS_list"]

    if len(tic_ids_ccd) != len(tic_ids_cmos) or len(rms_ccd) != len(rms_cmos):
        raise ValueError("Mismatched data lengths between JSON files. Ensure both files have the same TIC_IDs.")

    # Compute RMS ratio and collect Tmag and color values
    rms_ratio = []
    tmag_values = []
    color_values = []

    for i in range(len(tic_ids_ccd)):
        if rms_cmos[i] != 0:  # Avoid division by zero
            ratio = rms_ccd[i] / rms_cmos[i]
            rms_ratio.append(ratio)
            tmag_values.append(tmag_ccd[i])
            color_values.append(color[i])

    return tmag_values, rms_ratio, color_values


def plot_comparison(data_ccd_bright, data_ccd_dark, data_cmos_bright, data_cmos_dark):
    fig, axs = plt.subplots(2, 2, figsize=(12, 8), gridspec_kw={'height_ratios': [2, 1]})

    # Compute ratios for bright and dark nights
    tmag_bright, rms_ratio_bright, color_bright = compute_rms_ratios(data_ccd_bright, data_cmos_bright)
    tmag_dark, rms_ratio_dark, color_dark = compute_rms_ratios(data_ccd_dark, data_cmos_dark)

    # First row: scatter plots for bright and dark nights
    scatter_bright = axs[0, 0].scatter(tmag_bright, rms_ratio_bright, c=color_bright, cmap='coolwarm', vmin=0.5, vmax=1.5)
    axs[0, 0].axhline(y=1, color='black', linestyle='--')
    axs[0, 0].invert_xaxis()
    axs[0, 0].set_ylim(0, 2)
    axs[0, 0].set_xlim(14.1, 9)
    axs[0, 0].set_xlabel('Tmag')
    axs[0, 0].set_ylabel('CCD / CMOS RMS Ratio')

    scatter_dark = axs[0, 1].scatter(tmag_dark, rms_ratio_dark, c=color_dark, cmap='coolwarm', vmin=0.5, vmax=1.5)
    axs[0, 1].axhline(y=1, color='black', linestyle='--')
    axs[0, 1].invert_xaxis()
    axs[0, 1].set_ylim(0, 2)
    axs[0, 1].set_xlim(14.1, 9)
    axs[0, 1].set_xlabel('Tmag')

    # Add colorbars
    cbar_dark = fig.colorbar(scatter_dark, ax=axs[0, 1], orientation='vertical', fraction=0.046, pad=0.04)
    cbar_dark.set_label('$\mathdefault{G_{BP}-G_{RP}}$')
    cbar_bright = fig.colorbar(scatter_dark, ax=axs[0, 0], orientation='vertical', fraction=0.046, pad=0.04)

    bins = 50000  # Adjust bins as necessary for better visibility
    print(f'The mean ratio for bright is: {np.median(rms_ratio_bright)}')
    print(f'The mean ratio for dark is: {np.median(rms_ratio_dark)}')
    # Add mean lines to the histograms
    axs[1, 0].axvline(x=np.median(rms_ratio_bright), color='red', linestyle='--', label='Mean (Bright)')
    axs[1, 1].axvline(x=np.median(rms_ratio_dark), color='red', linestyle='--', label='Mean (Dark)')

    axs[1, 0].hist(rms_ratio_bright, bins=bins)
    axs[1, 0].set_xlabel('CCD / CMOS RMS Ratio')
    axs[1, 0].set_ylabel('Frequency')
    axs[1, 0].axvline(x=1, color='black', linestyle='--')
    axs[1, 0].set_xlim(0, 2)
    axs[1, 0].set_yscale('log')
    axs[1, 0].set_ylim(0, 500)

    axs[1, 1].hist(rms_ratio_dark, bins=bins)
    axs[1, 1].set_xlabel('CCD / CMOS RMS Ratio')
    axs[1, 1].axvline(x=1, color='black', linestyle='--')
    axs[1, 1].set_xlim(0, 2)
    axs[1, 1].set_yscale('log')
    axs[1, 1].set_ylim(0, 500)

    plt.tight_layout()
    plt.savefig('RMS_Ratio.pdf', dpi=300)
    plt.show()


def main():
    plot_images()
    data_ccd_bright, data_ccd_dark, data_cmos_bright, data_cmos_dark = load_json_files()
    plot_comparison(data_ccd_bright, data_ccd_dark, data_cmos_bright, data_cmos_dark)


if __name__ == "__main__":
    main()