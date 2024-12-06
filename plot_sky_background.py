#!/usr/bin/env python
import json
import numpy as np
import matplotlib.pyplot as plt
from plot_images import plot_images

plot_images()


def load_json(file_path):
    """
    Load JSON file and return data.
    """
    with open(file_path, 'r') as file:
        return json.load(file)


def main():
    # Load the JSON data
    sky_cmos_file = 'sky_flux_vs_temperature_CMOS_0705.json'
    sky_ccd_file = 'sky_flux_vs_temperature_CCD_0705.json'

    print(f"Loading Sky CMOS data from {sky_cmos_file}...")
    sky_cmos_data = load_json(sky_cmos_file)
    print(f"Loading Sky CCD data from {sky_ccd_file}...")
    sky_ccd_data = load_json(sky_ccd_file)

    # Extract converted flux values and normalize them by (400 * pi)
    sky_flux_cmos = [entry["Converted_Flux"] / (400 * np.pi) for entry in sky_cmos_data]
    sky_flux_ccd = [entry["Converted_Flux"] / (400 * np.pi) for entry in sky_ccd_data]

    print(f"Extracted {len(sky_flux_cmos)} CMOS flux values.")
    print(f"Extracted {len(sky_flux_ccd)} CCD flux values.")

    # Plot the histograms for sky flux
    plt.figure()
    plt.hist(
        sky_flux_cmos,
        bins=1000,
        color='blue',
        alpha=1,
        label='Sky Flux (CMOS)'
    )
    plt.hist(
        sky_flux_ccd,
        bins=1000,
        color='red',
        alpha=1,
        label='Sky Flux (CCD)'
    )
    plt.xlabel(r'Sky background ($\mathdefault{e^{-}s^{-1}\,arcsec^{-2}}$)')
    plt.ylabel('Frequency')
    plt.xlim(0.5, 2)
    plt.tight_layout()
    plt.savefig('sky_flux_histogram_0705.pdf', dpi=300)
    plt.show()


if __name__ == "__main__":
    main()
