#!/usr/bin/env python
import json
import numpy as np
import matplotlib.pyplot as plt
from astropy.stats import sigma_clip
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
    cmos_file = 'files/flux_vs_temperature_CMOS_0705.json'
    ccd_file = 'files/flux_vs_temperature_CCD_0705.json'

    cmos_file_sky = 'files/sky_flux_vs_temperature_CMOS_0705'
    ccd_file_sky = 'files/sky_flux_vs_temperature_CCD_0705.json'

    print(f"Loading CMOS data from {cmos_file}...")
    cmos_data = load_json(cmos_file)
    print(f"Loading CCD data from {ccd_file}...")
    ccd_data = load_json(ccd_file)

    # Create dictionaries keyed by TIC_ID for quick lookups
    cmos_dict = {entry["TIC_ID"]: entry for entry in cmos_data}
    ccd_dict = {entry["TIC_ID"]: entry for entry in ccd_data}

    # Match TIC_IDs and calculate the flux ratio (CMOS/CCD)
    temperatures = []
    flux_ratios = []
    tmags = []
    colors = []

    for tic_id, cmos_entry in cmos_dict.items():
        if tic_id in ccd_dict:
            ccd_entry = ccd_dict[tic_id]

            # Exclude stars with Tmag > 12.5
            if cmos_entry["Tmag"] > 12:
                continue

            # Exclude stars with color > 8000
            if cmos_entry["Teff"] > 8000:
                continue

            # Calculate the flux ratio
            cmos_flux = cmos_entry["Converted_Flux"]
            ccd_flux = ccd_entry["Converted_Flux"]
            flux_ratio = cmos_flux / ccd_flux

            # Append the data
            temperatures.append(cmos_entry["Teff"])
            flux_ratios.append(flux_ratio)
            tmags.append(cmos_entry["Tmag"])  # Use Tmag from CMOS (assuming consistent)
            colors.append(cmos_entry["COLOR"])

    # Convert to numpy arrays for easier handling
    temperatures = np.array(temperatures)
    flux_ratios = np.array(flux_ratios)
    tmags = np.array(tmags)
    colors = np.array(colors)

    # Perform sigma clipping
    # for i in range(3):  # Perform 3 iterations of sigma clipping
    #     clipped_indices = sigma_clip(flux_ratios, sigma=3, cenfunc='mean', stdfunc='std').mask
    #     temperatures = temperatures[~clipped_indices]
    #     flux_ratios = flux_ratios[~clipped_indices]
    #     tmags = tmags[~clipped_indices]
    #     colors = colors[~clipped_indices]

    # Fit a degree-4 polynomial
    polynomial_coeffs = np.polyfit(temperatures, flux_ratios, 4)
    polynomial_fit = np.poly1d(polynomial_coeffs)
    fitted_values = polynomial_fit(temperatures)

    # Print the polynomial function in readable form
    poly_terms = [f"{coeff:.3e}x^{i}" for i, coeff in enumerate(polynomial_coeffs[::-1])]
    poly_function = " + ".join(poly_terms)
    print(f"Fitted Polynomial Function: f(x) = {poly_function}")

    # Plotting
    print("Creating the plot...")
    plt.figure()

    # Scatter plot with temperatures and flux ratios
    scatter = plt.scatter(
        temperatures, flux_ratios, c=tmags, cmap='hot_r', edgecolor='k', alpha=1)

    # # Add the polynomial fit
    # plt.plot(
    #     np.sort(temperatures),
    #     polynomial_fit(np.sort(temperatures)),
    #     color='black',
    #     linestyle='--',
    #     label='Polynomial Fit (Degree 4)'
    # )

    # Add colorbar
    cbar = plt.colorbar(scatter, label='TESS Magnitude')
    # cbar.set_label('$\mathdefault{G_{BP}-G_{RP}}$')
    plt.xlabel('Teff (K)')
    plt.ylabel('CMOS/CCD Flux Ratio')
    # plt.ylim(0.8, 1.6)
    plt.ylim(0.51, 1.75)
    plt.tight_layout()
    plt.savefig('ratio_flux_0705.pdf', dpi=300)
    plt.show()


if __name__ == "__main__":
    main()