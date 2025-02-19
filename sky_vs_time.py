#!/usr/bin/env python
import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.io import fits
from scipy.interpolate import interp1d

from plot_images import plot_images, bin_time_flux_error

plot_images()


def get_phot_files():
    """
    Get photometry files for CMOS and CCD.
    """
    path_CMOS = '/Users/u5500483/Downloads/'
    path_CCD = '/Users/u5500483/Downloads/'

    phot_files = {"CMOS": None, "CCD": None}

    for filename in os.listdir(path_CMOS):
        if filename.startswith('phot') and filename.endswith('CMOS.fits'):
            phot_files["CMOS"] = os.path.join(path_CMOS, filename)

    for filename in os.listdir(path_CCD):
        if filename.startswith('phot') and filename.endswith('CCD.fits'):
            phot_files["CCD"] = os.path.join(path_CCD, filename)

    return phot_files


def extract_sky_background(phot_file, tic_id, aperture, gain, exposure, area):
    """
    Extracts time and sky background data for a specific TIC_ID from a photometry FITS file.
    """
    with fits.open(phot_file) as phot_hdul:
        phot_data = phot_hdul[1].data

        # Ensure TIC_ID exists
        if tic_id not in phot_data['TIC_ID']:
            return None, None  # Return empty if TIC_ID not found

        # Filter for the target TIC_ID
        mask = phot_data['TIC_ID'] == tic_id
        times = phot_data['jd_mid'][mask]  # Julian Date mid-exposure
        sky_background = phot_data[f'flux_w_sky_{aperture}'][mask] - phot_data[f'flux_{aperture}'][mask]

        # Sort by time
        sorted_indices = np.argsort(times)
        times = times[sorted_indices]  # Convert to BJD_TDB
        times = times - times[0]  # Convert to days relative to the first observation
        sky_background = sky_background[sorted_indices]

        # Normalize sky background per arcsec^2
        sky_background = ((sky_background * gain) / exposure) / area

        return times, sky_background


def main():
    target_tic_id = 214662790

    # Get photometry files
    phot_files = get_phot_files()
    print(f"CMOS file: {phot_files['CMOS']}")
    print(f"CCD file: {phot_files['CCD']}")

    # Set parameters
    CMOS_PARAMS = {"aperture": 5, "gain": 1.13, "exposure": 10.0, "area": (400 * np.pi)}
    CCD_PARAMS = {"aperture": 4, "gain": 2.0, "exposure": 10.0, "area": (400 * np.pi)}

    # Extract data
    times_cmos, sky_cmos = extract_sky_background(phot_files["CMOS"], target_tic_id, **CMOS_PARAMS)
    print(f'The length for the CMOS data is {len(times_cmos)}')
    times_ccd, sky_ccd = extract_sky_background(phot_files["CCD"], target_tic_id, **CCD_PARAMS)
    print(f'The length for the CCD data is {len(times_ccd)}')
    # Plot results
    plt.figure()

    if times_cmos is not None:
        plt.plot(times_cmos, sky_cmos, 'bo-', label='CMOS', alpha=0.7)

    if times_ccd is not None:
        plt.plot(times_ccd, sky_ccd, 'ro-', label='CCD', alpha=0.7)

    plt.xlabel('Time (days)')
    plt.ylabel(r'Sky background ($\mathdefault{e^{-}s^{-1}\,arcsec^{-2}}$)')
    plt.grid()
    save_path = '/Users/u5500483/Downloads/'
    plt.savefig(save_path + 'sky_vs_time_0705.pdf', dpi=300, bbox_inches='tight')
    plt.show()

    # Binning
    # Binning
    binned_time_cmos, binned_sky_cmos, _ = bin_time_flux_error(times_cmos, sky_cmos, np.zeros_like(sky_cmos), 30)
    binned_time_ccd, binned_sky_ccd, _ = bin_time_flux_error(times_ccd, sky_ccd, np.zeros_like(sky_ccd), 30)

    # Ensure both binned datasets have the same length
    min_length = min(len(binned_time_cmos), len(binned_time_ccd))

    # Interpolate to match the shorter dataset length
    if len(binned_time_cmos) > min_length:
        interp_func_cmos = interp1d(binned_time_cmos, binned_sky_cmos, kind='linear', fill_value="extrapolate")
        binned_time_cmos = np.linspace(binned_time_cmos[0], binned_time_cmos[-1], min_length)
        binned_sky_cmos = interp_func_cmos(binned_time_cmos)

    if len(binned_time_ccd) > min_length:
        interp_func_ccd = interp1d(binned_time_ccd, binned_sky_ccd, kind='linear', fill_value="extrapolate")
        binned_time_ccd = np.linspace(binned_time_ccd[0], binned_time_ccd[-1], min_length)
        binned_sky_ccd = interp_func_ccd(binned_time_ccd)

    print(f'Final Binned Data Length: CMOS = {len(binned_time_cmos)}, CCD = {len(binned_time_ccd)}')

    # Plot binned but ratio of CMOS to CCD
    plt.figure()

    plt.plot(binned_time_cmos, binned_sky_cmos / binned_sky_ccd, 'o-', color='purple', label='CMOS/CCD', alpha=0.7)

    plt.xlabel('Time (days)')
    plt.ylabel('CMOS/CCD Sky background ratio')
    plt.grid()
    plt.savefig(save_path + 'sky_ratio_vs_time_0705.pdf', dpi=300, bbox_inches='tight')
    plt.show()


if __name__ == "__main__":
    main()