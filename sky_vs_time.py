#!/usr/bin/env python
import numpy as np
import os
import matplotlib.pyplot as plt
from astropy.io import fits
from matplotlib import gridspec
from scipy.interpolate import interp1d
from astropy.time import Time
import matplotlib.dates as mdates  # Get tick locations from ax0 in datetime, convert to JD
from matplotlib.dates import DateFormatter
from plot_images import plot_images, bin_time_flux_error

plot_images()


def get_phot_files():
    """
    Get photometry files for CMOS and CCD for full and new moon
    """
    path = '/Users/u5500483/Downloads/'

    phot_files = {"CMOS": None, "CCD": None, "CMOS_FULL": None, "CCD_FULL": None}

    for filename in os.listdir(path):
        if filename.startswith('phot') and filename.endswith('CMOS_FULL.fits'):
            phot_files["CMOS_FULL"] = os.path.join(path, filename)

    for filename in os.listdir(path):
        if filename.startswith('phot') and filename.endswith('CCD_FULL.fits'):
            phot_files["CCD_FULL"] = os.path.join(path, filename)

    for filename in os.listdir(path):
        if filename.startswith('phot') and filename.endswith('CMOS.fits'):
            phot_files["CMOS"] = os.path.join(path, filename)

    for filename in os.listdir(path):
        if filename.startswith('phot') and filename.endswith('CCD.fits'):
            phot_files["CCD"] = os.path.join(path, filename)

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
        airmass = phot_data['airmass'][mask]

        # Sort by time
        sorted_indices = np.argsort(times)
        times = times[sorted_indices]  # Convert to BJD_TDB
        airmass = airmass[sorted_indices]  # sort it the same way
        from astropy.time import Time
        times_jd = times  # Keep original JD
        times_utc = Time(times_jd, format='jd').to_datetime()

        sky_background = sky_background[sorted_indices]

        # Normalize sky background per arcsec^2
        sky_background = ((sky_background * gain) / exposure) / area

        return times_jd, times_utc, sky_background, airmass


def main():
    target_tic_id = 214662790

    # Get photometry files
    phot_files = get_phot_files()
    print(f"CMOS file: {phot_files['CMOS']}")
    print(f"CCD file: {phot_files['CCD']}")
    print(f"CMOS_FULL file: {phot_files['CMOS_FULL']}")
    print(f"CCD_FULL file: {phot_files['CCD_FULL']}")

    # Set parameters
    CMOS_PARAMS = {"aperture": 5, "gain": 1.13, "exposure": 10.0, "area": (400 * np.pi)}
    CCD_PARAMS = {"aperture": 4, "gain": 2.0, "exposure": 10.0, "area": (400 * np.pi)}

    times_cmos_jd, times_cmos_utc, sky_cmos, airmass_cmos = (
        extract_sky_background(phot_files["CMOS"], target_tic_id, **CMOS_PARAMS))
    times_ccd_jd, times_ccd_utc, sky_ccd, airmass_ccd = (
        extract_sky_background(phot_files["CCD"], target_tic_id, **CCD_PARAMS))

    times_cmos_jd_FULL, times_cmos_utc_FULL, sky_cmos_FULL, airmass_cmos_FULL = (
        extract_sky_background(phot_files["CMOS_FULL"], target_tic_id, **CMOS_PARAMS))
    times_ccd_jd_FULL, times_ccd_utc_FULL, sky_ccd_FULL, airmass_ccd_FULL = (
        extract_sky_background(phot_files["CCD_FULL"], target_tic_id, **CCD_PARAMS))

    # Binning
    binned_time_cmos, binned_sky_cmos, _ = bin_time_flux_error(times_cmos_jd, sky_cmos, np.zeros_like(sky_cmos), 30)
    binned_time_ccd, binned_sky_ccd, _ = bin_time_flux_error(times_ccd_jd, sky_ccd, np.zeros_like(sky_ccd), 30)
    binned_time_cmos_FULL, binned_sky_cmos_FULL, _ = bin_time_flux_error(times_cmos_jd_FULL, sky_cmos_FULL,
                                                                         np.zeros_like(sky_cmos_FULL), 30)
    binned_time_ccd_FULL, binned_sky_ccd_FULL, _ = bin_time_flux_error(times_ccd_jd_FULL, sky_ccd_FULL,
                                                                       np.zeros_like(sky_ccd_FULL), 30)

    binned_time_cmos_utc = Time(binned_time_cmos, format='jd').to_datetime()
    binned_time_ccd_utc = Time(binned_time_ccd, format='jd').to_datetime()
    binned_time_cmos_FULL_utc = Time(binned_time_cmos_FULL, format='jd').to_datetime()
    binned_time_ccd_FULL_utc = Time(binned_time_ccd_FULL, format='jd').to_datetime()

    # Ensure both binned datasets have the same length
    min_length = min(len(binned_time_cmos), len(binned_time_ccd))
    min_length_full = min(len(binned_time_cmos_FULL), len(binned_time_ccd_FULL))

    # Interpolate to match the shorter dataset length
    if len(binned_time_cmos) > min_length:
        interp_func_cmos = interp1d(binned_time_cmos, binned_sky_cmos, kind='linear', fill_value="extrapolate")
        binned_time_cmos = np.linspace(binned_time_cmos[0], binned_time_cmos[-1], min_length)
        binned_sky_cmos = interp_func_cmos(binned_time_cmos)
        binned_time_cmos_utc = Time(binned_time_cmos, format='jd').to_datetime()

    if len(binned_time_ccd) > min_length:
        interp_func_ccd = interp1d(binned_time_ccd, binned_sky_ccd, kind='linear', fill_value="extrapolate")
        binned_time_ccd = np.linspace(binned_time_ccd[0], binned_time_ccd[-1], min_length)
        binned_sky_ccd = interp_func_cmos(binned_time_ccd)
        binned_time_ccd_utc = Time(binned_time_ccd, format='jd').to_datetime()

    if len(binned_time_cmos_FULL) > min_length_full:
        interp_func_cmos_FULL = interp1d(binned_time_cmos_FULL, binned_sky_cmos_FULL, kind='linear', fill_value="extrapolate")
        binned_time_cmos_FULL = np.linspace(binned_time_cmos_FULL[0], binned_time_cmos_FULL[-1], min_length_full)
        binned_sky_cmos_FULL = interp_func_cmos_FULL(binned_time_cmos_FULL)
        binned_time_cmos_FULL_utc = Time(binned_time_cmos_FULL, format='jd').to_datetime()

    if len(binned_time_ccd_FULL) > min_length_full:
        interp_func_ccd_FULL = interp1d(binned_time_ccd_FULL, binned_sky_ccd_FULL, kind='linear', fill_value="extrapolate")
        binned_time_ccd_FULL = np.linspace(binned_time_ccd_FULL[0], binned_time_ccd_FULL[-1], min_length_full)
        binned_sky_ccd_FULL = interp_func_cmos_FULL(binned_time_ccd_FULL)
        binned_time_ccd_FULL_utc = Time(binned_time_ccd_FULL, format='jd').to_datetime()

    print(f'Final Binned Data Length: CMOS = {len(binned_time_cmos)}, CCD = {len(binned_time_ccd)}')
    print(f'Final Binned Data Length: CMOS_FULL = {len(binned_time_cmos_FULL)}, CCD_FULL = {len(binned_time_ccd_FULL)}')

    # Create subplots: 1 column, 2 rows
    fig = plt.figure(figsize=(12, 8))
    gs = gridspec.GridSpec(2, 2, height_ratios=[1, 1])

    # --- First Plot: Sky background ---
    ax0 = fig.add_subplot(gs[0])
    if times_cmos_utc is not None:
        ax0.plot(times_cmos_utc, sky_cmos, 'bo-', label='CMOS', alpha=0.7)
    if times_ccd_utc is not None:
        ax0.plot(times_ccd_utc, sky_ccd, 'ro-', label='CCD', alpha=0.7)

    ax0.set_ylabel(r'Sky background ($\mathdefault{e^{-}s^{-1}\,arcsec^{-2}}$)')
    ax0.tick_params(axis='x', which='both', labelbottom=False, bottom=False)
    ax0.axvline(Time('2024-07-06 10:05:00', format='iso').to_datetime(), color='black', linestyle='--', linewidth=2)
    ax0.grid()

    # --- Airmass as secondary x-axis (top) ---
    ax2 = ax0.twiny()
    ax2.set_xlim(ax0.get_xlim())
    tick_datetimes = ax0.get_xticks()  # These are in Matplotlib date format
    tick_dates = [mdates.num2date(tick) for tick in tick_datetimes]  # Convert to datetime
    tick_jds = Time(tick_dates).jd  # Convert datetime to JD

    # Prepare airmass interpolator using JD
    unique_jd_mid = np.unique(times_cmos_jd)
    unique_airmass = [airmass_cmos[np.where(times_cmos_jd == jd)[0][0]] for jd in unique_jd_mid]
    airmass_interp = interp1d(unique_jd_mid, unique_airmass, fill_value="extrapolate")

    # Interpolate airmass at tick JD positions
    airmass_ticks = airmass_interp(tick_jds)

    # Set airmass labels on the secondary x-axis
    ax2.set_xticks(tick_datetimes)
    ax2.set_xticklabels([f"{am:.2f}" for am in airmass_ticks])
    ax2.set_xlabel('Airmass')

    # --- Second Plot: CMOS/CCD ratio ---
    ax1 = fig.add_subplot(gs[2], sharex=ax0)
    ax1.plot(binned_time_cmos_utc, binned_sky_cmos / binned_sky_ccd, 's-', color='purple', label='CMOS/CCD',
                alpha=0.7)
    ax1.axvline(Time('2024-07-06 10:05:00', format='iso').to_datetime(), color='black', linestyle='--', linewidth=2)
    ax1.set_xlabel('UTC Time (20240706)')
    ax1.set_ylabel('CMOS/CCD Ratio')
    ax1.grid()

    # Format shared x-axis as HH:MM:SS only on second plot
    ax1.xaxis.set_major_formatter(DateFormatter('%H:%M:%S'))
    ax1.tick_params(axis='x', rotation=45)

    # --- Fourth Plot: Sky background for full moon ---
    ax4 = fig.add_subplot(gs[1])
    if times_cmos_utc_FULL is not None:
        ax4.plot(times_cmos_utc_FULL, sky_cmos_FULL, 'bo-', label='CMOS', alpha=0.7)
    if times_ccd_utc_FULL is not None:
        ax4.plot(times_ccd_utc_FULL, sky_ccd_FULL, 'ro-', label='CCD', alpha=0.7)
    ax4.tick_params(axis='x', which='both', labelbottom=False, bottom=False)
    ax4.axvline(Time('2024-06-23 10:03:00', format='iso').to_datetime(), color='black', linestyle='--', linewidth=2)
    ax4.grid()

    # --- Third Plot: CMOS/CCD ratio for full moon ---
    ax3 = fig.add_subplot(gs[3], sharex=ax4)
    ax3.plot(binned_time_cmos_FULL_utc, binned_sky_cmos_FULL / binned_sky_ccd_FULL, 's-', color='purple', label='CMOS/CCD',
                alpha=0.7)
    ax3.axvline(Time('2024-06-23 10:03:00', format='iso').to_datetime(), color='black', linestyle='--', linewidth=2)
    ax3.set_xlabel('UTC Time (20240623)')
    ax3.grid()
    ax3.xaxis.set_major_formatter(DateFormatter('%H:%M:%S'))
    ax3.tick_params(axis='x', rotation=45)

    ax1.text(
        0.02, 0.03, 'No-moon', transform=ax1.transAxes,
        fontsize=15, fontweight='bold', ha='left', va='bottom',
        bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3')
    )
    ax2.text(
        0.02, 0.03, 'No-moon', transform=ax2.transAxes,
        fontsize=15, fontweight='bold', ha='left', va='bottom',
        bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3')
    )
    ax3.text(
        0.02, 0.03, 'Full-moon', transform=ax3.transAxes,
        fontsize=15, fontweight='bold', ha='left', va='bottom',
        bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3')
    )
    ax4.text(
        0.02, 0.03, 'Full-moon', transform=ax4.transAxes,
        fontsize=15, fontweight='bold', ha='left', va='bottom',
        bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3')
    )

    # --- Airmass as secondary x-axis (top) ---
    ax5 = ax4.twiny()
    ax5.set_xlim(ax4.get_xlim())
    tick_datetimes_FULL = ax4.get_xticks()  # These are in Matplotlib date format
    tick_dates_FULL = [mdates.num2date(tick) for tick in tick_datetimes_FULL]  # Convert to datetime
    tick_jds_FULL = Time(tick_dates_FULL).jd  # Convert datetime to JD

    # Prepare airmass interpolator using JD
    unique_jd_mid_FULL = np.unique(times_cmos_jd_FULL)
    unique_airmass_FULL = [airmass_cmos_FULL[np.where(times_cmos_jd_FULL == jd)[0][0]] for jd in unique_jd_mid_FULL]
    airmass_interp_FULL = interp1d(unique_jd_mid_FULL, unique_airmass_FULL, fill_value="extrapolate")

    # Interpolate airmass at tick JD positions
    airmass_ticks_FULL = airmass_interp_FULL(tick_jds_FULL)

    # Set airmass labels on the secondary x-axis
    ax5.set_xticks(tick_datetimes_FULL)
    ax5.set_xticklabels([f"{am:.2f}" for am in airmass_ticks_FULL])
    ax5.set_xlabel('Airmass')

    # Save and show
    save_path = '/Users/u5500483/Downloads/'
    plt.tight_layout()
    plt.savefig(save_path + 'sky_vs_time_combined.pdf', dpi=300, bbox_inches='tight')
    plt.show()


if __name__ == "__main__":
    main()