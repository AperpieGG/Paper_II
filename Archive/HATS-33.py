import matplotlib.pyplot as plt
import pandas as pd

from plot_images import plot_images, load_json, bin_time_flux_error
import numpy as np
from astropy.timeseries import LombScargle

# Plot the images (if needed)
plot_images()
#
# # Read the .dat file to a list of lists
# with open("./table3.dat", "r") as file:
#     datContent = [line.strip().split() for line in file.readlines()]
#
# # Filter rows where the first column is 'HATS-33' and the instrument (last column) is 'HS'
# filtered_rows_HS = [row for row in datContent if row[0] == "HATS-33" and row[5] == "HS"]
#
# # Filter rows where the first column is 'HATS-33' and the instrument (last column) is 'LCOGT1m/sinistro'
# filtered_rows_LCOGT = [row for row in datContent if row[0] == "HATS-33" and row[5] == "i"]
#
# # Extract time (column 2) and flux (column 3) for both instruments
# time_HS = np.array([float(row[1]) for row in filtered_rows_HS])  # Convert to float
# flux_HS = np.array([float(row[2]) for row in filtered_rows_HS])  # Convert to float
#
# time_LCOGT = np.array([float(row[1]) for row in filtered_rows_LCOGT])  # Convert to float
# flux_LCOGT = np.array([float(row[2]) for row in filtered_rows_LCOGT])  # Convert to float
#
# # Combine both datasets (time and flux)
# time_combined = np.concatenate([time_HS])
# flux_combined = np.concatenate([flux_HS])
#
#
# # Lomb-Scargle Periodogram function to plot and estimate the peak period
# def lomb_scargle_analysis(time, flux):
#     # Create the Lomb-Scargle periodogram
#     frequency, power = LombScargle(time, flux).autopower()
#     # Find the peak frequency (max power)
#     peak_frequency = frequency[np.argmax(power)]
#
#     # Convert frequency to period (in days)
#     print(f'The peak frequency is: {peak_frequency:.4f} 1/day')
#     peak_period = 1 / peak_frequency
#     print(f"Peak period: {peak_period:.4f} days")
#
#     # Plot the Lomb-Scargle periodogram
#     plt.figure(figsize=(10, 6), dpi=100)
#     plt.plot(frequency, power)
#     plt.xlabel('Frequency (1/day)')
#     plt.ylabel('Power')
#     plt.xlim(0,1.2)
#     plt.tight_layout()
#     plt.grid(True)
#
#     return peak_period
#
#
# # Plot the time vs flux for both instruments combined
# plt.figure(figsize=(10, 6), dpi=100)
# plt.plot(time_combined, flux_combined, 'bo', label='HS + LCOGT1m/sinistro')
# plt.xlabel('Time (BJD)', fontsize=12)
# plt.ylabel('Observed Flux (mag)', fontsize=12)
# plt.grid(True)
# plt.legend()
# plt.tight_layout()
# plt.show()
#
#
# def plot_phase_folded_light_curve_double(time, flux, period):
#     # Convert time to phase
#     foldTimes = (time / period) % 1
#
#     # Duplicate for two cycles (0–2)
#     new_flux = np.hstack([flux, flux])
#     new_phase = np.hstack([foldTimes, foldTimes + 1])
#
#     # Sort for cleaner plot
#     sorted_indices = np.argsort(new_phase)
#     new_phase = new_phase[sorted_indices]
#     new_flux = new_flux[sorted_indices]
#
#     # Plot
#     plt.figure(figsize=(10, 6), dpi=100)
#     plt.scatter(new_phase, new_flux, color='blue', s=5, label=f"Period: {period:.4f} days")
#     plt.xlim(0, 1)
#     plt.ylim(-0.02, 0.02)
#     plt.gca().invert_yaxis()  # Invert if working with magnitudes
#     plt.xlabel('Phase')
#     plt.ylabel('Observed Flux (mag)')
#     plt.title('Phase-Folded Light Curve over Two Cycles')
#     plt.legend()
#     plt.tight_layout()
#     plt.show()
#
#
# # Run the period analysis and get the period
# peak_period = lomb_scargle_analysis(time_combined, flux_combined)
#
# # Plot the phase-folded light curve
# plot_phase_folded_light_curve_double(time_combined, flux_combined, peak_period)

file = '/Users/u5500483/Downloads/TIC-360742636_REF_IMAGE_ID_81120240715235614.fits'
from astropy.io import fits


def read_phot_file(filename):
    """
    Read the photometry file.

    Parameters
    ----------
    filename : str
        Photometry file to read.

    Returns
    -------
    astropy.table.table.Table
        Table containing the photometry data.
    """
    # Read the photometry file here using fits or any other appropriate method
    try:
        with fits.open(filename) as ff:
            # Access the data in the photometry file as needed
            tab = ff[3].data
            return tab
    except Exception as e:
        print(f"Error reading photometry file {filename}: {e}")
        return None


data = read_phot_file(file)
# find the headers
print(data.columns)
# take the column FLUX and BJD
flux = data['FLUX_SYSREM']
flux_err = data['FLUX_ERR']

mean_flux = np.mean(flux)
# normalize the flux
flux = flux / mean_flux
time = data['BJD']

flux_err = flux * np.sqrt((flux / flux_err) ** 2)

output_file_path = 'allesfit/NGTS.csv'
detrended_lc_df = pd.DataFrame({
    '# time': time,
    'flux': flux,
    'flux_err': flux_err
})
detrended_lc_df.to_csv(output_file_path, index=False)









