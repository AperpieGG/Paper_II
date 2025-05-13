import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import numpy as np
from scipy.integrate import trapz
import argparse
from plot_images import plot_images
plot_images()
# Constants
h = 6.626e-34  # Planck's constant (J·s)
c = 3.0e8  # Speed of light (m/s)
k = 1.38e-23  # Boltzmann's constant (J/K)
telescope_diameter = 0.2  # Telescope diameter in meters
telescope_area = np.pi * (telescope_diameter / 2) ** 2  # Collecting area in m²

# Load data from CSV
data_path = "/Users/u5500483/Documents/GitHub/Paper_II/files/goto_filters.csv"
data = pd.read_csv(data_path, header=1)

# Extract relevant columns
wv_ccd = np.array(data['WV_CCD']) * 1e-9  # Convert nm to meters
qe_ccd = np.array(data['QE_CCD']) / 100  # Convert % to fraction
wv_cmos = np.array(data['WV_CMOS']) * 1e-9
qe_cmos = np.array(data['QE_CMOS']) / 100
wv_ngts = np.array(data['WV_NGTS']) * 1e-9
ngts_tra = np.array(data['NGTS_TRA']) / 100
wv_l_band = np.array(data['Wavelength (nm)']) * 1e-9
l_band_tra = np.array(data['%T_L']) / 100
wv_tel = np.array(data['WV_TEL']) * 1e-9
tel_tra = np.array(data['TEL_TRA']) / 100
wv_paranal = np.array(data['WV_PARANAL']) * 1e-9
paranal_tra = np.array(data['TRAN_PARANAL'])

# Ensure telescope transmission drops to zero below 465 nm
tel_tra = np.where(wv_tel < 476e-9, 0, tel_tra)

# Define wavelength range for interpolation (in meters)
wavelengths = np.linspace(400, 1000, 1000) * 1e-9  # Convert nm to meters

# Interpolation functions
qe_ccd_interp = interp1d(wv_ccd, qe_ccd, bounds_error=False, fill_value=0)(wavelengths)
qe_cmos_interp = interp1d(wv_cmos, qe_cmos, bounds_error=False, fill_value=0)(wavelengths)
ngts_interp = interp1d(wv_ngts, ngts_tra, bounds_error=False, fill_value=0)(wavelengths)
l_band_interp = interp1d(wv_l_band, l_band_tra, bounds_error=False, fill_value=0)(wavelengths)
tel_interp = interp1d(wv_tel, tel_tra, bounds_error=False, fill_value=0)(wavelengths)
atm_interp = interp1d(wv_paranal, paranal_tra, bounds_error=False, fill_value=0)(wavelengths)


# Planck function to compute blackbody spectrum
def blackbody_lambda(wavelength, temperature):
    return (2 * h * c ** 2 / wavelength ** 5) / (np.exp(h * c / (wavelength * k * temperature)) - 1)


# Function to compute electrons per second
def compute_electrons_per_second(temperature, camera, filter_type):
    # Compute Blackbody Spectrum
    stellar_flux = blackbody_lambda(wavelengths, temperature)

    # Select the correct QE curve
    qe_interp = qe_cmos_interp if camera == "CMOS" else qe_ccd_interp

    # Select the correct filter transmission
    filter_interp = ngts_interp if filter_type == "NGTS" else l_band_interp

    # Compute total transmission
    total_transmission = qe_interp * filter_interp * atm_interp * tel_interp

    # Compute Photon Energy
    photon_energy = (h * c) / wavelengths

    # Compute Detected Electrons per Second
    detected_photon_flux = (
            stellar_flux * total_transmission * telescope_area * wavelengths / photon_energy
    )

    # Integrate over Wavelengths to get total e-/s
    electrons_per_second = trapz(detected_photon_flux, wavelengths)

    return electrons_per_second, detected_photon_flux


# Argument Parser
parser = argparse.ArgumentParser(description="Calculate electrons per second for a given camera and filter.")
parser.add_argument("--cam", type=str, required=True, choices=["CMOS", "CCD"], help="Camera type (CMOS or CCD)")
parser.add_argument("--filter", type=str, required=True, choices=["L", "NGTS"], help="Filter type (L or NGTS)")
parser.add_argument("--Teff", type=float, required=True, help="Effective temperature of the star (Kelvin)")
args = parser.parse_args()

# Compute electron count
electrons_per_second, detected_photon_flux = compute_electrons_per_second(args.Teff, args.cam, args.filter)

# Generate Plots
fig, axs = plt.subplots(3, 1, figsize=(10, 14))

# Plot the blackbody spectrum
axs[0].plot(wavelengths * 1e9, blackbody_lambda(wavelengths, args.Teff), label=f"Blackbody Spectrum (T={args.Teff}K)",
            color='black')
axs[0].set_ylabel("Flux (arbitrary units)")
axs[0].set_title("Blackbody Radiation Curve")
axs[0].legend()
axs[0].grid()

# Plot the transmission curves
axs[1].plot(wavelengths * 1e9, qe_ccd_interp if args.cam == "CCD" else qe_cmos_interp, label=f"{args.cam} QE",
            color="blue", linestyle="--")
axs[1].plot(wavelengths * 1e9, ngts_interp if args.filter == "NGTS" else l_band_interp, label=f"{args.filter} Filter",
            color="red", linestyle="--")
axs[1].plot(wv_paranal * 1e9, paranal_tra, label="Atmosphere Transmission", color="green", linestyle="--")
axs[1].plot(wavelengths * 1e9, tel_interp, label="Telescope", color="purple", linestyle="--")
axs[1].set_ylabel("Transmission")
axs[1].set_title("Transmission Factors")
axs[1].legend()
axs[1].grid()

# Plot detected photon flux
axs[2].plot(wavelengths * 1e9, detected_photon_flux, label="Detected Photons per Second", color="orange")
axs[2].set_xlabel("Wavelength (nm)")
axs[2].set_ylabel("Detected Photon Rate")
axs[2].set_title("Final Detected Photon Flux")
axs[2].legend()
axs[2].grid()
plt.savefig('/Users/u5500483/Downloads/filter_analysis.pdf', dpi=300)
plt.tight_layout()
plt.show()

# Print the final estimated electron count
print(
    f"Electrons per second detected for Teff={args.Teff}K, Camera={args.cam}, Filter={args.filter}: {electrons_per_second:.2e} e-/s")
