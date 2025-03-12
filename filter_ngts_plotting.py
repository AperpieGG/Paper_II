import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import gridspec
from scipy.interpolate import UnivariateSpline, interp1d
import numpy as np
from plot_images import plot_images

plot_images()
path = '/Users/u5500483/Documents/GitHub/Paper_II/files/'
# Load the CSV file with headers on the second line
file_path = path + 'goto_filters.csv'  # Replace with the actual file path
data = pd.read_csv(file_path, header=1)
print(data.head())

qe_ccd = data['QE_CCD']  # This is the quantum efficiency of the CCD
wv_ccd = data['WV_CCD']  # This is the wavelength of the CCD
ngts_tra = data['NGTS_TRA']  # This is the transmission of the NGTS
wv_ngts = data['WV_NGTS']  # This is the wavelength of the NGTS


def smoothen_curve(x, y, smooth):
    spline = UnivariateSpline(x, y, s=smooth)

    smooth_x = np.linspace(x.min(), x.max(), 1000)
    smooth_y = spline(smooth_x)

    return smooth_x, smooth_y


# Ensure no NaN values in wv_cmos and qe_cmos
data = data.dropna(subset=['WV_CMOS', 'QE_CMOS'])

# Sort wv_cmos and qe_cmos in ascending order
data = data.sort_values(by='WV_CMOS')

# Extract sorted wv_cmos and qe_cmos
wv_cmos = data['WV_CMOS']
qe_cmos = data['QE_CMOS']

# Smooth the curve after sorting
wv_cmos, qe_cmos = smoothen_curve(wv_cmos, qe_cmos, 400)

# Interpolate CCD QE data onto the CMOS wavelength grid
interp_ccd = interp1d(wv_ccd, qe_ccd, kind='linear', bounds_error=False, fill_value=np.nan)
qe_ccd_interp = interp_ccd(wv_cmos)

# Compute QE ratio (CMOS/CCD)
qe_ratio = qe_cmos / qe_ccd_interp

# Create figure with GridSpec for subplot arrangement
fig = plt.figure(figsize=(6, 6))
gs = gridspec.GridSpec(2, 1, height_ratios=[3, 1])

# First subplot: Quantum Efficiency of sCMOS and CCD
ax0 = fig.add_subplot(gs[0])
ax0.plot(wv_cmos, qe_cmos, label="sCMOS", color="blue")
ax0.plot(wv_ccd, qe_ccd, label="CCD", color="red", linestyle="--")
ax0.fill_between(wv_ngts, 0, 100, where=(ngts_tra > 90), color="grey", alpha=0.5)

ax0.set_ylabel("Quantum Efficiency %")
ax0.set_xlim(400, 1000)
ax0.set_ylim(0, 100)
ax0.grid(True)
ax0.legend(loc='upper center', bbox_to_anchor=(0.5, 1.12), ncol=5, frameon=False, fontsize=16)
# Second subplot: QE Ratio CMOS/CCD
ax1 = fig.add_subplot(gs[1], sharex=ax0)
ax1.plot(wv_cmos[qe_ratio >= 1], qe_ratio[qe_ratio >= 1], color="blue")
ax1.plot(wv_cmos[qe_ratio < 1], qe_ratio[qe_ratio < 1], color="red")
ax1.axhline(y=1, color='black', linestyle=':', linewidth=1)  # Reference line at QE Ratio = 1
print(f'The max and min QE ratio is: {np.max(qe_ratio)} and {np.min(qe_ratio)} at {wv_cmos[np.argmax(qe_ratio)]} and {wv_cmos[np.argmin(qe_ratio)]} nm respectively')
ax1.set_xlabel("Wavelength (nm)")
ax1.set_ylabel("QE Ratio")
ax1.set_xlim(400, 1000)
ax1.fill_between(wv_ngts, 0, 100, where=(ngts_tra > 90), color="grey", alpha=0.5)
ax1.set_ylim(0, 2)
ax1.grid(True)

# Adjust layout and show the plot
plt.tight_layout()
save_path = '/Users/u5500483/Downloads/'
plt.savefig(save_path + 'filter_analysis.pdf', dpi=300)
plt.show()
