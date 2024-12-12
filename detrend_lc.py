import json
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from wotan import transit_mask, flatten
from plot_images import plot_images, bin_time_flux_error, load_json

# Display any necessary images
plot_images()

# Load the data
file_path = '/target_light_curve_402026209_CMOS_0802.json'
data = load_json(file_path)

# Extract relevant data
time = np.array(data['Time_BJD'])
flux = np.array(data['Relative_Flux'])
flux_err = np.array(data['Relative_Flux_err'])

# Bin the data
time_binned, flux_binned, fluxerr_binned = bin_time_flux_error(time, flux, flux_err, 30)

# Automatically create a transit mask using Wotan
period = 1.338231602  # Example period; replace with actual value
duration = 0.0898  # Example transit duration in days; adjust if needed
T0 = 2454823.58592  # Example mid-transit time; adjust as needed

# Generate the transit mask
mask = transit_mask(time=time_binned, period=period, duration=duration, T0=T0)

# Apply Wotan flattening with transit masking
flattened_flux, trend_flux = flatten(
    time_binned, flux_binned, method='cosine',
    window_length=3*duration, return_trend=True,
    robust=True, mask=mask
)

# Restore transit points by combining masked and flattened light curve
restored_flux = flux_binned / trend_flux

# Plot original light curve with trend
plt.figure()
plt.plot(time_binned, flux_binned, 'ro', label='Original Light Curve')
plt.plot(time_binned, trend_flux, 'b-', label='Fitted Trend (Wotan)')
plt.xlabel('Time (BJD)')
plt.ylabel('Relative Flux')
plt.title('Light Curve with Wotan Detrending')
plt.ylim(0.965, 1.02)
plt.tight_layout()
plt.show()

# Plot the detrended light curve with restored transit
plt.figure()
plt.plot(time_binned, restored_flux, 'go', label='Detrended Light Curve (Transit Restored)')
plt.axhline(1, color='b', linestyle='--', label='Baseline')
plt.xlabel('Time (BJD)')
plt.ylabel('Relative Flux')
plt.title('Detrended LC')
plt.ylim(0.965, 1.02)
plt.tight_layout()
plt.show()

# Save the detrended light curve to a CSV file
output_file_path = '/allesfitter/NGTS.csv'

# subtract time by the first time value

# Create a DataFrame with the detrended light curve
detrended_lc_df = pd.DataFrame({
    '# time': time_binned,
    'flux': restored_flux,
    'flux_err': fluxerr_binned
})

# Save to CSV
detrended_lc_df.to_csv(output_file_path, index=False)
print(f"CSV file saved at: {output_file_path}")