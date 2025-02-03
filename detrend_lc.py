#!/usr/bin/env python
import json
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import argparse
from wotan import transit_mask, flatten
from plot_images import plot_images, bin_time_flux_error, load_json, bin_by_time_interval


def main(period, duration, T0, target_file, save_option):
    # Display any necessary images
    plot_images()

    # Load the data from the target JSON file
    data = load_json(target_file)

    # Extract relevant data
    # 1561 for CCD, 2000 for CMOS for WASP-30
    time = np.array(data['Time_BJD'])
    flux = np.array(data['Relative_Flux'])
    flux_err = np.array(data['Relative_Flux_err'])
    print(f'Time length: {len(time)}, Flux length: {len(flux)}, Flux Error length: {len(flux_err)}')

    # Bin the data, for CMOS is 30 for 5 min and for CCD is 23 (10 + 3 seconds readout)
    # time_binned, flux_binned, fluxerr_binned = bin_time_flux_error(time, flux, flux_err, 4)
    # print(f'Binned Time length: {len(time_binned)}, Binned Flux length: {len(flux_binned)}, '
    #       f'Binned Flux Error length: {len(fluxerr_binned)}')

    time_binned, flux_binned, fluxerr_binned = time, flux, flux_err

    # Generate the transit mask
    mask = transit_mask(time=time_binned, period=period, duration=duration, T0=T0)

    # Apply Wotan flattening with transit masking
    flattened_flux, trend_flux = flatten(
        time_binned, flux_binned, method='lowess',
        window_length=3 * duration, return_trend=True,
        robust=True, mask=mask
    )

    # Restore transit points by combining masked and flattened light curve
    restored_flux = flux_binned / trend_flux

    # Create a figure with side-by-side plots
    fig, axs = plt.subplots(1, 2, figsize=(10, 4), dpi=100, sharex=True, sharey=True)
    # Plot original light curve with trend
    axs[0].plot(time_binned, flux_binned, 'r.', label='Original LC')
    axs[0].plot(time_binned, trend_flux, 'b-', label='Trend (Wotan)')
    axs[0].set_xlabel('Time (BJD)')
    axs[0].set_ylabel('Relative Flux')

    # Plot the detrended light curve with restored transit
    axs[1].plot(time_binned, restored_flux, 'g.', label='Detrended LC')
    axs[1].axhline(1, color='b', linestyle='--', label='Baseline')
    axs[1].set_xlabel('Time (BJD)')

    plt.tight_layout()
    plt.show()

    # Save the detrended light curve to a CSV file if requested
    if save_option.lower() == 'yes':
        output_file_path = 'allesfit/NGTS.csv'
        detrended_lc_df = pd.DataFrame({
            '# time': time_binned,
            'flux': flux_binned,
            'flux_err': fluxerr_binned
        })
        detrended_lc_df.to_csv(output_file_path, index=False)
        print(f"CSV file saved at: {output_file_path}")
    else:
        print("CSV file saving skipped.")

    # Find the lowest points
    lowest_points = restored_flux[restored_flux < 0.98]  # Example threshold

    # Calculate the average of the lowest points
    if len(lowest_points) > 0:
        average_lowest_point = np.mean(lowest_points)
        print(f"Average of Lowest Points: {average_lowest_point:.6f}")
    else:
        print("No points below the specified threshold.")


if __name__ == "__main__":
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Process Light Curve with Wotan")
    parser.add_argument('--period', type=float, required=True, help="Transit period in days")
    parser.add_argument('--duration', type=float, required=True, help="Transit duration in days")
    parser.add_argument('--T0', type=float, required=True, help="Mid-transit time (BJD)")
    parser.add_argument('--target_file', type=str, required=True, help="Path to the target JSON file")
    # parser.add_argument('--bin_fac', type=int, default=30, help="Bin factor for time and flux")
    parser.add_argument('--save_option', type=str, default='no', choices=['yes', 'no'], help="Save CSV file (yes/no)")

    args = parser.parse_args()

    # Execute the main function with parsed arguments
    main(args.period, args.duration, args.T0, args.target_file, args.save_option)

