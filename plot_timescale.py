#!/usr/bin/env python
import json
from matplotlib import pyplot as plt, ticker
from plot_images import plot_images

plt.rc('font', family='Times New Roman')
plot_images()


# Load the JSON data
def load_json(file_path):
    with open(file_path, 'r') as file:
        return json.load(file)


def extract_data(file_path):
    """
    Extract times, RMS values, and models for both CMOS and CCD.
    """
    data = load_json(file_path)

    # Extract data for CMOS
    times_CMOS = data['file1']['times']
    avg_rms_CMOS = data['file1']['avg_rms']
    RMS_model_CMOS = data['file1']['rms_model']

    # Extract data for CCD
    times_CCD = data['file2']['times']
    avg_rms_CCD = data['file2']['avg_rms']
    RMS_model_CCD = data['file2']['rms_model']

    return times_CMOS, avg_rms_CMOS, RMS_model_CMOS, times_CCD, avg_rms_CCD, RMS_model_CCD


def main():
    # File paths
    file_paths = [
        'files/rms_vs_time_0705_10-11.json',
        'files/rms_vs_time_0705_11-12.json',
        'files/rms_vs_time_0705_12-13.json',
        'files/rms_vs_time_0705_13-14.json',
    ]
    labels = [r"$\mathdefault{10<\mathdefault{T}<11}$", r"$\mathdefault{11<\mathdefault{T}<12}$",
              r"$\mathdefault{12<\mathdefault{T}<13}$", r"$\mathdefault{13<\mathdefault{T}<14}$"]

    # Create a single figure with 1 row and 4 columns
    fig, axes = plt.subplots(1, 4, figsize=(12, 8), sharey=True, constrained_layout=True)

    for i, file_path in enumerate(file_paths):
        times_CMOS, avg_rms_CMOS, RMS_model_CMOS, times_CCD, avg_rms_CCD, RMS_model_CCD = extract_data(file_path)

        # Plot for CMOS
        axes[i].plot(times_CMOS, avg_rms_CMOS, 'o', color='blue', label='CMOS Data')
        axes[i].plot(times_CMOS, RMS_model_CMOS, '--', color='blue', label='CMOS Model')

        # Plot for CCD
        axes[i].plot(times_CCD, avg_rms_CCD, 'o', color='red', label='CCD Data')
        axes[i].plot(times_CCD, RMS_model_CCD, '--', color='red', label='CCD Model')

        # Add a text box in the top right corner
        axes[i].text(
            0.95, 0.98, labels[i], transform=axes[i].transAxes,
            fontsize=15, verticalalignment='top', horizontalalignment='right',
            bbox=dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3')
        )

        # Set logarithmic scales
        axes[i].set_xscale('log')
        axes[i].set_yscale('log')

        axes[i].set_xlabel('Exposure Time (s)')
        if i == 0:
            axes[i].set_ylabel('RMS (ppm)')

        # Format the y-axis tick labels
        axes[i].yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=False))
        axes[i].yaxis.set_minor_formatter(ticker.ScalarFormatter(useMathText=False))
        axes[i].tick_params(axis='y', which='minor', length=4)

    # Save and show the plot
    plt.savefig('RMS_vs_Time.pdf', dpi=300)
    plt.show()


if __name__ == "__main__":
    main()