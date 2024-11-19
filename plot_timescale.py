import json
import os
from matplotlib import pyplot as plt, ticker
from plot_images import plot_images

plot_images()


def load_rms_mags_data():
    """
    Load RMS and magnitude data from JSON file.
    """
    filename_1 = 'files/rms_vs_timescale_0714.json'
    with open(filename_1, 'r') as file:
        data = json.load(file)

    return data


def plot_rms_timescale(data):
    """
    Plot RMS noise model on a given axis.
    """
    times1 = data['file1']['times']
    avg_rms1 = data['file1']['avg_rms']
    RMS_model1 = data['file1']['rms_model']

    times2 = data['file2']['times']
    avg_rms2 = data['file2']['avg_rms']
    RMS_model2 = data['file2']['rms_model']

    fig, axs = plt.subplots(1, 2, figsize=(6, 6), sharey=True)

    axs[0].plot(times1, avg_rms1, 'o', color='black')
    axs[0].plot(times1, RMS_model1, '--', color='black')
    axs[0].axvline(x=900, color='red', linestyle='-')
    axs[0].set_xscale('log')
    axs[0].set_yscale('log')
    axs[0].set_xlabel('Exposure time (s)')
    axs[0].set_ylabel('RMS (ppm)')

    axs[1].plot(times2, avg_rms2, 'o', color='black')
    axs[1].plot(times2, RMS_model2, '--', color='black')
    axs[1].axvline(x=900, color='red', linestyle='-')
    axs[1].set_xscale('log')
    axs[1].set_xlabel('Exposure time (s)')

    plt.gca().yaxis.set_major_formatter(ticker.ScalarFormatter(useMathText=False))
    plt.gca().yaxis.set_minor_formatter(ticker.ScalarFormatter(useMathText=False))
    plt.gca().tick_params(axis='y', which='minor', length=4)
    plt.tight_layout()
    plt.savefig('rms_vs_timescale.pdf', dpi=300)
    plt.show()


if __name__ == "__main__":
    data = load_rms_mags_data()
    plot_rms_timescale(data)
