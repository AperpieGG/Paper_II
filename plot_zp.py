import matplotlib.pyplot as plt
import numpy as np

from plot_images import plot_images, load_json


def main():
    plot_images()

    # Load the data from the target JSON file
    zp_list = load_json('zp5_list.json')

    # plot zp_list on a histogram
    plt.hist(zp_list, bins=100, label=f"Avg: {np.nanmean(zp_list):.2f}")
    plt.axvline(np.nanmean(zp_list), color='black', linestyle='--', label=f"Mean: {np.nanmean(zp_list):.2f}")
    # plt.axvline(np.nanmode(zp_list), color='g', linestyle='--', label=f"Median: {np.nanmedian(zp_list):.2f}")
    plt.xlabel('Zero Point')
    plt.ylabel('Frequency')
    # plt.yscale('log')
    # plt.legend(loc='upper right')
    plt.xlim(18, 21)
    plt.savefig(f'zp5.pdf', dpi=300)
    print(np.mean(zp_list))
    plt.show()


main()
