import matplotlib.pyplot as plt
import numpy as np

from plot_images import plot_images, load_json


def main():
    plot_images()

    # Load the data from the target JSON file
    data_5 = load_json('zp5_list.json')
    zp_list_5 = data_5['zp_list']
    color_list_5 = data_5['color_list']

    data_4 = load_json('zp4_list.json')
    zp_list_4 = data_4['zp_list']
    color_list_4 = data_4['color_list']

    # plot zp_list on a histogram
    plt.hist(zp_list_5, bins=100, label=f"Avg: {np.nanmean(zp_list_5):.2f}", color='blue', alpha=0.5)
    plt.hist(zp_list_4, bins=100, label=f"Avg: {np.nanmean(zp_list_4):.2f}", color='red', alpha=0.5)
    plt.axvline(np.nanmean(zp_list_5), color='navy', linestyle='--', label=f"Mean: {np.nanmean(zp_list_5):.2f}")
    plt.axvline(np.nanmean(zp_list_4), color='maroon', linestyle='--', label=f"Mean: {np.nanmean(zp_list_4):.2f}")
    # plt.axvline(np.nanmode(zp_list), color='g', linestyle='--', label=f"Median: {np.nanmedian(zp_list):.2f}")
    plt.xlabel('Zero Point')
    plt.ylabel('Frequency')
    plt.yscale('log')
    # plt.legend(loc='upper right')
    plt.xlim(17, 21)
    plt.savefig(f'zp5.pdf', dpi=300)
    print(np.mean(zp_list_5))
    plt.show()


main()
