import matplotlib.pyplot as plt
import numpy as np

from plot_images import plot_images, load_json


def main():
    plot_images()

    # Load the data from the target JSON file
    data_5 = load_json('zp5_list.json')
    zp_list_5 = data_5['zp_list']
    color_list_5 = data_5['color_list']
    fluxes_list_5 = data_5['flux_list']
    mags_list_5 = data_5['tmag_list']

    data_4 = load_json('zp4_list.json')
    zp_list_4 = data_4['zp_list']
    color_list_4 = data_4['color_list']
    fluxes_list_4 = data_4['flux_list']
    mags_list_4 = data_4['tmag_list']

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
    print(f'Avg zp5: {np.mean(zp_list_5)}')
    print(f'Avg zp4: {np.mean(zp_list_4)}')

    # choose the value with the highest frequency
    from collections import Counter
    print(f'Mode zp5: {Counter(zp_list_5).most_common(1)}')
    print(f'Mode zp4: {Counter(zp_list_4).most_common(1)}')
    plt.show()

    # color-coded with color index
    # exclude stars with Tmag < 9.5
    tmag_filtered_list_5 = []
    flux_filtered_list_5 = []
    color_filtered_list_5 = []
    print(f'Number of stars: {len(mags_list_5)}')
    for i in range(len(mags_list_5)):
        if mags_list_5[i] > 9.5 and color_list_5[i] < 1.3:
            tmag_filtered_list_5.append(mags_list_5[i])
            flux_filtered_list_5.append(2.5 * np.log10(fluxes_list_5[i] / 10))
            color_filtered_list_5.append(color_list_5[i])

    # Perform linear fit
    print(f'Number of stars: {len(tmag_filtered_list_5)}')
    slope, intercept = np.polyfit(tmag_filtered_list_5, flux_filtered_list_5, 1)
    print(f'The zeropoint CMOS is {intercept}')

    # Plot the data
    plt.scatter(
        tmag_filtered_list_5,
        flux_filtered_list_5,
        c=color_filtered_list_5,
        cmap='coolwarm',
        vmin=0.5,
        vmax=1.5,
    )
    plt.colorbar(label=r'$\mathdefault{G_{BP}-G_{RP}}$')

    # Add fitted line
    fit_x = np.array(tmag_filtered_list_5)
    fit_y = intercept + slope * fit_x  # Linear fit on log scale
    plt.plot(fit_x, fit_y, color='black', label=f'Fit: y = {slope:.2f}x + {intercept:.2f}')

    # Adjust axes
    plt.yscale('linear')  # Linear for easier comparison since log10 is precomputed
    plt.ylabel(r'$\mathdefault{2.5 \log_{10}(Flux/10)}$')
    plt.xlabel('Tmag')
    plt.gca().invert_xaxis()
    plt.xlim(9, 16.5)

    # Add legend
    plt.legend()
    plt.show()

    # color-coded with color index
    # exclude stars with Tmag < 9.5
    tmag_filtered_list_4 = []
    flux_filtered_list_4 = []
    color_filtered_list_4 = []

    for i in range(len(mags_list_4)):
        if mags_list_4[i] > 9.5 and color_list_4[i] < 1.3:
            tmag_filtered_list_4.append(mags_list_4[i])
            flux_filtered_list_4.append(2.5 * np.log10(fluxes_list_4[i] / 10))
            color_filtered_list_4.append(color_list_4[i])

    # Perform linear fit
    slope, intercept = np.polyfit(tmag_filtered_list_4, flux_filtered_list_4, 1)
    print(f'The zeropoint CCD is {intercept}')

    # Plot the data
    plt.scatter(
        tmag_filtered_list_4,
        flux_filtered_list_4,
        c=color_filtered_list_4,
        cmap='coolwarm',
        vmin=0.5,
        vmax=1.5,
    )
    plt.colorbar(label=r'$\mathdefault{G_{BP}-G_{RP}}$')

    # Add fitted line
    fit_x = np.array(tmag_filtered_list_4)
    fit_y = intercept + slope * fit_x  # Linear fit on log scale
    plt.plot(fit_x, fit_y, color='black', label=f'Fit: y = {slope:.2f}x + {intercept:.2f}')

    # Adjust axes
    plt.yscale('linear')  # Linear for easier comparison since log10 is precomputed
    plt.ylabel(r'$\mathdefault{2.5 \log_{10}(Flux/10)}$')
    plt.xlabel('Tmag')
    plt.gca().invert_xaxis()
    plt.xlim(9, 16.5)

    # Add legend
    plt.legend()
    plt.show()


main()
