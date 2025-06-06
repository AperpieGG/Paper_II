from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
from astropy.visualization import ZScaleInterval
from matplotlib.patches import Circle, FancyArrow
from matplotlib.colors import LogNorm
from plot_images import plot_images
from matplotlib import ticker


plot_images()


# Define functions
def read_image(filename, subtract_value, multiply_value, micron):
    """
    Read and preprocess an image from a FITS file.
    """
    with fits.open(filename) as hdul:
        image = hdul[0].data.astype(float)
        print(f'Mean pixel value: {np.mean(image)}')
        preprocessed_image = ((image - subtract_value) * multiply_value) / 10 / (micron**2)
        print(f'Preprocessed mean pixel value: {np.mean(preprocessed_image)}')
    return preprocessed_image


def add_compass(ax, origin_x, origin_y, length=4):
    """
    Add North and East arrows to the plot.
    """
    # Add North arrow
    ax.arrow(origin_x, origin_y, 0, length, head_width=1, head_length=1,
             fc='white', ec='white', lw=1)
    ax.text(origin_x, origin_y + length + 2, "N", color='white', fontsize=10, ha='center')

    # Add East arrow
    ax.arrow(origin_x, origin_y, -length, 0, head_width=1, head_length=1,
             fc='white', ec='white', lw=1)
    ax.text(origin_x - length - 3, origin_y, "E", color='white', fontsize=10, va='center')


def find_position():
    """
    Return x-y positions for both CMOS and CCD images.
    """
    x_cmos = 1025.0089593013547
    y_cmos = 895.7486945591139
    x_ccd = 1033.2574339079276
    y_ccd = 1131.0896942990596
    return x_cmos, y_cmos, x_ccd, y_ccd


def plot_image(ax, image_data, x, y, title, aperture, radius=25.5, origin_setting='lower', show_compass=False, pixel_scale=None):
    """
    Plot an image with specified parameters and overlay aperture, annulus, and scale bar.
    """
    # Define cropping limits
    x_min = max(int(x - radius), 0)
    x_max = min(int(x + radius), image_data.shape[1])
    y_min = max(int(y - radius), 0)
    y_max = min(int(y + radius), image_data.shape[0])

    # Crop the image
    cropped_image_data = image_data[y_min:y_max, x_min:x_max]

    # Plot the image
    extent = [x - radius, x + radius, y - radius, y + radius]
    im = ax.imshow(cropped_image_data, cmap='viridis', origin=origin_setting, extent=extent, norm=norm)

    ax.set_title(title)
    ax.set_xlabel('X Pixel')
    if title != "CMOS Frame":
        ax.set_ylabel('Y Pixel')

    # Plot aperture
    circle = Circle((x, y), radius=aperture, edgecolor='lime', facecolor='none', lw=1)
    ax.add_patch(circle)

    annulus = Circle((x, y), radius=15, edgecolor='lime', facecolor='none', lw=1, linestyle='dashed')
    dannulus = Circle((x, y), radius=20, edgecolor='lime', facecolor='none', lw=1, linestyle='dashed')
    ax.add_patch(annulus)
    ax.add_patch(dannulus)

    # Add compass
    if show_compass:
        add_compass(ax, x + radius - 2, y - radius + 2)

    # Add scale bar (1 arcminute)
    if pixel_scale:
        arcmin_pixels = 60 / pixel_scale  # 60 arcsec in 1 arcminute
        bar_length = arcmin_pixels
        bar_x = x - radius + 2
        bar_y = y - radius + 2
        ax.plot([bar_x, bar_x + bar_length], [bar_y, bar_y], color='white', lw=3)
        ax.text(bar_x + bar_length / 2, bar_y + 0.5, r"1$'$", color='white', ha='center', va='bottom', fontsize=14)

        return im


# Main Execution
if __name__ == "__main__":
    # File paths
    path = '/Users/u5500483/Downloads/'
    ccd_file = path + "IMAGE81120240801232349.fits"
    cmos_file = path + "NG1858-4651_TIC-269217040_811_S43-20240801232406275.fits"

    # Preprocess images
    ccd_image = read_image(ccd_file, subtract_value=1586, multiply_value=2, micron=5.01)
    cmos_image = read_image(cmos_file, subtract_value=150, multiply_value=1.13, micron=4.01)
    print(f'Averge pixel value of CCD image: {np.median(ccd_image)}')
    print(f'Averge pixel value of CMOS image: {np.median(cmos_image)}')
    # Positions for the stars
    x_cmos, y_cmos, x_ccd, y_ccd = find_position()
    norm = LogNorm(vmin=np.min(cmos_image), vmax=100)

    # Plot images side by side
    fig, axs = plt.subplots(1, 2, figsize=(8, 3.6))

    # Plot CCD image with origin='upper'
    im_ccd = plot_image(axs[0], ccd_image, x_ccd, y_ccd, aperture=4,  title="CCD Frame",
                        origin_setting='upper', show_compass=True, pixel_scale=5.01)

    # Plot CMOS image
    im_cmos = plot_image(axs[1], cmos_image, x_cmos, y_cmos, aperture=5, title="CMOS Frame",
                         origin_setting='lower', show_compass=True, pixel_scale=4.01)
    # Create colorbar matching the size of the images
    cbar = fig.colorbar(im_cmos, ax=axs.ravel().tolist(), orientation='vertical', fraction=0.021, pad=0.04)
    cbar.set_label("Pixel value ($\mathdefault{e^{-}s^{-1}\,arcsec^{-2}}$)")

    # Customizing the colorbar labels in ke/s
    cbar.set_label("Pixel value ($\mathdefault{e^{-}s^{-1}\,arcsec^{-2}}$)")
    # cbar.set_ticks([1, 2, 4, 5, 10, 20, 50, 100])
    # cbar.ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda x, _: f"{int(x)}"))
    plt.savefig(path + 'CMOS_CCD_Frames.pdf', dpi=300)
    plt.show()
