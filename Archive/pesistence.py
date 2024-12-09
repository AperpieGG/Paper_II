import argparse
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from datetime import datetime, timedelta
from set_temperature import Temperature
from pyAndorSDK import AndorSDK3, CameraException, ErrorCodes

def initialize_camera():
    cam = AndorSDK3().GetCamera(0)
    return cam

def save_images(list_images):
    for i, image in enumerate(list_images):
        hdu = fits.PrimaryHDU(image)
        hdu.writeto(f'./image_{i}.fits', overwrite=True)

def acquire_bright_images(cam):
    cam.ShutterMode = 'Closed'
    num_frames = 1
    exposure = 60
    cam.FrameCount = num_frames
    cam.ExposureTime = exposure
    bright_images = []

    Imageseries_bright = cam.acquire_series(timeout=50000)
    brightimage = np.zeros([cam.AOIHeight, cam.AOIWidth])
    for i in range(num_frames):
        brightimage += Imageseries_bright[i].image.astype(np.float32)
        bright_images.append(brightimage[1:50, 1:50])
    return bright_images

def acquire_dark_images(cam, num_frames):
    cam.ShutterMode = 'Closed'
    exposure = 1
    cam.FrameCount = num_frames
    cam.ExposureTime = exposure
    dark_images = []

    Imageseries_dark = cam.acquire_series(timeout=50000)
    darkimage = np.zeros([cam.AOIHeight, cam.AOIWidth])
    for i in range(num_frames):
        darkimage += Imageseries_dark[i].image.astype(np.float32)
        dark_images.append(darkimage[1:50, 1:50])
    return dark_images

from mpl_toolkits.axes_grid1.inset_locator import inset_axes


def main():
    parser = argparse.ArgumentParser(description='Acquire bright and dark images')
    parser.add_argument('temperature', type=float, help='Temperature to set the camera')
    args = parser.parse_args()

    cam = initialize_camera()

    target_temp = args.temperature
    Temperature(cam).SetTemperature(target_temp)

    # Acquire images
    dark_images_before = acquire_dark_images(cam, num_frames=5)
    bright_images = acquire_bright_images(cam)
    dark_images_after = acquire_dark_images(cam, num_frames=5)

    save_images(dark_images_before)
    save_images(bright_images)
    save_images(dark_images_after)

    # Prepare to plot mean signal vs. relative time
    mean_signals = []
    relative_time = []

    # Time settings
    dark_exposure_time = 1.0  # Each dark image has a 1-second exposure
    bright_exposure_time = 20.0  # Bright image has a 20-second exposure
    time_between_images = 0.04  # Time between images is 0.04 seconds
    total_time_between_dark_images = dark_exposure_time + time_between_images
    total_time_between_bright_images = bright_exposure_time + time_between_images

    # Before bright image
    for i in range(len(dark_images_before)):
        mean_signal = np.mean(dark_images_before[i])
        mean_signals.append(mean_signal)
        relative_time.append(-total_time_between_dark_images * (len(dark_images_before) - i))

    # Bright image at time = 0
    for i in range(len(bright_images)):
        mean_signal = np.mean(bright_images[i])
        mean_signals.append(mean_signal)
        relative_time.append(0)  # Bright image taken at t = 0 seconds

    # After bright image
    for i in range(len(dark_images_after)):
        mean_signal = np.mean(dark_images_after[i])
        mean_signals.append(mean_signal)
        # Account for the time of the bright exposure
        relative_time.append(total_time_between_dark_images * (i + 1) + bright_exposure_time)

    # Plot mean signals as a function of relative time
    fig, ax = plt.subplots(figsize=(12, 8))
    ax.plot(relative_time[:len(dark_images_before)], mean_signals[:len(dark_images_before)], 'o-', color='blue', label='Dark Images')  # Plot dark images in blue
    ax.plot(relative_time[len(dark_images_before):len(dark_images_before) + len(bright_images)],
            mean_signals[len(dark_images_before):len(dark_images_before) + len(bright_images)], 'o-', color='red', label='Bright Image')  # Plot bright image in red
    ax.plot(relative_time[len(dark_images_before) + len(bright_images):],
            mean_signals[len(dark_images_before) + len(bright_images):], 'o-', color='blue')  # Continue with dark images
    ax.set_xlabel('Relative Time (s)')
    ax.set_ylabel('Mean Signal')
    ax.set_title('Mean Signal vs. Relative Time')
    ax.grid(True)
    ax.axvline(0, color='red', linestyle='--', label='Bright Image')
    ax.legend()

    # Adding inset zoom images before and after the bright image
    # Plot zoomed imshow for the last dark image before the bright image
    axins_before = inset_axes(ax, width="20%", height="20%", loc='upper left', borderpad=3)
    im_before = axins_before.imshow(dark_images_before[-1], cmap='hot', origin='lower')
    axins_before.set_xticks([])
    axins_before.set_yticks([])
    axins_before.set_title("Before Bright")

    # Add colorbar to the inset
    cbar_before = plt.colorbar(im_before, ax=axins_before)
    cbar_before.set_label('Intensity')

    # Plot zoomed imshow for the first dark image after the bright image
    axins_after = inset_axes(ax, width="20%", height="20%", loc='upper right', borderpad=3)
    im_after = axins_after.imshow(dark_images_after[0], cmap='hot', origin='lower')
    axins_after.set_xticks([])
    axins_after.set_yticks([])
    axins_after.set_title("After Bright")

    # Add colorbar to the inset
    cbar_after = plt.colorbar(im_after, ax=axins_after)
    cbar_after.set_label('Intensity')

    # Add annotation arrows pointing to the zoomed-in images
    # For the image before the bright image
    ax.annotate('Image Before Bright',
                xy=(relative_time[len(dark_images_before) - 1], mean_signals[len(dark_images_before) - 1]),
                xytext=(-150, -100),
                textcoords='offset points',
                arrowprops=dict(facecolor='black', shrink=0.05),
                fontsize=10, color='blue')

    # For the image after the bright image
    ax.annotate('Image After Bright',
                xy=(relative_time[len(dark_images_before) + len(bright_images)], mean_signals[len(dark_images_before) + len(bright_images)]),
                xytext=(150, -100),
                textcoords='offset points',
                arrowprops=dict(facecolor='black', shrink=0.05),
                fontsize=10, color='blue')

    plt.show()

if __name__ == "__main__":
    main()