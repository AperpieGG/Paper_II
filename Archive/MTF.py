import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter
from scipy.fftpack import fft, fftfreq
import cv2


def create_slanted_edge_image(size, angle, edge_position):
    """
    Create a slanted edge image.
    :param size: tuple of (height, width)
    :param angle: angle of the slanted edge in degrees
    :param edge_position: position of the edge (0 to 1) relative to the width
    :return: slanted edge image
    """
    height, width = size
    edge_image = np.zeros((height, width))
    x = np.arange(width)
    y = np.tan(np.radians(angle)) * (-x) + edge_position * width
    for i in range(height):
        edge_image[i, :] = x < y[i]
    return edge_image


# Parameters
image_size = (512, 512)
edge_angle = 5  # degrees
edge_position = 0.5  # center

# Create slanted edge image
edge_image = create_slanted_edge_image(image_size, edge_angle, edge_position)


def add_noise_and_blur(image, sigma_blur=1.0, sigma_noise=0.01):
    """
    Add Gaussian blur and noise to an image.
    :param image: input image
    :param sigma_blur: standard deviation for Gaussian blur
    :param sigma_noise: standard deviation for Gaussian noise
    :return: noisy and blurred image
    """
    blurred_image = gaussian_filter(image, sigma=sigma_blur)
    noisy_image = blurred_image + sigma_noise * np.random.randn(*image.shape)
    return noisy_image


# Add noise and blur
simulated_image = add_noise_and_blur(edge_image)

# Display the simulated image
plt.imshow(simulated_image, cmap='gray')
plt.title('Simulated Image')
plt.show()


def compute_mtf(image, edge_angle):
    """
    Compute the MTF from a slanted edge image.
    :param image: slanted edge image
    :param edge_angle: angle of the slanted edge in degrees
    :return: spatial frequencies and MTF
    """
    # Rotate the image to make the edge vertical
    rotation_matrix = cv2.getRotationMatrix2D((image.shape[1] / 2, image.shape[0] / 2), -edge_angle, 1)
    rotated_image = cv2.warpAffine(image, rotation_matrix, (image.shape[1], image.shape[0]))

    # Crop the central region
    crop_size = (100, image.shape[1])
    center = (rotated_image.shape[0] // 2, rotated_image.shape[1] // 2)
    cropped_image = rotated_image[center[0] - crop_size[0] // 2:center[0] + crop_size[0] // 2, :]

    # Compute the edge spread function (ESF)
    esf = np.mean(cropped_image, axis=0)

    # Compute the line spread function (LSF)
    lsf = np.diff(esf)

    # Compute the MTF
    mtf = np.abs(fft(lsf))
    mtf = mtf / np.max(mtf)

    # Compute the spatial frequencies
    freqs = fftfreq(len(lsf), d=1)

    return freqs[:len(freqs) // 2], mtf[:len(mtf) // 2]


# Compute MTF
freqs, mtf = compute_mtf(simulated_image, edge_angle)

# Plot MTF
plt.plot(freqs, mtf, 'b-')
plt.xlabel('Spatial Frequency')
plt.ylabel('MTF')
plt.title('MTF Plot')
plt.grid(True)
plt.show()
