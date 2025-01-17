import math


def calculate_ratio(Rb_to_Rstar, a_to_Rstar):
    # Calculate the numerator and denominator
    numerator = 1 + Rb_to_Rstar  # R_star + R_b = R_star * (1 + Rb/R_star)
    denominator = a_to_Rstar  # a = a_to_Rstar * R_star

    # Calculate the ratio
    ratio = numerator / denominator
    return ratio


# Given values
Rb_to_Rstar = 0.1091  # Rb / R_star
a_to_Rstar = 6.59  # a / R_star

# Calculate cosine of 89.06 degrees
angle_in_degrees = 88.0
cosine_value = math.cos(math.radians(angle_in_degrees))

# Print the result

# Calculate and print the ratio
result = calculate_ratio(Rb_to_Rstar, a_to_Rstar)
print(f"The ratio (R_star + R_b) / a_b is approximately: {result:.3f}")
print(f"The cos(i) of {angle_in_degrees} degrees is approximately: {cosine_value:.6f}")
