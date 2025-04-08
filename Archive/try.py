import json

# Load the JSON file
file_path = "../target_light_curve_9725627_CCD_0828.json"  # Replace with your JSON file path

with open(file_path, "r") as file:
    data = json.load(file)

# Extract the parameters
time_bjd = data["Time_BJD"]
relative_flux = data["Relative_Flux"]
relative_flux_err = data["Relative_Flux_err"]

# Threshold value
threshold = 2460551.8224

# Filter the data
filtered_time_bjd = []
filtered_relative_flux = []
filtered_relative_flux_err = []

for t, f, e in zip(time_bjd, relative_flux, relative_flux_err):
    if t <= threshold:
        filtered_time_bjd.append(t)
        filtered_relative_flux.append(f)
        filtered_relative_flux_err.append(e)

# Update the JSON structure
filtered_data = {
    "Time_BJD": filtered_time_bjd,
    "Relative_Flux": filtered_relative_flux,
    "Relative_Flux_err": filtered_relative_flux_err
}

# Save the filtered data back to a new JSON file
output_file_path = "../target_light_curve_9725627_CCD_0828_c.json"  # Replace with your desired output file path

with open(output_file_path, "w") as file:
    json.dump(filtered_data, file, indent=4)

print(f"Filtered data saved to {output_file_path}.")