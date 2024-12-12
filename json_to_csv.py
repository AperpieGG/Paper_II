import json
import pandas as pd

# File path
path = '/Users/u5500483/Documents/GitHub/Paper_II/data_transit/target_light_curve_144065872_CMOS_0714.json'


def load_json(file_path):
    """
    Load JSON file and return data.
    """
    with open(file_path, 'r') as file:
        return json.load(file)


def convert_json_to_csv(file_path):
    data = load_json(file_path)

    # Determine the reference time from the first data point
    reference_time = data['Time_BJD'][0]
    print(f"Reference time (first BJD value): {reference_time}")

    # Convert BJD to days
    time_days = [bjd - 2460506 for bjd in data['Time_BJD']]
    relative_flux = data['Relative_Flux']
    relative_flux_err = data['Relative_Flux_err']

    # Create DataFrame with the correct column names
    df = pd.DataFrame({
        '# time': time_days,
        'flux': relative_flux,
        'flux_err': relative_flux_err
    })

    # Save as CSV with the correct file name
    path_save = '/Users/u5500483/Documents/GitHub/Paper_II/allesfitter/'
    df.to_csv(path_save + 'NGTS.csv', index=False)
    print("CSV file saved as NGTS.csv")


# Run the conversion
convert_json_to_csv(path)
