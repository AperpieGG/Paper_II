import csv
import os
from matplotlib import pyplot as plt
from plot_images import plot_images

plot_images()

TARGETS = [
    'HATS-33_0715_BIG', 'KELT-10_0801', 'KELT-10_0826', 'TOI-905_0713', 'TOI-2109_0717', 'TOI-4463_0827',
    'WASP-4_0802', 'WASP-30_0828', 'WASP-95_0714', 'WASP-95_0807', 'WASP-97_0804', 'WASP-97_0806'
]

PARAMS_TO_EXTRACT = ['b_rr', 'b_rsuma', 'b_cosi', 'b_epoch', 'b_period']
PARAMS_TO_PLOT = ['b_rr', 'b_rsuma', 'b_cosi', 'b_period']

param_titles = {
    'b_rr': r'$R_b / R_\star$',
    'b_rsuma': r'$(R_\star + R_b) / a_b$',
    'b_cosi': r'$\cos{i_b}$',
    'b_period': r'$P_b$ [days]'
}

cam_colours = {'CMOS': 'blue', 'CCD': 'green'}
cam_markers = {'CMOS': 'o', 'CCD': 's'}


def parse_bounds(bounds_str):
    parts = bounds_str.strip().split()
    if len(parts) == 3 and parts[0].lower() == 'normal':
        return float(parts[1]), float(parts[2])
    else:
        return None, None


def load_csv_params(csv_path):
    results = {}
    with open(csv_path, 'r') as file:
        reader = csv.reader(file)
        for row in reader:
            if not row or row[0].startswith('#'):
                continue
            name, value, fit, bounds, *_ = (row + [''] * 6)[:6]
            name = name.strip()
            if name in PARAMS_TO_EXTRACT:
                try:
                    value = float(value.strip())
                    mean, std = parse_bounds(bounds)
                    results[name] = (value, std)
                except:
                    continue
    return results


def load_ns_table_csv(ns_path):
    result = {}
    with open(ns_path, 'r') as f:
        reader = csv.reader(f)
        for row in reader:
            if not row or row[0].startswith('#'):
                continue
            name, median, lower, upper, *_ = (row + [''] * 5)[:5]
            name = name.strip()
            if name not in PARAMS_TO_EXTRACT:
                continue
            if median.strip() in ['(fixed)', ''] or lower.strip() in ['(fixed)', ''] or upper.strip() in ['(fixed)',
                                                                                                          '']:
                continue
            try:
                median_val = float(median.strip())
                lower_err = float(lower.strip())
                upper_err = float(upper.strip())
                std_estimate = (lower_err + upper_err) / 2
                result[name] = (median_val, std_estimate)
            except ValueError:
                continue
    return result


# Create plot once
fig, axes = plt.subplots(2, 2, figsize=(10, 8), constrained_layout=True)

# Now loop over each camera and plot onto the same subplots
for cam in ['CMOS', 'CCD']:
    print(f"\n--- Processing camera: {cam} ---")

    for TARGET in TARGETS:
        path = f'/Users/u5500483/Documents/GitHub/Paper_II/allesfit/Final_Version_Results/allesfit_{cam}_{TARGET}/'
        csv_path = os.path.join(path, 'params.csv')
        ns_table_path = os.path.join(path, 'results/ns_table.csv')

        if not os.path.exists(csv_path) or not os.path.exists(ns_table_path):
            print(f"{TARGET} ({cam}): Missing file(s). Skipping.")
            continue

        csv_params = load_csv_params(csv_path)
        ns_params = load_ns_table_csv(ns_table_path)

        for idx, param in enumerate(PARAMS_TO_PLOT):
            if param not in csv_params or param not in ns_params:
                continue

            x_val, x_std = csv_params[param]
            y_val, y_std = ns_params[param]
            if x_val is None or y_val is None:
                continue

            row, col = divmod(idx, 2)
            ax = axes[row][col]
            ax.errorbar(
                x_val, y_val,
                xerr=x_std if x_std else 0,
                yerr=y_std if y_std else 0,
                fmt=cam_markers[cam],
                color=cam_colours[cam],
                label=cam if TARGET == TARGETS[0] else "",  # Only label once per cam
                capsize=4,
                alpha=0.8
            )
            ax.set_xlabel(f"Literature: {param_titles.get(param, param)}")
            ax.set_ylabel(f"Measured: {param_titles.get(param, param)}")
            ax.grid(True)

# Add red y = x line and legend to all axes
for idx, param in enumerate(PARAMS_TO_PLOT):
    row, col = divmod(idx, 2)
    ax = axes[row][col]
    ax.autoscale()
    xlims = ax.get_xlim()
    ylims = ax.get_ylim()
    lims = [min(xlims[0], ylims[0]), max(xlims[1], ylims[1])]
    ax.plot(lims, lims, 'r--')
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.legend()

plt.show()
