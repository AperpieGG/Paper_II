import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from plot_images import plot_images

plot_images()

# Load CSV files
path = '/allesfit/'
file1 = path + 'comparison_from_literature.csv'
file2 = path + 'comparison_from_measured.csv'

df1 = pd.read_csv(file1)
df2 = pd.read_csv(file2)

# Define columns to compare
columns = [
    ('Rp/Rs', 'Rp/Rs_err'),
    ('Depth (%)', 'Depth_err'),
    ('T0 (BJD)', 'T0_err'),
    ('Duration (days)', 'Duration_err')
]

# Set up subplot
fig, axs = plt.subplots(2, 2, figsize=(10, 8))
axs = axs.flatten()

# Style and colour mapping
styles = {'CMOS': ('o', 'blue'), 'CCD': ('s', 'red')}

# Plotting
for i, (col, err_col) in enumerate(columns):
    ax = axs[i]
    all_vals = []  # for determining min and max across all values
    for source in ['CMOS', 'CCD']:
        marker, colour = styles[source]
        subset = df2[df2['Source'] == source]
        for _, row in subset.iterrows():
            target = row['Target']
            match = df1[df1['Target'] == target]
            if not match.empty:
                x = match.iloc[0][col]
                xerr = match.iloc[0][err_col]
                y = row[col]
                yerr = row[err_col]
                ax.errorbar(
                    x, y,
                    xerr=xerr, yerr=yerr,
                    fmt=marker,
                    color=colour,
                    markersize=5
                )
                all_vals.extend([x, y])

    # Set axis labels
    ax.set_xlabel(f'{col} (Literature)')
    ax.set_ylabel(f'{col} (Measured)')
    ax.grid(True)

    # Plot y = x reference line
    min_val = min(all_vals)
    max_val = max(all_vals)
    if not columns[i] == ('T0 (BJD)', 'T0_err'):
        ax.plot([min_val, max_val], [min_val, max_val], 'k--', linewidth=1)

# Create custom legend handles
legend_elements = [
    Line2D([0], [0], marker='o', color='blue', label='CMOS', linestyle='None'),
    Line2D([0], [0], marker='s', color='red', label='CCD', linestyle='None')
]

# Add one shared legend at the top centre
fig.legend(handles=legend_elements, loc='upper center', ncol=2, frameon=False, fontsize='medium')

plt.tight_layout(rect=[0, 0, 1, 0.95])  # Leave space at top for legend
plt.show()