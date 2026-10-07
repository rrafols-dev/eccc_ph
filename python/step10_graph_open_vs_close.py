import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from root import ROOT

# ===========================================================================
# Set paths
# ===========================================================================
OUTPUTS = ROOT + 'outputs_csv/mig_full/'
OUTPUTS2 = ROOT + 'tabFig/'

# Load baseline data - closed and open economy
baseline_closed = pd.read_csv(OUTPUTS + 'resultspy_new_2010.csv')
baseline_open = pd.read_csv(OUTPUTS + 'resultspy_new_2010_ROW.csv')

# Extract fundamentals
barAs_closed = baseline_closed['barAs'].values
barAu_closed = baseline_closed['barAu'].values
barTs_closed = baseline_closed['barTs'].values
barTu_closed = baseline_closed['barTu'].values

barAs_open = baseline_open['barAs'].values
barAu_open = baseline_open['barAu'].values
barTs_open = baseline_open['barTs'].values
barTu_open = baseline_open['barTu'].values

# ===========================================================================
# Create Figure: Closed vs Open Economy Fundamentals
# ===========================================================================
def plot_scatter(ax, x, y, title, xlabel='Closed Economy', ylabel='Open Economy'):
    # Filter out NaN values
    valid_mask = ~(np.isnan(x) | np.isnan(y))
    x_valid = x[valid_mask]
    y_valid = y[valid_mask]
    
    ax.scatter(x_valid, y_valid, s=20, alpha=0.5, color='gray')
    
    # Add 45-degree line
    min_val = min(x_valid.min(), y_valid.min())
    max_val = max(x_valid.max(), y_valid.max())
    ax.plot([min_val, max_val], [min_val, max_val], 'b-', linewidth=1.5)
    
    # Calculate correlation on valid values only
    if len(x_valid) > 0:
        corr = np.corrcoef(x_valid, y_valid)[0, 1]
        ax.text(0.05, 0.95, f'Correlation: {corr:.4f}', 
                transform=ax.transAxes, fontsize=9, verticalalignment='top',
                bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.set_title(title, fontsize=11)
    ax.grid(True, alpha=0.3)
    
    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# Create 2x2 figure
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Top left: Amenity, Skilled
plot_scatter(axes[0, 0], barAs_closed, barAs_open, r'$\bar{A}_H$ (Amenity, Skilled)')

# Top right: Amenity, Unskilled
plot_scatter(axes[0, 1], barAu_closed, barAu_open, r'$\bar{A}_U$ (Amenity, Unskilled)')

# Bottom left: Productivity, Skilled
plot_scatter(axes[1, 0], barTs_closed, barTs_open, r'$\bar{T}_H$ (Productivity, Skilled)')

# Bottom right: Productivity, Unskilled
plot_scatter(axes[1, 1], barTu_closed, barTu_open, r'$\bar{T}_U$ (Productivity, Unskilled)')

plt.tight_layout()
plt.savefig(OUTPUTS2 + 'fundamentals_closed_vs_open.png', dpi=300, bbox_inches='tight')
plt.show()

print("Figure saved successfully!")

# Print correlations with NaN filtering
valid_mask_As = ~(np.isnan(barAs_closed) | np.isnan(barAs_open))
valid_mask_Au = ~(np.isnan(barAu_closed) | np.isnan(barAu_open))
valid_mask_Ts = ~(np.isnan(barTs_closed) | np.isnan(barTs_open))
valid_mask_Tu = ~(np.isnan(barTu_closed) | np.isnan(barTu_open))

print("\nCorrelations:")
print(f"  barAs: {np.corrcoef(barAs_closed[valid_mask_As], barAs_open[valid_mask_As])[0, 1]:.4f}")
print(f"  barAu: {np.corrcoef(barAu_closed[valid_mask_Au], barAu_open[valid_mask_Au])[0, 1]:.4f}")
print(f"  barTs: {np.corrcoef(barTs_closed[valid_mask_Ts], barTs_open[valid_mask_Ts])[0, 1]:.4f}")
print(f"  barTu: {np.corrcoef(barTu_closed[valid_mask_Tu], barTu_open[valid_mask_Tu])[0, 1]:.4f}")