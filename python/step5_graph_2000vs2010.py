import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from root import ROOT

# Paths
OUTPUTS = ROOT + 'outputs_csv/mig_full/'
OUTPUTS2 = ROOT + 'tabFig/'

# Load inverted results from both years
results_2000 = pd.read_csv(OUTPUTS + 'resultspy_new_2000.csv')
results_2010 = pd.read_csv(OUTPUTS + 'resultspy_new_2010.csv')

# Extract fundamentals
barAs_2000 = results_2000['barAs'].values
barAu_2000 = results_2000['barAu'].values
barTs_2000 = results_2000['barTs'].values
barTu_2000 = results_2000['barTu'].values

barAs_2010 = results_2010['barAs'].values
barAu_2010 = results_2010['barAu'].values
barTs_2010 = results_2010['barTs'].values
barTu_2010 = results_2010['barTu'].values

# ===========================================================================
# Create Figure A3: Location Fundamentals From Model Inversion, 2000 vs 2010
# ===========================================================================
def plot_scatter(ax, x, y, title, xlabel='2010', ylabel='2000'):
    ax.scatter(x, y, s=15, alpha=0.4, color='gray')
    
    # Calculate correlation
    corr = np.corrcoef(x, y)[0, 1]
    
    # Add 45-degree line
    min_val = min(x.min(), y.min())
    max_val = max(x.max(), y.max())
    ax.plot([min_val, max_val], [min_val, max_val],  linewidth=2)
    
    ax.set_xlabel(xlabel, fontsize=11)
    ax.set_ylabel(ylabel, fontsize=11)
    ax.set_title(title, fontsize=12)
    
    # Add correlation text
    ax.text(0.05, 0.95, f'Corr: {corr:.2f}***', 
            transform=ax.transAxes, fontsize=10, 
            verticalalignment='top',
            bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))
    
    # Remove top and right spines
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# Create 2x2 figure
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Top left: Amenity, Low-skilled
plot_scatter(axes[0, 0], barAu_2010, barAu_2000, 'Amenity, Low-skilled')

# Top right: Amenity, Skilled
plot_scatter(axes[0, 1], barAs_2010, barAs_2000, 'Amenity, Skilled')

# Bottom left: Productivity, Low-skilled
plot_scatter(axes[1, 0], barTu_2010, barTu_2000, 'Productivity, Low-skilled')

# Bottom right: Productivity, Skilled
plot_scatter(axes[1, 1], barTs_2010, barTs_2000, 'Productivity, Skilled')

plt.tight_layout()
plt.savefig(OUTPUTS2 + 'fundamentals_2000_2010.png', dpi=300, bbox_inches='tight')
plt.show()
