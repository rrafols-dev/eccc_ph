import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from cf_shocks import compute_fundamental_shocks
from root import ROOT

# ===========================================================================
# Set universal parameters
# ===========================================================================
PARAMS = {
    'sigma': 4,
    'theta': 1.5,
    'kappa': 1.5,
    'eta': -0.05,
    'mus': 0.02,
    'muu': 0.02,
    'beta_hat': -1.199,
    'betas_hat': -1.171,
    'betau_hat': -1.315,
    'islands': 0.274,
    'islandu': 0.145,
    'provs': 1.37,
    'provu': 1.146,
    'homes': 2.911,
    'homeu': 2.579,
    'longs': 0.164,
    'longu': 0.291,
    'lats': -0.02,
    'latu': -0.098,
    'N': 1600,
    'shock_amen': 1,  #counterfactual channels allowed to react to climate
    'shock_prod': 1,
    'simYear': 2100,
    'tempAs':  -3.48,  #climate elasticities (w/ regionFE)
    'tempAu': -0.3, 
    'tempTs': -4.36,
    'tempTu': -9.30,
    'migsimple': 0 
}

#Universal Inputs
CSVFILES = ROOT + 'inputs_csv/'
OUTPUTS = ROOT + 'outputs_csv/mig_full/'
OUTPUTS2 = ROOT + 'tabFig/'
BASELINE_FILE = 'resultspy_new_2010.csv'
climate_data = np.array(pd.read_csv(CSVFILES + 'climate.csv'))
climate_bounds = pd.read_csv(CSVFILES + 'climate_bounds.csv')

data_slr = np.array(pd.read_csv(CSVFILES + 'slr_area.csv'))
area = data_slr[:, [1]]

print("Loading universal distance and geography matrices...")
data_dist = np.array(pd.read_csv(CSVFILES + 'distance.csv'))
data_dist = np.delete(data_dist, 1597, axis=0)
data_dist = np.delete(data_dist, 1598, axis=1)
data_dist = data_dist[:, 1:]

diff_long = np.array(pd.read_csv(CSVFILES + 'diff_long.csv'))
diff_long = np.delete(diff_long, 1597, axis=0)
diff_long = np.delete(diff_long, 1598, axis=1)
diff_long = diff_long[:, 1:]

diff_lat = np.array(pd.read_csv(CSVFILES + 'diff_lat.csv'))
diff_lat = np.delete(diff_lat, 1597, axis=0)
diff_lat = np.delete(diff_lat, 1598, axis=1)
diff_lat = diff_lat[:, 1:]

sameprov = np.array(pd.read_csv(CSVFILES + 'sameprov.csv'))
sameprov = np.delete(sameprov, 1597, axis=0)
sameprov = np.delete(sameprov, 1598, axis=1)
sameprov = sameprov[:, 1:]

sameisland = np.array(pd.read_csv(CSVFILES + 'sameisland.csv'))
sameisland = np.delete(sameisland, 1597, axis=0)
sameisland = np.delete(sameisland, 1598, axis=1)
sameisland = sameisland[:, 1:]

id_2 = data_slr[:,[0]] 
N = PARAMS['N']
I = np.eye(N)



# ===========================================================================
# Calculate shocks to exogenous objects
# ===========================================================================
print("Computing hats...")
# Compute fundamental shocks
barAs_hat, barAu_hat, barTs_hat, barTu_hat = compute_fundamental_shocks(
    PARAMS['shock_amen'], PARAMS['shock_prod'], PARAMS['simYear'],
    PARAMS['tempAs'], PARAMS['tempAu'], PARAMS['tempTs'], PARAMS['tempTu'],
    climate_data, OUTPUTS, BASELINE_FILE
)


# Load baseline 2010 data
baseline_2010 = pd.read_csv(OUTPUTS + BASELINE_FILE)
barAs_2010 = baseline_2010['barAs'].values
barAu_2010 = baseline_2010['barAu'].values
barTs_2010 = baseline_2010['barTs'].values
barTu_2010 = baseline_2010['barTu'].values

# Compute 2100 levels from hats
barAs_2100 = barAs_2010 * barAs_hat.flatten()
barAu_2100 = barAu_2010 * barAu_hat.flatten()
barTs_2100 = barTs_2010 * barTs_hat.flatten()
barTu_2100 = barTu_2010 * barTu_hat.flatten()

# ===========================================================================
# Create Figure 8: Structural Parameters, 2100 vs 2010
# ===========================================================================
def plot_scatter(ax, x, y, title, xlabel='2010', ylabel='2100'):
    ax.scatter(x, y, s=20, alpha=0.5, color='gray')
    
    # Add 45-degree line
    min_val = min(x.min(), y.min())
    max_val = max(x.max(), y.max())
    ax.plot([min_val, max_val], [min_val, max_val], 'b-', linewidth=1.5)
    
    ax.set_xlabel(xlabel, fontsize=10)
    ax.set_ylabel(ylabel, fontsize=10)
    ax.set_title(title, fontsize=11)
    ax.grid(True, alpha=0.3)
    
    # Remove top and right spines (make it L-shaped!)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

# Create 2x2 figure
fig, axes = plt.subplots(2, 2, figsize=(12, 10))

# Top left: Amenity, Low-skilled
plot_scatter(axes[0, 0], barAu_2010, barAu_2100, 'Amenity, Low-skilled')

# Top right: Amenity, Skilled
plot_scatter(axes[0, 1], barAs_2010, barAs_2100, 'Amenity, Skilled')

# Bottom left: Productivity, Low-skilled
plot_scatter(axes[1, 0], barTu_2010, barTu_2100, 'Productivity, Low-skilled')

# Bottom right: Productivity, Skilled
plot_scatter(axes[1, 1], barTs_2010, barTs_2100, 'Productivity, Skilled')

plt.tight_layout()
plt.savefig(OUTPUTS2 + 'structural_parameters.png', dpi=300, bbox_inches='tight')
plt.show()

print("Figure saved successfully!")