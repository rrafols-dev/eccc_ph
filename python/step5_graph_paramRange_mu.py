# process_sweep_mu_results.py
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import to_rgba
from root import ROOT

OUTPUTS2 = ROOT + 'tabFig/'
OUTPUTS = ROOT + 'outputs_csv/mig_full/'
CSVFILES = ROOT + 'inputs_csv/'

# Load baseline once
baseline = pd.read_csv(f"{OUTPUTS}resultspy_new_2010.csv")
baseline_Ls = baseline['Ls'].values
baseline_Lu = baseline['Lu'].values
totalL = baseline_Ls.sum() + baseline_Lu.sum()

# Parameter values for mu (agglomeration parameter is actually alpha in the model!)
MU_VALUES = np.arange(0.02, 0.0801, 0.005)
n = len(MU_VALUES)

results = []

for i in range(n):
    seq = i + 1
    
    file_pattern = f'cf_SLR_muRange_seq{seq}.csv'
    df = pd.read_csv(f"{OUTPUTS}{file_pattern}")
    
    # Weight welfare by baseline labor, then sum across locations
    Ws_agg = (baseline_Ls * df['Ws_hat']).sum()
    Wu_agg = (baseline_Lu * df['Wu_hat']).sum()
    W_agg = Ws_agg + Wu_agg
    
    # Weight output by total baseline labor, then sum
    Y_agg = (df['Y_hat'] * (baseline_Ls + baseline_Lu)).sum()
    
    # Rescale to per-capita
    results.append({
        'seq': seq,
        'mu': MU_VALUES[i],
        'W': W_agg / totalL,
        'Ws': Ws_agg / baseline_Ls.sum(),
        'Wu': Wu_agg / baseline_Lu.sum(),
        'Y': Y_agg / totalL
    })

# Save aggregated results
results_df = pd.DataFrame(results)
results_df.to_csv(f"{OUTPUTS}aggregate_sweep_mu_results.csv", index=False)

print(f"Processed {len(results)} simulations")
print(results_df.head(10))

# Load results
df = pd.read_csv(f"{OUTPUTS}aggregate_sweep_mu_results.csv")

# Calculate dynamic y-axis limits
y_min = df[['W', 'Ws', 'Wu']].min().min()
y_max = df[['W', 'Ws', 'Wu']].max().max()
y_margin = (y_max - y_min) * 0.05
y_limits = (y_min - y_margin, y_max + y_margin)

fig, ax = plt.subplots(1, 1, figsize=(8, 5))

data = df.sort_values('mu')

# Plot welfare measures - opaque lines, transparent markers
ax.plot(data['mu'], data['W'], marker='o', label='Aggregate', linewidth=2.2, 
        color='crimson', markersize=4, 
        markerfacecolor=to_rgba('crimson', 0.5),
        markeredgecolor=to_rgba('crimson', 0.5))
ax.plot(data['mu'], data['Ws'], marker='s', label='Skilled', linewidth=1.5, 
        color='navy', markersize=4,
        markerfacecolor=to_rgba('navy', 0.5),
        markeredgecolor=to_rgba('navy', 0.5))
ax.plot(data['mu'], data['Wu'], marker='^', label='Unskilled', linewidth=1.5, 
        color='cornflowerblue', markersize=4,
        markerfacecolor=to_rgba('cornflowerblue', 0.5),
        markeredgecolor=to_rgba('cornflowerblue', 0.5))

# Add reference line at 1.0
ax.axhline(y=1.0, color='gray', linestyle='--', alpha=0.5)


# Set fixed y-axis limits
ax.set_ylim(y_limits)
ax.set_ylim(0.82, 0.87)

# Remove top and right spines for "L" shape
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Labels
ax.set_xlabel('$\\alpha$', fontsize=14)
ax.set_ylabel('Welfare', fontsize=12)
ax.legend()

plt.tight_layout()
plt.savefig(f"{OUTPUTS2}welfare_paramMu.png", dpi=300, bbox_inches='tight')
plt.close()
print(f"Plot saved to: {OUTPUTS2}welfare_paramMu.png")

print("Plot saved!")
