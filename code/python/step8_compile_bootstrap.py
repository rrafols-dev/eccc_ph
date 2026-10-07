import numpy as np
import pandas as pd
import glob
import matplotlib.pyplot as plt
import seaborn as sns
from root import ROOT

# ===========================================================================
# Configuration
# ===========================================================================
OUTPUTS_BASE = ROOT + 'outputs_csv/mig_full/'
BASELINE_FILE = ROOT + 'outputs_csv/mig_full/u_resultspy_new_2010.csv'
OUTPUTS2 = ROOT + 'tabFig/'

# Define all scenarios to analyze
SCENARIOS = [    
    {'name': 'SLRwTemp', 'shock_amen': 1, 'shock_prod': 1, 'label': 'Baseline: All Damages'},
    {'name': 'SLRwTemp_coastProtect', 'shock_amen': 1, 'shock_prod': 1, 'label': 'Baseline + Coastal Protection'},
    {'name': 'SLRwTemp_newCity_boost15', 'shock_amen': 1, 'shock_prod': 1, 'label': 'Baseline + New Inland City'},
]

# ===========================================================================
# Helper function to calculate welfare from a single CSV
# ===========================================================================
def calculate_welfare_from_csv(csv_path, baseline_Ls, baseline_Lu, totalL):
    """Calculate welfare metrics from simulation output CSV"""
    df = pd.read_csv(csv_path)
    
    Ws_hat = df['Ws_hat'].values
    Wu_hat = df['Wu_hat'].values
    Y_hat = df['Y_hat'].values
    
    # Weight by baseline population
    Ws_weighted = Ws_hat * baseline_Ls
    Wu_weighted = Wu_hat * baseline_Lu
    Y_weighted = Y_hat * (baseline_Ls + baseline_Lu)
    
    # Sum and convert to per-capita
    Ws_total = Ws_weighted.sum()
    Wu_total = Wu_weighted.sum()
    W_total = Ws_total + Wu_total
    
    Ws = Ws_total / baseline_Ls.sum()
    Wu = Wu_total / baseline_Lu.sum()
    W = W_total / totalL
    Y = Y_weighted.sum() / totalL
    
    return W, Ws, Wu, Y

# ===========================================================================
# Load baseline data (common across all scenarios)
# ===========================================================================
print("Loading baseline data...")
baseline = pd.read_csv(BASELINE_FILE)
baseline_Ls = baseline['u_Lso'].values
baseline_Lu = baseline['u_Luo'].values
totalL = baseline_Ls.sum() + baseline_Lu.sum()

print(f"Total baseline population: {totalL:,.0f}")
print(f"  Skilled: {baseline_Ls.sum():,.0f}")
print(f"  Unskilled: {baseline_Lu.sum():,.0f}")

# ===========================================================================
# Process each scenario
# ===========================================================================
all_results = []

for scenario in SCENARIOS:
    scenario_name = scenario['name']
    print("\n" + "="*70)
    print(f"Processing: {scenario_name}")
    print("="*70)
    
    # Get point estimate (main result)
    point_file = OUTPUTS_BASE + f'cf_{scenario_name}_sim5.csv'
    print(f"Point estimate: {point_file}")
    
    try:
        W_point, Ws_point, Wu_point, Y_point = calculate_welfare_from_csv(
            point_file, baseline_Ls, baseline_Lu, totalL
        )
        point_available = True
        print(f"  W={W_point:.4f}, Ws={Ws_point:.4f}, Wu={Wu_point:.4f}, Y={Y_point:.4f}")
    except FileNotFoundError:
        print(f"  WARNING: Point estimate file not found!")
        point_available = False
        W_point = Ws_point = Wu_point = Y_point = None
    
    # Get bootstrap results
    bootstrap_folder = OUTPUTS_BASE + f'bootstrap/{scenario_name}/'
    bootstrap_files = sorted(glob.glob(bootstrap_folder + f'cf_{scenario_name}_sim5_boot*.csv'))
    
    print(f"Bootstrap files found: {len(bootstrap_files)}")
    
    if len(bootstrap_files) == 0:
        print(f"  WARNING: No bootstrap files found in {bootstrap_folder}")
        continue
    
    # Calculate welfare for each bootstrap draw
    bootstrap_W = []
    bootstrap_Ws = []
    bootstrap_Wu = []
    bootstrap_Y = []
    
    for i, boot_file in enumerate(bootstrap_files):
        if i % 50 == 0:
            print(f"  Processing bootstrap {i}/{len(bootstrap_files)}...")
        
        try:
            W, Ws, Wu, Y = calculate_welfare_from_csv(
                boot_file, baseline_Ls, baseline_Lu, totalL
            )
            bootstrap_W.append(W)
            bootstrap_Ws.append(Ws)
            bootstrap_Wu.append(Wu)
            bootstrap_Y.append(Y)
        except Exception as e:
            print(f"  Error in {boot_file}: {e}")
            continue
    
    print(f"  Successful bootstrap calculations: {len(bootstrap_W)}/{len(bootstrap_files)}")
    
    # Calculate statistics
    if len(bootstrap_W) > 0:
        bootstrap_W = np.array(bootstrap_W)
        bootstrap_Ws = np.array(bootstrap_Ws)
        bootstrap_Wu = np.array(bootstrap_Wu)
        bootstrap_Y = np.array(bootstrap_Y)
        
        # Store results
        all_results.append({
            'scenario': scenario_name,
            'label': scenario['label'],
            'W_point': W_point,
            'Ws_point': Ws_point,
            'Wu_point': Wu_point,
            'Y_point': Y_point,
            'W_mean': bootstrap_W.mean(),
            'W_ci_lower': np.percentile(bootstrap_W, 2.5),
            'W_ci_upper': np.percentile(bootstrap_W, 97.5),
            'Ws_mean': bootstrap_Ws.mean(),
            'Ws_ci_lower': np.percentile(bootstrap_Ws, 2.5),
            'Ws_ci_upper': np.percentile(bootstrap_Ws, 97.5),
            'Wu_mean': bootstrap_Wu.mean(),
            'Wu_ci_lower': np.percentile(bootstrap_Wu, 2.5),
            'Wu_ci_upper': np.percentile(bootstrap_Wu, 97.5),
            'Y_mean': bootstrap_Y.mean(),
            'Y_ci_lower': np.percentile(bootstrap_Y, 2.5),
            'Y_ci_upper': np.percentile(bootstrap_Y, 97.5),
        })

# ===========================================================================
# Create summary table
# ===========================================================================
df_results = pd.DataFrame(all_results)
print("\n" + "="*70)
print("SUMMARY TABLE")
print("="*70)
print(df_results.to_string(index=False))

# Save to CSV
output_file = OUTPUTS_BASE + 'bootstrap_summary_all_scenarios.csv'
df_results.to_csv(output_file, index=False)
print(f"\nSummary saved to: {output_file}")

# ===========================================================================
# Create interval plot - Single panel with all 4 outcomes per scenario
# ===========================================================================
print("\nCreating interval plot...")

fig, ax = plt.subplots(1, 1, figsize=(14, 8))
#fig.suptitle('Climate Change Impact: Point Estimates with 95% Bootstrap Confidence Intervals', fontsize=14, fontweight='bold')

# Define outcomes with colors and markers
outcomes = [
    ('Y', 'Aggregate Output', 'k', '^', 0.3),    
    ('Wu', 'Unskilled Welfare', 'cornflowerblue', 's', 0.1),
    ('Ws', 'Skilled Welfare', 'navy', 'o', -0.1),    
    ('W', 'Aggregate Welfare', 'crimson', 'D',- 0.3),

]

scenarios = df_results['label'].values
n_scenarios = len(scenarios)

# Create y-positions for scenarios
y_base = np.arange(n_scenarios)

# Plot each outcome type
for outcome_code, outcome_label, color, marker, offset in outcomes:
    points = df_results[f'{outcome_code}_point'].values
    ci_lower = df_results[f'{outcome_code}_ci_lower'].values
    ci_upper = df_results[f'{outcome_code}_ci_upper'].values
    
    # Offset y-positions for this outcome
    y_pos = y_base + offset
    
    # Plot confidence intervals as horizontal bars
    for i in range(len(scenarios)):
        ax.plot([ci_lower[i], ci_upper[i]], [y_pos[i], y_pos[i]], 
               color=color, linewidth=3, alpha=0.6)
    
    # Plot point estimates
    ax.scatter(points, y_pos, color=color, s=50, zorder=5, 
              marker=marker, linewidths=1.5, label=outcome_label)

# Add vertical line at 1.0 (no change)
ax.axvline(x=1.0, color='darkgray', linestyle='--', linewidth=1, alpha=0.5)

# Formatting
ax.set_yticks(y_base)
ax.set_yticklabels(scenarios, fontsize=13, fontweight='bold')
ax.set_xlabel('Counterfactual Ratio', fontsize=11)
ax.grid(axis='x', alpha=0.3)
ax.legend(loc='best', fontsize=13, frameon=True)

# Set y-limits to give space
ax.set_ylim(-0.5, n_scenarios - 0.5)

# Set x-axis limits
ax.set_xlim(0.75, 1.00)

# Place legend inside plot at top center (reverse order)
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles[::-1], labels[::-1], loc='upper center', fontsize=10, 
          frameon=False, fancybox=True, shadow=False, ncol=4)


# Add buffer at top for legend (use NEGATIVE value because axis is inverted)
ax.set_ylim(n_scenarios - 0.5, -0.7)

plt.tight_layout()

# Save plot
plot_file = OUTPUTS2 + 'bootstrap_intervals_all_scenarios.png'
plt.savefig(plot_file, dpi=300, bbox_inches='tight')
print(f"Plot saved to: {plot_file}")

plt.show()

print("\n" + "="*70)
print("All scenarios processed!")
print("="*70)
