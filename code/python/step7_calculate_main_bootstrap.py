import numpy as np
import pandas as pd
import shutil
from cf_sims import run_simulation
from cf_frictions import compute_cost_matrices
from cf_shocks import compute_fundamental_shocks
from root import ROOT

# ===========================================================================
# Bootstrap setup
# ===========================================================================
# Read bootstrap draws (created from your Stata estimates)
bootstrap_draws = pd.read_csv('bootstrap_draws.csv')
B = len(bootstrap_draws)  # Number of bootstrap iterations

print(f"Running {B} bootstrap iterations...")
print("="*70)

# ===========================================================================
# Set universal parameters (non-climate parameters stay fixed)
# ===========================================================================
PARAMS_BASE = {
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
    'shock_amen': 1,
    'shock_prod': 1,
    'simYear': 2100,
    'migsimple': 0
}

# Climate scenario setup: run twice as {1,0} and {0,1} 
hb = 1  # SET THIS: 1 for upperbound
lb = 0  # SET THIS: 1 for lowerbound

# Universal Inputs
CSVFILES = ROOT + 'inputs_csv/'
OUTPUTS_BASE = ROOT + 'outputs_csv/'
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
N = PARAMS_BASE['N']
I = np.eye(N)

# Set output directory based on scenario
OUTPUTS = OUTPUTS_BASE + 'mig_full/'
scenario_name = 'SLRwTemp'


# Create scenario-specific bootstrap folder
if PARAMS_BASE['migsimple']==0 and PARAMS_BASE['simYear']==2100 and lb==1 and hb==0:
    OUTPUTS2 = OUTPUTS + f'lowerbound/bootstrap/{scenario_name}/'
    climate_data[:, 3] = climate_bounds['avetemp_lb'].values
    OUTPUTS = OUTPUTS2
elif PARAMS_BASE['migsimple']==0 and PARAMS_BASE['simYear']==2100 and lb==0 and hb==1:
    OUTPUTS2 = OUTPUTS + f'upperbound/bootstrap/{scenario_name}/'
    climate_data[:, 3] = climate_bounds['avetemp_ub'].values
    OUTPUTS = OUTPUTS2
else:
    OUTPUTS = OUTPUTS + f'bootstrap/{scenario_name}/'

# Create bootstrap output directory
import os
os.makedirs(OUTPUTS, exist_ok=True)

# Copy baseline file
source = OUTPUTS_BASE + 'mig_full/' + BASELINE_FILE
destination = OUTPUTS + BASELINE_FILE
shutil.copy(source, destination)

# ===========================================================================
# Calculate bilateral frictions (same for all bootstrap iterations)
# ===========================================================================
print("Computing bilateral cost matrices...")
T, Ms, Mu = compute_cost_matrices(
    data_dist, sameisland, sameprov, diff_lat, diff_long,
    PARAMS_BASE['beta_hat'], PARAMS_BASE['betas_hat'], PARAMS_BASE['betau_hat'], 
    PARAMS_BASE['islands'], PARAMS_BASE['islandu'], PARAMS_BASE['provs'], PARAMS_BASE['provu'], 
    PARAMS_BASE['lats'], PARAMS_BASE['latu'], PARAMS_BASE['longs'], PARAMS_BASE['longu'], 
    PARAMS_BASE['homes'], PARAMS_BASE['homeu'], PARAMS_BASE['sigma'], PARAMS_BASE['theta'], 
    PARAMS_BASE['migsimple'], I,
)

# ===========================================================================
# Only sim5 (the one you care about)
# ===========================================================================
SIM5_COLUMN = 5

# Storage for bootstrap results
bootstrap_results = []

# ===========================================================================
# Bootstrap loop
# ===========================================================================
for b in range(B):
    if b % 10 == 0:
        print(f"\nBootstrap iteration {b}/{B}")
    
    # Get climate elasticities for this draw
    eta_As = bootstrap_draws.loc[b, 'eta_As']
    eta_Au = bootstrap_draws.loc[b, 'eta_Au']
    eta_Ts = bootstrap_draws.loc[b, 'eta_Ts']
    eta_Tu = bootstrap_draws.loc[b, 'eta_Tu']
    
    # Create PARAMS for this iteration (climate elasticities in percentage points)
    PARAMS = PARAMS_BASE.copy()
    PARAMS['tempAs'] = eta_As * 100  # Convert to percentage
    PARAMS['tempAu'] = eta_Au * 100
    PARAMS['tempTs'] = eta_Ts * 100
    PARAMS['tempTu'] = eta_Tu * 100
    
    # Compute fundamental shocks with these parameters
    barAs_hat, barAu_hat, barTs_hat, barTu_hat = compute_fundamental_shocks(
        PARAMS['shock_amen'], PARAMS['shock_prod'], PARAMS['simYear'],
        PARAMS['tempAs'], PARAMS['tempAu'], PARAMS['tempTs'], PARAMS['tempTu'],
        climate_data, OUTPUTS, BASELINE_FILE
    )
    
    # Prepare universal data
    UNIVERSAL_DATA = {
        'id_2': id_2,    
        'I': I,
        'T': T,
        'Ms': Ms,
        'Mu': Mu,
        'barAs_hat': barAs_hat,
        'barTs_hat': barTs_hat,
        'barAu_hat': barAu_hat,
        'barTu_hat': barTu_hat,
        'data_slr': data_slr
    }
    
    # Extract sim5 area
    area_new = data_slr[:, [SIM5_COLUMN]]
    area_hat = area_new / area
    
    # Run simulation
    output_file = f'cf_{scenario_name}_sim5_boot{b:04d}.csv'
    
    try:
        run_simulation(
            params=PARAMS,
            universal_data=UNIVERSAL_DATA,
            csvfiles=CSVFILES,
            outputs=OUTPUTS,
            baseline_file=BASELINE_FILE,
            area_hat=area_hat,
            output_file=output_file
        )
        
        # Just track parameters for this draw
        bootstrap_results.append({
            'draw': b,
            'eta_As': eta_As,
            'eta_Au': eta_Au,
            'eta_Ts': eta_Ts,
            'eta_Tu': eta_Tu,
            'status': 'success'
        })
        
    except Exception as e:
        print(f"Error in bootstrap iteration {b}: {e}")
        bootstrap_results.append({
            'draw': b,
            'eta_As': eta_As,
            'eta_Au': eta_Au,
            'eta_Ts': eta_Ts,
            'eta_Tu': eta_Tu,
            'status': f'failed: {str(e)}'
        })

# ===========================================================================
# Save bootstrap tracking file
# ===========================================================================
df_results = pd.DataFrame(bootstrap_results)
df_results.to_csv(OUTPUTS + 'bootstrap_parameters.csv', index=False)

print("\n" + "="*70)
print("Bootstrap simulations complete!")
print("="*70)

# Summary
successful = df_results['status'] == 'success'
print(f"\nSuccessful: {successful.sum()}/{B}")
print(f"Failed: {(~successful).sum()}/{B}")

if (~successful).sum() > 0:
    print("\nFailed iterations:")
    print(df_results[~successful][['draw', 'status']])

print(f"\nSimulation outputs saved to: {OUTPUTS}")
print(f"Parameter tracking saved to: {OUTPUTS}bootstrap_parameters.csv")
print("\nNext step: Run welfare calculation script to process all CSVs")
