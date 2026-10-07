import numpy as np
import pandas as pd
from cf_frictions import compute_cost_matrices
from cf_shocks_newCity import compute_fundamental_shocks
from cf_sims_breaks import run_simulation
import os
from root import ROOT

# ===========================================================================
# Bootstrap setup
# ===========================================================================
bootstrap_draws = pd.read_csv('bootstrap_draws.csv')
B = len(bootstrap_draws)

print(f"Running {B} bootstrap iterations...")
print("="*70)

# ===========================================================================
# Set policy flags 
# ===========================================================================
# Run twice as as {1,0} and {0,1}
coastProtect = 1  # SET: 0 or 1
newCity = 0      # SET: 0 or 1


# If newCity=1, set boost factor
BOOST_FACTOR = 1.5  # SET: 1.5, 2.0, or 10.0 (only used if newCity=1)
BOOST_IDS = [1227, 1536, 1235, 1241]
### END HERE: MANUAL REVIEW ###

# ===========================================================================
# Set universal parameters (non-climate stay fixed)
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

# Universal Inputs
CSVFILES = ROOT + 'inputs_csv/'
OUTPUTS_BASE = ROOT + 'outputs_csv/mig_full/'
BASELINE_FILE = 'resultspy_new_2010.csv'

climate_data = np.array(pd.read_csv(CSVFILES + 'climate.csv'))
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


if coastProtect == 1 and newCity == 0:
    scenario_name = 'SLRwTemp_coastProtect'
elif coastProtect == 0 and newCity == 1:
    scenario_name = f'SLRwTemp_newCity_boost{int(BOOST_FACTOR*10)}'
else:
    raise ValueError("Invalid policy combination")

print(f"\nScenario: {scenario_name}")
print("="*70)

# Create bootstrap folder
OUTPUTS = OUTPUTS_BASE + f'bootstrap/{scenario_name}/'
os.makedirs(OUTPUTS, exist_ok=True)

# Copy baseline file
import shutil
source = OUTPUTS_BASE + BASELINE_FILE
destination = OUTPUTS + BASELINE_FILE
shutil.copy(source, destination)

# ===========================================================================
# Calculate bilateral frictions
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

Ms_use = Ms
Mu_use = Mu

# ===========================================================================
# Only sim5
# ===========================================================================
SIM5_COLUMN = 5
area_new = data_slr[:, [SIM5_COLUMN]]
area_hat = area_new / area

# Coastal protection: no area loss for these locations
if coastProtect == 1:
    mask = np.isin(id_2, [9991, 9992, 9993])
    area_hat[mask] = 1.0

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
    
    # Create PARAMS for this iteration
    PARAMS = PARAMS_BASE.copy()
    PARAMS['tempAs'] = eta_As * 100
    PARAMS['tempAu'] = eta_Au * 100
    PARAMS['tempTs'] = eta_Ts * 100
    PARAMS['tempTu'] = eta_Tu * 100
    
    # Compute fundamental shocks
    if newCity == 1:
        barAs_hat, barAu_hat, barTs_hat, barTu_hat = compute_fundamental_shocks(
            PARAMS['shock_amen'], PARAMS['shock_prod'], PARAMS['simYear'],
            PARAMS['tempAs'], PARAMS['tempAu'], PARAMS['tempTs'], PARAMS['tempTu'],
            climate_data, OUTPUTS, BASELINE_FILE, id_2,
            boost_ids=BOOST_IDS,
            boost_factor=BOOST_FACTOR
        )
    else:
        barAs_hat, barAu_hat, barTs_hat, barTu_hat = compute_fundamental_shocks(
            PARAMS['shock_amen'], PARAMS['shock_prod'], PARAMS['simYear'],
            PARAMS['tempAs'], PARAMS['tempAu'], PARAMS['tempTs'], PARAMS['tempTu'],
            climate_data, OUTPUTS, BASELINE_FILE, id_2
        )
    
    # If coastal protection, no climate shocks for protected locations
    if coastProtect == 1:
        mask = np.isin(id_2, [9991, 9992, 9993])
        barAs_hat[mask] = 1.0
        barAu_hat[mask] = 1.0
        barTs_hat[mask] = 1.0
        barTu_hat[mask] = 1.0
    
    # Prepare universal data
    UNIVERSAL_DATA = {
        'id_2': id_2,    
        'I': I,
        'T': T,
        'Ms': Ms_use,
        'Mu': Mu_use,
        'barAs_hat': barAs_hat,
        'barTs_hat': barTs_hat,
        'barAu_hat': barAu_hat,
        'barTu_hat': barTu_hat,
        'data_slr': data_slr
    }
    
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

successful = df_results['status'] == 'success'
print(f"\nSuccessful: {successful.sum()}/{B}")
print(f"Failed: {(~successful).sum()}/{B}")

if (~successful).sum() > 0:
    print("\nFailed iterations:")
    print(df_results[~successful][['draw', 'status']])

print(f"\nSimulation outputs saved to: {OUTPUTS}")
print(f"Parameter tracking saved to: {OUTPUTS}bootstrap_parameters.csv")
print("\nNext step: Run welfare calculation script to process all CSVs")
