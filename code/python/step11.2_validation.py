import numpy as np
import pandas as pd
from cf_sims_exportM import run_simulation
from root import ROOT

# ===========================================================================
# Set universal parameters
# ===========================================================================
PARAMS = {
    'sigma': 4,
    'theta': 1.5,
    # 'theta_s': 1.2,
    # 'theta_u': 1.1, #1.05 if same [old]
    'theta_s': .95, 
    'theta_u': 1.0,
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
    'shock_amen': 0,  #counterfactual channels allowed to react to climate
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
climate_data = np.array(pd.read_csv(CSVFILES + 'climate.csv'))
climate_data10 = np.array(pd.read_csv(CSVFILES + 'climateRF_2010.csv'))
climate_data00 = np.array(pd.read_csv(CSVFILES + 'climateRF_2000.csv'))
climate_data10 = np.delete(climate_data10, 1597, axis=0)
climate_data00 = np.delete(climate_data00, 1597, axis=0)
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
# Calculate bilateral frictions
# ===========================================================================
print("Computing bilateral cost matrices...")


theta = PARAMS['theta'] 
theta_s = PARAMS['theta_s'] 
theta_u = PARAMS['theta_u'] 

beta_hat_scaled = PARAMS['beta_hat'] / (PARAMS['sigma'] - 1)

betas_hat_scaled = PARAMS['betas_hat'] / theta_s
Ms = (betas_hat_scaled * np.arcsinh(data_dist)
      + PARAMS['islands'] / theta_s * sameisland
      + PARAMS['provs'] / theta_s * sameprov
      + PARAMS['lats'] / theta_s * diff_lat
      + PARAMS['longs'] / theta_s * diff_long
      + PARAMS['homes'] / theta_s * I)
Ms = np.exp(Ms)
Ms = Ms / Ms.diagonal()[:, None]

betau_hat_scaled = PARAMS['betau_hat'] / theta_u
Mu = (betau_hat_scaled * np.arcsinh(data_dist)
      + PARAMS['islandu'] / theta_u * sameisland
      + PARAMS['provu'] / theta_u * sameprov
      + PARAMS['latu'] / theta_u * diff_lat
      + PARAMS['longu'] / theta_u * diff_long
      + PARAMS['homeu'] / theta_u * I)
Mu = np.exp(Mu)
Mu = Mu / Mu.diagonal()[:, None]

T = np.exp(beta_hat_scaled * np.arcsinh(data_dist))
T = T / T.diagonal()[:, None]

# ===========================================================================
# Calculate shocks to exogenous objects
# ===========================================================================

# Load 2010 fundamentals directly
base_2010 = pd.read_csv(OUTPUTS + 'resultspy_new_2010.csv')
barAs_new = base_2010['barAs'].to_numpy().reshape(-1, 1)
barAu_new = base_2010['barAu'].to_numpy().reshape(-1, 1)
barTs_new = base_2010['barTs'].to_numpy().reshape(-1, 1)
barTu_new = base_2010['barTu'].to_numpy().reshape(-1, 1)

# Load 2000 fundamentals for computing hats
base_2000 = pd.read_csv(OUTPUTS + 'resultspy_new_2000.csv')
barAs = base_2000['barAs'].to_numpy().reshape(-1, 1)
barAu = base_2000['barAu'].to_numpy().reshape(-1, 1)
barTs = base_2000['barTs'].to_numpy().reshape(-1, 1)
barTu = base_2000['barTu'].to_numpy().reshape(-1, 1)

# Compute hats (2010 / 2000)
barAs_hat = barAs_new / barAs
barAu_hat = barAu_new / barAu
barTs_hat = barTs_new / barTs
barTu_hat = barTu_new / barTu


# ===========================================================================
# Declare inputs for counterfactuals
# ===========================================================================
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



# ===========================================================================
# Define simulation scenarios
# ===========================================================================
# Define your simulation scenarios

simulationsSLR = [
    {
        'name': 'Scenario_col1',
        'slr_column': 1,
        'output_file': 'validateTHETA.csv'
    },
]


# ===========================================================================
# Run all simulations
# ===========================================================================

for sim in simulationsSLR:
    print(f"\n{'='*70}")
    print(f"Running simulation: {sim['name']}")
    print(f"{'='*70}")
    
    # Extract the area column for this simulation
    area_new = data_slr[:, [sim['slr_column']]]
    area_hat = area_new / area
    
    run_simulation(
        params=PARAMS,
        universal_data=UNIVERSAL_DATA,
        csvfiles=CSVFILES,
        outputs=OUTPUTS,
        baseline_file='resultspy_new_2000.csv',
        area_hat=area_hat,
        output_file=sim['output_file']
    )
    
    print(f"\nCompleted: {sim['name']}")

print("\n" + "="*70)
print("All simulations completed!")
print("="*70)