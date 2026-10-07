import numpy as np
import pandas as pd
import shutil
from cf_sims import run_simulation
from cf_frictions import compute_cost_matrices
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
    'migsimple': 0 #can be zero
}
#activate lowerbound and higherbound climate scenario by making a 1-0 combo:
hb = 1
lb = 1

#Universal Inputs
CSVFILES = ROOT + 'inputs_csv/'
OUTPUTS = ROOT + 'outputs_csv/mig_full/'
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

if PARAMS['migsimple']==0:
    OUTPUTS = ROOT + 'outputs_csv/mig_full/'


if PARAMS['migsimple']==0 and PARAMS['simYear']==2050:
    # Copy results 2010 baseline to new folder
    OUTPUTS2 = ROOT + 'outputs_csv/mig_full/2050/'
    source = OUTPUTS + BASELINE_FILE
    destination = OUTPUTS2 + BASELINE_FILE
    shutil.copy(source, destination)
    
    #update where outputs will be saved
    OUTPUTS = OUTPUTS2

if PARAMS['migsimple']==0 and PARAMS['simYear']==2100 and lb==1 and hb==0:
    # Copy results 2010 baseline to new folder
    OUTPUTS2 = ROOT + 'outputs_csv/mig_full/lowerbound/'
    source = OUTPUTS + BASELINE_FILE
    destination = OUTPUTS2 + BASELINE_FILE
    shutil.copy(source, destination)
    
    climate_data[:, 3] = climate_bounds['avetemp_lb'].values
    
    #update where outputs will be saved
    OUTPUTS = OUTPUTS2

if PARAMS['migsimple']==0 and PARAMS['simYear']==2100 and lb==0 and hb==1:
    # Copy results 2010 baseline to new folder
    OUTPUTS2 = ROOT + 'outputs_csv/mig_full/upperbound/'
    source = OUTPUTS + BASELINE_FILE
    destination = OUTPUTS2 + BASELINE_FILE
    shutil.copy(source, destination)
    
    climate_data[:, 3] = climate_bounds['avetemp_ub'].values
    
    #update where outputs will be saved
    OUTPUTS = OUTPUTS2




# ===========================================================================
# Calculate bilateral frictions
# ===========================================================================
print("Computing bilateral cost matrices...")
T, Ms, Mu = compute_cost_matrices(
    data_dist, sameisland, sameprov, diff_lat, diff_long,
    PARAMS['beta_hat'], PARAMS['betas_hat'], PARAMS['betau_hat'], 
    PARAMS['islands'], PARAMS['islandu'], PARAMS['provs'], PARAMS['provu'], 
    PARAMS['lats'], PARAMS['latu'], PARAMS['longs'], PARAMS['longu'], 
    PARAMS['homes'], PARAMS['homeu'], PARAMS['sigma'], PARAMS['theta'], 
    PARAMS['migsimple'], I,
)


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

if PARAMS['shock_amen']==0 and PARAMS['shock_prod']==0 :
    simulationsSLR = [
        {
            'name': 'Scenario_col1',
            'slr_column': 1,
            'output_file': 'cf_SLRonly_sim1.csv'
        },
        {
            'name': 'Scenario_col2',
            'slr_column': 2,
            'output_file': 'cf_SLRonly_sim2.csv'
        },
        {
            'name': 'Scenario_col3',
            'slr_column': 3,
            'output_file': 'cf_SLRonly_sim3.csv'
        },
        {
            'name': 'Scenario_col4',
            'slr_column': 4,
            'output_file': 'cf_SLRonly_sim4.csv'
        },
        {
            'name': 'Scenario_col5',
            'slr_column': 5,
            'output_file': 'cf_SLRonly_sim5.csv'
        },
        {
            'name': 'Scenario_col6',
            'slr_column': 6,
            'output_file': 'cf_SLRonly_sim6.csv'
        },
        {
            'name': 'Scenario_col7',
            'slr_column': 7,
            'output_file': 'cf_SLRonly_sim7.csv'
        },    
        {
            'name': 'Scenario_col8',
            'slr_column': 8,
            'output_file': 'cf_SLRonly_sim8.csv'
        },
    ]



if PARAMS['shock_amen']==1 and PARAMS['shock_prod']==1 :
    simulationsSLR = [
        {
            'name': 'Scenario_col1',
            'slr_column': 1,
            'output_file': 'cf_SLRwTemp_sim1.csv'
        },
        {
            'name': 'Scenario_col2',
            'slr_column': 2,
            'output_file': 'cf_SLRwTemp_sim2.csv'
        },
        {
            'name': 'Scenario_col3',
            'slr_column': 3,
            'output_file': 'cf_SLRwTemp_sim3.csv'
        },
        {
            'name': 'Scenario_col4',
            'slr_column': 4,
            'output_file': 'cf_SLRwTemp_sim4.csv'
        },
        {
            'name': 'Scenario_col5',
            'slr_column': 5,
            'output_file': 'cf_SLRwTemp_sim5.csv'
        },
        {
            'name': 'Scenario_col6',
            'slr_column': 6,
            'output_file': 'cf_SLRwTemp_sim6.csv'
        },
        {
            'name': 'Scenario_col7',
            'slr_column': 7,
            'output_file': 'cf_SLRwTemp_sim7.csv'
        },    
        {
            'name': 'Scenario_col8',
            'slr_column': 8,
            'output_file': 'cf_SLRwTemp_sim8.csv'
        },
    ]




if PARAMS['shock_amen']==0 and PARAMS['shock_prod']==1 :
    simulationsSLR = [
        {
            'name': 'Scenario_col1',
            'slr_column': 1,
            'output_file': 'cf_SLRwTemp_prodOnly_sim1.csv'
        },
        {
            'name': 'Scenario_col2',
            'slr_column': 2,
            'output_file': 'cf_SLRwTemp_prodOnly_sim2.csv'
        },
        {
            'name': 'Scenario_col3',
            'slr_column': 3,
            'output_file': 'cf_SLRwTemp_prodOnly_sim3.csv'
        },
        {
            'name': 'Scenario_col4',
            'slr_column': 4,
            'output_file': 'cf_SLRwTemp_prodOnly_sim4.csv'
        },
        {
            'name': 'Scenario_col5',
            'slr_column': 5,
            'output_file': 'cf_SLRwTemp_prodOnly_sim5.csv'
        },
        {
            'name': 'Scenario_col6',
            'slr_column': 6,
            'output_file': 'cf_SLRwTemp_prodOnly_sim6.csv'
        },
        {
            'name': 'Scenario_col7',
            'slr_column': 7,
            'output_file': 'cf_SLRwTemp_prodOnly_sim7.csv'
        },    
        {
            'name': 'Scenario_col8',
            'slr_column': 8,
            'output_file': 'cf_SLRwTemp_prodOnly_sim8.csv'
        },
    ]

if PARAMS['shock_amen']==1 and PARAMS['shock_prod']==0 :
    simulationsSLR = [
        {
            'name': 'Scenario_col1',
            'slr_column': 1,
            'output_file': 'cf_SLRwTemp_amenOnly_sim1.csv'
        },
        {
            'name': 'Scenario_col2',
            'slr_column': 2,
            'output_file': 'cf_SLRwTemp_amenOnly_sim2.csv'
        },
        {
            'name': 'Scenario_col3',
            'slr_column': 3,
            'output_file': 'cf_SLRwTemp_amenOnly_sim3.csv'
        },
        {
            'name': 'Scenario_col4',
            'slr_column': 4,
            'output_file': 'cf_SLRwTemp_amenOnly_sim4.csv'
        },
        {
            'name': 'Scenario_col5',
            'slr_column': 5,
            'output_file': 'cf_SLRwTemp_amenOnly_sim5.csv'
        },
        {
            'name': 'Scenario_col6',
            'slr_column': 6,
            'output_file': 'cf_SLRwTemp_amenOnly_sim6.csv'
        },
        {
            'name': 'Scenario_col7',
            'slr_column': 7,
            'output_file': 'cf_SLRwTemp_amenOnly_sim7.csv'
        },    
        {
            'name': 'Scenario_col8',
            'slr_column': 8,
            'output_file': 'cf_SLRwTemp_amenOnly_sim8.csv'
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
        baseline_file='resultspy_new_2010.csv',
        area_hat=area_hat,
        output_file=sim['output_file']
    )
    
    print(f"\nCompleted: {sim['name']}")

print("\n" + "="*70)
print("All simulations completed!")
print("="*70)
