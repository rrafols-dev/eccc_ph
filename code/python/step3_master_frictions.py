import numpy as np
import pandas as pd
from cf_sims_breaks import run_simulation
from cf_frictions import compute_cost_matrices
from cf_shocks import compute_fundamental_shocks
from root import ROOT

# ===========================================================================
# Set universal parameters
# ===========================================================================
no_migCost = 2 # 1 = only home bias and distance ; 2 = 5% Reduction cost
epsilon = .95 # 95 = 5 percent reduction 

    
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
N = PARAMS['N']
I = np.eye(N)



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


if no_migCost==1:
    T_costless, Ms_costless, Mu_costless = compute_cost_matrices(
        data_dist, sameisland, sameprov, diff_lat, diff_long,
        PARAMS['beta_hat'], PARAMS['betas_hat'], PARAMS['betau_hat'], 
        0, 0, 0, 0, 0, 0, 0, 0, 
        PARAMS['homes'], PARAMS['homeu'], PARAMS['sigma'], PARAMS['theta'], 
        PARAMS['migsimple'] , I,
    )


if no_migCost==2:
    Ms_costless2 = Ms ** epsilon
    Mu_costless2 = Mu ** epsilon
    np.fill_diagonal(Ms_costless2,1)
    np.fill_diagonal(Mu_costless2,1)
    
    



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
if no_migCost==2:
    UNIVERSAL_DATA = {
        'id_2': id_2,    
        'I': I,
        'T': T,
        'Ms': Ms_costless2,
        'Mu': Mu_costless2,
        'barAs_hat': barAs_hat,
        'barTs_hat': barTs_hat,
        'barAu_hat': barAu_hat,
        'barTu_hat': barTu_hat,
        'data_slr': data_slr
    }

if no_migCost==1:
    UNIVERSAL_DATA = {
        'id_2': id_2,    
        'I': I,
        'T': T,
        'Ms': Ms_costless,
        'Mu': Mu_costless,
        'barAs_hat': barAs_hat,
        'barTs_hat': barTs_hat,
        'barAu_hat': barAu_hat,
        'barTu_hat': barTu_hat,
        'data_slr': data_slr
    }
    

    


# ===========================================================================
# Define simulation scenarios (SLR + w/temp effects)
# ===========================================================================
# Define your simulation scenarios

if no_migCost==2 :
    simulationsSLR = [
        {
            'name': 'Scenario_col5',
            'slr_column': 5,
            'output_file': 'cf_SLRwTemp_reducedMig_sim5.csv'
        },
    ]


if no_migCost==1  :
    simulationsSLR = [
        {
            'name': 'Scenario_col5',
            'slr_column': 5,
            'output_file': 'cf_SLRwTemp_noMigCost_sim5.csv'
        },
    ]








# ===========================================================================
# Run main simulation 
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
    