import numpy as np
import pandas as pd
from cf_sims_theta import run_simulation
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
    'migsimple': 0 # always 0 for the final set of results in the paper 
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
    PARAMS['migsimple'], I
)



# ===========================================================================
# Parameter sweep: theta
# ===========================================================================
SLR_COLUMN = 5
THETA_VALUES = np.arange(1.1, 1.9, .1) 

# Precompute constant area_hat (SLR column 5)
area_new = data_slr[:, [SLR_COLUMN]]
area_hat = area_new / area

def _run_one(args):
    seq, val_theta = args

    local_params = PARAMS.copy()
    local_params['theta'] = float(val_theta)
    
    print(f"\n[Run {seq}] Processing parameter values:")
    print(f"   theta = {local_params['theta']:.2f}")

    # Recompute cost matrices with new theta
    T_new, Ms_new, Mu_new = compute_cost_matrices(
        data_dist, sameisland, sameprov, diff_lat, diff_long,
        local_params['beta_hat'], local_params['betas_hat'], local_params['betau_hat'], 
        local_params['islands'], local_params['islandu'], local_params['provs'], local_params['provu'], 
        local_params['lats'], local_params['latu'], local_params['longs'], local_params['longu'], 
        local_params['homes'], local_params['homeu'], local_params['sigma'], local_params['theta'], 
        local_params['migsimple'], I
    )

    # Recompute hats for THIS (tempTu, tempTs)
    barAs_hat, barAu_hat, barTs_hat, barTu_hat = compute_fundamental_shocks(
        local_params['shock_amen'], local_params['shock_prod'], local_params['simYear'],
        local_params['tempAs'], local_params['tempAu'], local_params['tempTs'], local_params['tempTu'],
        climate_data, OUTPUTS, BASELINE_FILE,
    )

    # Build UNIVERSAL_DATA for this run (same variable name as before)
    UNIVERSAL_DATA = {
        'id_2': id_2,
        'I': I,
        'T': T_new,  # Use new T
        'Ms': Ms_new,  # Use new Ms
        'Mu': Mu_new,  # Use new Mu
        'barAs_hat': barAs_hat,
        'barAu_hat': barAu_hat,
        'barTs_hat': barTs_hat,
        'barTu_hat': barTu_hat,
        'data_slr': data_slr
    }

    out_name = f"cf_SLR_thetaRange_seq{seq}.csv"

    print(f"\n{'-'*70}")
    print(f"Run #{seq}: theta={val_theta:.2f}")
    print(f"{'-'*70}")

    run_simulation(
        params=local_params,
        universal_data=UNIVERSAL_DATA,
        csvfiles=CSVFILES,
        outputs=OUTPUTS,
        baseline_file=BASELINE_FILE,
        area_hat=area_hat,
        output_file=out_name
    )

    return out_name

if __name__ == "__main__":
    n = len(THETA_VALUES)
    pairs = [(i + 1, THETA_VALUES[i]) for i in range(n)]
    
    print(f"Running {len(pairs)} jobs sequentially...")
    
    # Sequential instead of parallel
    results = []
    for pair in pairs:
        result = _run_one(pair)
        results.append(result)

    print("\n" + "="*70)
    print("All simulations completed! Files written:")
    for r in results:
        print(" -", r)
    print("="*70)
