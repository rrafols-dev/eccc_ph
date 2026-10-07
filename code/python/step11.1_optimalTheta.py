import numpy as np
import pandas as pd
from cf_sims_exportM import run_simulation
from root import ROOT

# ===========================================================================
# THETA ESTIMATION - MATCH MIGRATION RATES (WITH FOREIGN)
# ===========================================================================

CSVFILES = ROOT + 'inputs_csv/'
OUTPUTS = ROOT + 'outputs_csv/mig_full/'
N = 1600

BASE_PARAMS = {
    'sigma': 4,
    'kappa': 1.5,
    'eta': -0.05,
    'mus': 0.02,
    'muu': 0.02,
    'beta_hat_dist': .081,
    'betas_hat_dist': -1.926,
    'betau_hat_dist': -1.975,
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
    'simYear': 2010,
    'tempAs': -3.48,
    'tempAu': -0.3,
    'tempTs': -4.36,
    'tempTu': -9.30,
    'migsimple': 0
}

# ===========================================================================
# Load actual migration data
# foreign row = index 1597 (o_id==8887)
# foreign col = index 1597 (wa_skilled8887 is col 1598 in df, -1 for o_id = 1597 in numpy)
# ===========================================================================
print("="*70)
print("LOADING ACTUAL MIGRATION DATA")
print("="*70)

FOREIGN_ROW = 1597
FOREIGN_COL = 1597

migsk_2010_full = np.array(pd.read_csv(CSVFILES + 'Lod_sk_2010.csv'))[:, 1:]  # (1601, 1601)
migun_2010_full = np.array(pd.read_csv(CSVFILES + 'Lod_un_2010.csv'))[:, 1:]

# Save foreign col BEFORE deletion: emigrants from each origin to abroad, shape (1601,)
emigrants_sk_full = np.nan_to_num(migsk_2010_full[:, FOREIGN_COL].copy(), nan=0.0)
emigrants_un_full = np.nan_to_num(migun_2010_full[:, FOREIGN_COL].copy(), nan=0.0)

# Save foreign row BEFORE deletion: immigrants from abroad into each destination, shape (1601,)
immigrants_sk_full = np.nan_to_num(migsk_2010_full[FOREIGN_ROW, :].copy(), nan=0.0)
immigrants_un_full = np.nan_to_num(migun_2010_full[FOREIGN_ROW, :].copy(), nan=0.0)

# Drop foreign row and col -> domestic only (1600, 1600)
migsk_2010 = np.delete(np.delete(migsk_2010_full, FOREIGN_ROW, axis=0), FOREIGN_COL, axis=1)
migun_2010 = np.delete(np.delete(migun_2010_full, FOREIGN_ROW, axis=0), FOREIGN_COL, axis=1)

# Trim foreign vectors: drop the foreign entry itself -> shape (1600,)
emigrants_sk  = np.delete(emigrants_sk_full,  FOREIGN_ROW)
emigrants_un  = np.delete(emigrants_un_full,  FOREIGN_ROW)
immigrants_sk = np.delete(immigrants_sk_full, FOREIGN_COL)
immigrants_un = np.delete(immigrants_un_full, FOREIGN_COL)

# Full origin population = domestic row sum + emigrants to abroad
wa_sk = migsk_2010.sum(axis=1) + emigrants_sk   # shape (1600,)
wa_un = migun_2010.sum(axis=1) + emigrants_un

# Actual stayers
stayers_actual_sk = np.diag(migsk_2010)
stayers_actual_un = np.diag(migun_2010)

# Actual rates (foreign in numerators, wa in denominator)
outflows_actual_sk = migsk_2010.sum(axis=1) - stayers_actual_sk + emigrants_sk
outflows_actual_un = migun_2010.sum(axis=1) - stayers_actual_un + emigrants_un
inflows_actual_sk  = migsk_2010.sum(axis=0) - stayers_actual_sk + immigrants_sk
inflows_actual_un  = migun_2010.sum(axis=0) - stayers_actual_un + immigrants_un

outmig_rate_actual_sk = (outflows_actual_sk / wa_sk) * 100
outmig_rate_actual_un = (outflows_actual_un / wa_un) * 100
inmig_rate_actual_sk  = (inflows_actual_sk  / wa_sk) * 100
inmig_rate_actual_un  = (inflows_actual_un  / wa_un) * 100
stayrate_actual_sk    = (stayers_actual_sk   / wa_sk) * 100
stayrate_actual_un    = (stayers_actual_un   / wa_un) * 100

print(f"\nACTUAL RATES (mean %):")
print(f"  Skilled   - Stay: {stayrate_actual_sk.mean():.2f}, Outmig: {outmig_rate_actual_sk.mean():.2f}, Inmig: {inmig_rate_actual_sk.mean():.2f}")
print(f"  Unskilled - Stay: {stayrate_actual_un.mean():.2f}, Outmig: {outmig_rate_actual_un.mean():.2f}, Inmig: {inmig_rate_actual_un.mean():.2f}")

# ===========================================================================
# Load universal data
# ===========================================================================
def load_and_trim(fname):
    arr = np.array(pd.read_csv(CSVFILES + fname))
    arr = np.delete(arr, FOREIGN_ROW, axis=0)
    arr = np.delete(arr, FOREIGN_COL, axis=1)
    return arr[:, 1:]

data_dist  = load_and_trim('distance.csv')
diff_long  = load_and_trim('diff_long.csv')
diff_lat   = load_and_trim('diff_lat.csv')
sameprov   = load_and_trim('sameprov.csv')
sameisland = load_and_trim('sameisland.csv')

data_slr = np.array(pd.read_csv(CSVFILES + 'slr_area.csv'))
id_2 = data_slr[:, [0]]
I = np.eye(N)

# ===========================================================================
# Simulation function
# ===========================================================================

def simulate_at_theta_pair(theta_s, theta_u):
    print(f"\n{'─'*70}")
    print(f"Testing θ_s = {theta_s:.2f}, θ_u = {theta_u:.2f}")
    print(f"{'─'*70}")

    params = BASE_PARAMS.copy()
    params['theta'] = theta_s
    params['theta_s'] = theta_s
    params['theta_u'] = theta_u

    beta_hat_scaled = params['beta_hat'] / (params['sigma'] - 1)

    betas_hat_scaled = params['betas_hat'] / theta_s
    Ms = (betas_hat_scaled * np.arcsinh(data_dist)
          + params['islands'] / theta_s * sameisland
          + params['provs']   / theta_s * sameprov
          + params['lats']    / theta_s * diff_lat
          + params['longs']   / theta_s * diff_long
          + params['homes']   / theta_s * I)
    Ms = np.exp(Ms)
    Ms = Ms / Ms.diagonal()[:, None]

    betau_hat_scaled = params['betau_hat'] / theta_u
    Mu = (betau_hat_scaled * np.arcsinh(data_dist)
          + params['islandu'] / theta_u * sameisland
          + params['provu']   / theta_u * sameprov
          + params['latu']    / theta_u * diff_lat
          + params['longu']   / theta_u * diff_long
          + params['homeu']   / theta_u * I)
    Mu = np.exp(Mu)
    Mu = Mu / Mu.diagonal()[:, None]

    T = np.exp(beta_hat_scaled * np.arcsinh(data_dist))
    T = T / T.diagonal()[:, None]

    base_2010 = pd.read_csv(OUTPUTS + 'resultspy_new_2010.csv')
    base_2000 = pd.read_csv(OUTPUTS + 'resultspy_new_2000.csv')

    barAs_hat = base_2010['barAs'].values.reshape(-1,1) / base_2000['barAs'].values.reshape(-1,1)
    barAu_hat = base_2010['barAu'].values.reshape(-1,1) / base_2000['barAu'].values.reshape(-1,1)
    barTs_hat = base_2010['barTs'].values.reshape(-1,1) / base_2000['barTs'].values.reshape(-1,1)
    barTu_hat = base_2010['barTu'].values.reshape(-1,1) / base_2000['barTu'].values.reshape(-1,1)

    universal_data = {
        'id_2': id_2, 'I': I, 'T': T, 'Ms': Ms, 'Mu': Mu,
        'barAs_hat': barAs_hat, 'barAu_hat': barAu_hat,
        'barTs_hat': barTs_hat, 'barTu_hat': barTu_hat,
        'data_slr': data_slr
    }

    run_simulation(
        params=params,
        universal_data=universal_data,
        csvfiles=CSVFILES,
        outputs=OUTPUTS,
        baseline_file='resultspy_new_2000.csv',
        area_hat=np.ones((N, 1)),
        output_file=f'theta_test_s{theta_s:.2f}_u{theta_u:.2f}.csv'
    )

    # Load predicted shares (1600 x 1600, domestic only)
    predsk = pd.read_csv(OUTPUTS + f'theta_test_s{theta_s:.2f}_u{theta_u:.2f}_migration_skilled.csv').iloc[:, 1:].values
    predun = pd.read_csv(OUTPUTS + f'theta_test_s{theta_s:.2f}_u{theta_u:.2f}_migration_unskilled.csv').iloc[:, 1:].values

    # Scale by full origin population
    predsk_flows = predsk * wa_sk[:, np.newaxis]
    predun_flows = predun * wa_un[:, np.newaxis]

    stayers_pred_sk = np.diag(predsk_flows)
    stayers_pred_un = np.diag(predun_flows)

    # Correct rates: add back foreign in numerators
    dom_out_sk = predsk_flows.sum(axis=1) - stayers_pred_sk
    dom_out_un = predun_flows.sum(axis=1) - stayers_pred_un
    dom_in_sk  = predsk_flows.sum(axis=0) - stayers_pred_sk
    dom_in_un  = predun_flows.sum(axis=0) - stayers_pred_un

    outmig_rate_pred_sk = ((dom_out_sk + emigrants_sk)  / wa_sk) * 100
    outmig_rate_pred_un = ((dom_out_un + emigrants_un)  / wa_un) * 100
    inmig_rate_pred_sk  = ((dom_in_sk  + immigrants_sk) / wa_sk) * 100
    inmig_rate_pred_un  = ((dom_in_un  + immigrants_un) / wa_un) * 100
    stayrate_pred_sk    = (stayers_pred_sk / wa_sk) * 100
    stayrate_pred_un    = (stayers_pred_un / wa_un) * 100

    # MAE on rates (in percentage points)
    def mape(actual, pred):
        return np.mean(np.abs(pred - actual))

    def safe_corr(a, b):
        mask = np.isfinite(a) & np.isfinite(b)
        if mask.sum() < 2: return np.nan
        return np.corrcoef(a[mask], b[mask])[0, 1]

    err_stay_sk   = mape(stayrate_actual_sk,    stayrate_pred_sk)
    err_stay_un   = mape(stayrate_actual_un,    stayrate_pred_un)
    err_outmig_sk = mape(outmig_rate_actual_sk, outmig_rate_pred_sk)
    err_outmig_un = mape(outmig_rate_actual_un, outmig_rate_pred_un)
    err_inmig_sk  = mape(inmig_rate_actual_sk,  inmig_rate_pred_sk)
    err_inmig_un  = mape(inmig_rate_actual_un,  inmig_rate_pred_un)

    weighted_error = 0.5 * (err_outmig_sk + err_inmig_sk) / 2 \
                   + 0.5 * (err_outmig_un + err_inmig_un) / 2

    corr_stay_sk   = safe_corr(stayrate_actual_sk,    stayrate_pred_sk)
    corr_stay_un   = safe_corr(stayrate_actual_un,    stayrate_pred_un)
    corr_outmig_sk = safe_corr(outmig_rate_actual_sk, outmig_rate_pred_sk)
    corr_outmig_un = safe_corr(outmig_rate_actual_un, outmig_rate_pred_un)
    corr_inmig_sk  = safe_corr(inmig_rate_actual_sk,  inmig_rate_pred_sk)
    corr_inmig_un  = safe_corr(inmig_rate_actual_un,  inmig_rate_pred_un)

    print(f"  SKILLED (rates %):")
    print(f"    Stay   - Actual: {stayrate_actual_sk.mean():.2f}, Pred: {stayrate_pred_sk.mean():.2f}, MAPE: {err_stay_sk:.2f}%, Corr: {corr_stay_sk:.4f}")
    print(f"    Outmig - Actual: {outmig_rate_actual_sk.mean():.2f}, Pred: {outmig_rate_pred_sk.mean():.2f}, MAPE: {err_outmig_sk:.2f}%, Corr: {corr_outmig_sk:.4f}")
    print(f"    Inmig  - Actual: {inmig_rate_actual_sk.mean():.2f}, Pred: {inmig_rate_pred_sk.mean():.2f}, MAPE: {err_inmig_sk:.2f}%, Corr: {corr_inmig_sk:.4f}")
    print(f"  UNSKILLED (rates %):")
    print(f"    Stay   - Actual: {stayrate_actual_un.mean():.2f}, Pred: {stayrate_pred_un.mean():.2f}, MAPE: {err_stay_un:.2f}%, Corr: {corr_stay_un:.4f}")
    print(f"    Outmig - Actual: {outmig_rate_actual_un.mean():.2f}, Pred: {outmig_rate_pred_un.mean():.2f}, MAPE: {err_outmig_un:.2f}%, Corr: {corr_outmig_un:.4f}")
    print(f"    Inmig  - Actual: {inmig_rate_actual_un.mean():.2f}, Pred: {inmig_rate_pred_un.mean():.2f}, MAPE: {err_inmig_un:.2f}%, Corr: {corr_inmig_un:.4f}")
    print(f"  WEIGHTED MAPE: {weighted_error:.2f}%")

    return {
        'err_stay_sk': err_stay_sk,   'err_stay_un': err_stay_un,
        'err_outmig_sk': err_outmig_sk, 'err_outmig_un': err_outmig_un,
        'err_inmig_sk': err_inmig_sk,   'err_inmig_un': err_inmig_un,
        'weighted_error': weighted_error,
        'corr_stay_sk': corr_stay_sk,   'corr_stay_un': corr_stay_un,
        'corr_outmig_sk': corr_outmig_sk, 'corr_outmig_un': corr_outmig_un,
        'corr_inmig_sk': corr_inmig_sk,   'corr_inmig_un': corr_inmig_un,
        'mean_stay_pred_sk': stayrate_pred_sk.mean(),
        'mean_stay_pred_un': stayrate_pred_un.mean(),
        'mean_outmig_pred_sk': outmig_rate_pred_sk.mean(),
        'mean_outmig_pred_un': outmig_rate_pred_un.mean(),
        'mean_inmig_pred_sk': inmig_rate_pred_sk.mean(),
        'mean_inmig_pred_un': inmig_rate_pred_un.mean(),
    }

# ===========================================================================
# Grid search
# ===========================================================================
print("\n" + "="*70)
print("GRID SEARCH - MATCH MIGRATION RATES (WITH FOREIGN)")
print("="*70)

theta_values = np.array([0.7, 0.8, 0.9, 0.95, 1.0, 1.1, 1.2, 1.3, 1.4, 1.5])

results_list = []

for theta_s in theta_values:
    for theta_u in theta_values:
        if theta_s == 1.5 and theta_u == 1.5:
            print("Skipping θ_s=1.50, θ_u=1.50 (known divergence)")
            continue
        print(f"\n>>> θ_s={theta_s:.2f}, θ_u={theta_u:.2f}")
        try:
            result = simulate_at_theta_pair(theta_s, theta_u)
            result['theta_s'] = theta_s
            result['theta_u'] = theta_u
            result['converged'] = True
            results_list.append(result)
        except Exception as e:
            print(f"  FAILED: {str(e)}")
            results_list.append({
                'theta_s': theta_s, 'theta_u': theta_u,
                'weighted_error': np.nan, 'converged': False
            })

# ===========================================================================
# Results
# ===========================================================================
print("\n" + "="*70)
print("RESULTS")
print("="*70)

results_df = pd.DataFrame(results_list)
converged = results_df[results_df['converged'] == True].copy()

if len(converged) > 0:
    best_idx = converged['weighted_error'].idxmin()
    best = converged.loc[best_idx]

    print(f"\nOPTIMAL PARAMETERS:")
    print(f"  θ_s = {best['theta_s']:.2f}")
    print(f"  θ_u = {best['theta_u']:.2f}")
    print(f"  Weighted MAPE: {best['weighted_error']:.2f}%")

    print("\n" + "="*70)
    print("TOP 10 BEST FITS")
    print("="*70)
    print("\nθ_s   | θ_u   | MAE    | Out_s(A/P)    | In_s(A/P)     | Out_u(A/P)    | In_u(A/P)")
    print("------|-------|--------|---------------|---------------|---------------|---------------")
    top10 = converged.nsmallest(10, 'weighted_error')
    for _, row in top10.iterrows():
        print(f"{row['theta_s']:5.2f} | {row['theta_u']:5.2f} | {row['weighted_error']:6.3f} | "
              f"{outmig_rate_actual_sk.mean():5.2f}/{row['mean_outmig_pred_sk']:5.2f}      | "
              f"{inmig_rate_actual_sk.mean():5.2f}/{row['mean_inmig_pred_sk']:5.2f}      | "
              f"{outmig_rate_actual_un.mean():5.2f}/{row['mean_outmig_pred_un']:5.2f}      | "
              f"{inmig_rate_actual_un.mean():5.2f}/{row['mean_inmig_pred_un']:5.2f}")
else:
    print("No converged results!")

results_df.to_csv(OUTPUTS + 'theta_rates_calibration.csv', index=False)
print(f"\nSaved to: {OUTPUTS}theta_rates_calibration.csv")
print("="*70)