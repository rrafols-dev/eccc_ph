import numpy as np
import pandas as pd

def run_simulation(params, universal_data, csvfiles, outputs, baseline_file, area_hat, output_file):
    """
    Run exact-hat counterfactual simulation for spatial economics model
    
    Args:
        params: Dictionary of model parameters
        universal_data: Dictionary containing pre-computed matrices (T, Ms, Mu, etc.)
        csvfiles: Path to input CSV files
        outputs: Path to output directory
        baseline_file: Name of baseline results file
        slr_file: Name of sea level rise area file
        output_file: Name for output results file
    """
    
    sigma = params['sigma']
    theta_u = params['theta_u']
    theta_s = params['theta_s']
    kappa = params['kappa']
    eta = params['eta']
    mus = params['mus']
    muu = params['muu']
    N = params['N']


    id_2 = universal_data['id_2']    
    I = universal_data['I']
    T = universal_data['T']
    Ms = universal_data['Ms']
    Mu = universal_data['Mu']
    barAs_hat = universal_data['barAs_hat']
    barAu_hat = universal_data['barAu_hat']
    barTs_hat = universal_data['barTs_hat']
    barTu_hat = universal_data['barTu_hat']    
    

    print(f"Loading baseline: {baseline_file}")
    base = pd.read_csv(outputs + baseline_file)
    
    Ws = base['Ws'].to_numpy().reshape(-1, 1)
    Wu = base['Wu'].to_numpy().reshape(-1, 1)
    ws = base['ws'].to_numpy().reshape(-1, 1)
    wu = base['wu'].to_numpy().reshape(-1, 1)
    wages = base['wages'].to_numpy().reshape(-1, 1)
    wageu = base['wageu'].to_numpy().reshape(-1, 1)
    As = base['As'].to_numpy().reshape(-1, 1)
    Au = base['Au'].to_numpy().reshape(-1, 1)
    Ts = base['Ts'].to_numpy().reshape(-1, 1)
    Tu = base['Tu'].to_numpy().reshape(-1, 1)
    P = base['Price'].to_numpy().reshape(-1, 1)
    p = base['price'].to_numpy().reshape(-1, 1)
    Ls = base['Ls'].to_numpy().reshape(-1, 1)
    Lu = base['Lu'].to_numpy().reshape(-1, 1)
    Lso = base['Lso'].to_numpy().reshape(-1, 1)
    Luo = base['Luo'].to_numpy().reshape(-1, 1)
    xi = base['xi'].to_numpy().reshape(-1, 1)
    Ls_cfx = base['Ls_pred'].to_numpy().reshape(-1, 1)
    Lu_cfx = base['Lu_pred'].to_numpy().reshape(-1, 1)
    

    # ---------------------------------------------------------------------------
    # SET SHOCKS
    # ---------------------------------------------------------------------------       
    Ms_hat = Ms
    Mu_hat = Mu
    T_hat = T
    
    # ---------------------------------------------------------------------------
    # BASELINE POWER TERMS (INVARIANT TO LOOP UPDATING)
    # ---------------------------------------------------------------------------
    Ws_theta = Ws ** theta_s
    Wu_theta = Wu ** theta_u
    ws_theta = ws ** theta_s
    wu_theta = wu ** theta_u
    
    Ls = Ls_cfx 
    Lu = Lu_cfx
    
    den_LsL = Ms @ (Lso / Ws_theta)
    den_LuL = Mu @ (Luo / Wu_theta)
    P_den = T @ (p ** (1 - sigma))
    
    # ---------------------------------------------------------------------------
    # INITIALIZE ENDOGENOUS HATS
    # ---------------------------------------------------------------------------
    
    wages_hat = np.ones((N, 1))
    wageu_hat = np.ones((N, 1))
    P_hat = np.ones((N, 1))
    p_hat = np.ones((N, 1))
    Ws_hat = np.ones((N, 1))
    Wu_hat = np.ones((N, 1))
    Ls_hat = np.ones((N, 1))
    Lu_hat = np.ones((N, 1))
    ws_hat = np.ones((N, 1))
    wu_hat = np.ones((N, 1))
    
    print(f"Area shock: min={np.min(area_hat):.4f}, max={np.max(area_hat):.4f}")
    print(f"Locations unaffected: {np.sum(area_hat >= 0.99)}")
    print(f"Locations losing >10% area: {np.sum(area_hat < 0.90)}")
    print(f"Locations losing >30% area: {np.sum(area_hat < 0.70)}")
    
    diff = 1
    tol = 1e-4
    max_iter = 5000000
    iteration = 0 

    
    print("Starting iteration...")

    while diff > tol and iteration < max_iter:
        print('Iter', iteration, '--', diff)
        
        # ----------------------------------------------------------------------
        # Update density hat
        # ----------------------------------------------------------------------
        dens_hat = (Ls_hat * Ls + Lu_hat * Lu) / (Ls + Lu)
        dens_hat = dens_hat / area_hat
    
        # ----------------------------------------------------------------------
        # Update total productivity and amenity hats (with agglomeration/congestion)
        # ----------------------------------------------------------------------
        #productivities 
        Ts_hat = barTs_hat * (dens_hat ** mus)
        Tu_hat = barTu_hat * (dens_hat ** muu)
        #amenities 
        As_hat = barAs_hat * (dens_hat ** eta)
        Au_hat = barAu_hat * (dens_hat ** eta)
    
        # ----------------------------------------------------------------------
        # Update price hat (CES unit cost function)
        # ----------------------------------------------------------------------
    
        p_hat_new = (
            xi * ((wages_hat / Ts_hat) ** (1 - kappa))
            + (1 - xi) * ((wageu_hat / Tu_hat) ** (1 - kappa))
        ) ** (1 / (1 - kappa))

        
        # ----------------------------------------------------------------------
        # Price Index
        # ---------------------------------------------------------------------- 
        P_num = T_hat @ ((p * p_hat_new) ** (1 - sigma))
        P_hat_new = (P_num / P_den) ** (1 / (1 - sigma))
    
        # ----------------------------------------------------------------------
        # Migration equilibrium: Update welfare hats
        # ----------------------------------------------------------------------
     
        # Frechet welfare aggregator hats ..omega_g = (A_g hat)*(w_g hat)/(P hat)
        x_s = (As_hat * wages_hat) / P_hat_new
        x_u = (Au_hat * wageu_hat) / P_hat_new
    
        Ws_hat_new = ((Ms_hat @ (ws_theta * (x_s ** theta_s))) / Ws_theta) ** (1 / theta_s)
        Wu_hat_new = ((Mu_hat @ (wu_theta * (x_u ** theta_u))) / Wu_theta) ** (1 / theta_u)
    
    
        # ----------------------------------------------------------------------
        # Reallocation of Labor
        # ----------------------------------------------------------------------
        # --- Population reallocation (from migration) ---
        num_Ls = Ms_hat @ (Lso * (Ws_hat_new ** (-theta_s)) / Ws_theta)
        num_Lu = Mu_hat @ (Luo * (Wu_hat_new ** (-theta_u)) / Wu_theta)
        
        # Destination attractiveness  
        omega_s_theta = x_s ** theta_s
        omega_u_theta = x_u ** theta_u
        
        if iteration == 0:
            print(f"omega_s_theta: min={np.min(omega_s_theta):.6f}, max={np.max(omega_s_theta):.6f}")
            print(f"x_s: min={np.min(x_s):.6f}, max={np.max(x_s):.6f}")
            print(f"num_Ls: min={np.min(num_Ls):.6f}, max={np.max(num_Ls):.6f}")
            print(f"den_LsL: min={np.min(den_LsL):.6f}, max={np.max(den_LsL):.6f}")
            print(f"dens_hat at iter 0: min={np.min(dens_hat):.6f}, max={np.max(dens_hat):.6f}")
        
        Ls_hat_new = omega_s_theta * (num_Ls / den_LsL)
        Lu_hat_new = omega_u_theta * (num_Lu / den_LuL)
        
        # Enforce closed economy: Σ_d L'_d = Σ_d L_d for each skill type
        # Skilled workers
        total_Ls_baseline = np.sum(Ls)
        total_Ls_cf = np.sum(Ls_hat_new * Ls)
        Ls_adjustment = total_Ls_baseline / total_Ls_cf
        Ls_hat_new = Ls_hat_new * Ls_adjustment
        
        # Unskilled workers  
        total_Lu_baseline = np.sum(Lu)
        total_Lu_cf = np.sum(Lu_hat_new * Lu)
        Lu_adjustment = total_Lu_baseline / total_Lu_cf
        Lu_hat_new = Lu_hat_new * Lu_adjustment
        
        # ----------------------------------------------------------------------
        # Update wage from FOC
        # ----------------------------------------------------------------------
    
        # Use Yhat from Production CES Production Fn and not Income [EQ12]
        Y_hat = (xi * ((Ts_hat * Ls_hat_new) ** ((kappa-1)/kappa))
        + (1-xi) * ((Tu_hat * Lu_hat_new) ** ((kappa-1)/kappa))) ** (kappa/(kappa-1))
        
        # Individual wage FOCs to back out 
        wages_hat_new = p_hat_new * (Y_hat ** (1/kappa)) * (Ts_hat ** ((kappa-1)/kappa)) * (Ls_hat_new ** (-1/kappa))
        wageu_hat_new = p_hat_new * (Y_hat ** (1/kappa)) * (Tu_hat ** ((kappa-1)/kappa)) * (Lu_hat_new ** (-1/kappa))
        
        diff = (
            np.linalg.norm(wages_hat_new - wages_hat)
            + np.linalg.norm(wageu_hat_new - wageu_hat)
            + np.linalg.norm(P_hat_new - P_hat)
            + np.linalg.norm(p_hat_new - p_hat)
            + np.linalg.norm(Ws_hat_new - Ws_hat)
            + np.linalg.norm(Wu_hat_new - Wu_hat)
            + np.linalg.norm(Ls_hat_new - Ls_hat)
            + np.linalg.norm(Lu_hat_new - Lu_hat)
        )
        
        wages_hat = wages_hat_new
        wageu_hat = wageu_hat_new
        P_hat = P_hat_new
        p_hat = p_hat_new
        Ws_hat = Ws_hat_new
        Wu_hat = Wu_hat_new
        Ls_hat = Ls_hat_new
        Lu_hat = Lu_hat_new
        
        iteration += 1
    
    if iteration == max_iter:
        print(f"Warning: Maximum iterations ({max_iter}) reached. diff = {diff:.2e}")
    else:
        print(f"Converged in {iteration} iterations. Final diff = {diff:.2e}")
    

    ws_hat = (As_hat * wages_hat) / P_hat_new
    wu_hat = (Au_hat * wageu_hat) / P_hat_new
    
    Y2 = (wages * Ls) + (wageu * Lu)
    Ys_hat = (Ls_hat * Ls) * (wages_hat * wages) / (wages * Ls) 
    Yu_hat = (Lu_hat * Lu) * (wageu_hat * wageu) / (wageu * Lu) 
    Y_num2 = (Ls_hat * Ls) * (wages_hat * wages) + (Lu_hat * Lu) * (wageu_hat * wageu)
    Y_hat2 = Y_num2 / Y2


    
    names = [
        "id_2","As_hat","Au_hat","barAs_hat","barAu_hat",
        "Ts_hat","Tu_hat","barTs_hat","barTu_hat",
        "ws_hat","wu_hat","Ws_hat","Wu_hat",
        "P_hat","p_hat","wages_hat","wageu_hat",
        "Ls_hat","Lu_hat","Ys_hat","Yu_hat","Y_hat", "Y_hat2"
    ]
    
    env = locals()
    cols = [np.asarray(env[n]).reshape(-1) for n in names]
    M = np.column_stack(cols)
    df = pd.DataFrame(M, columns=names)
    df.to_csv(outputs + output_file, index=False)
        
    # ---------------------------------------------------------------------------
    # COMPUTE BILATERAL MIGRATION SHARES
    # ---------------------------------------------------------------------------
    # Skilled
    # Compute final levels
    ws_final = ws * ws_hat
    wu_final = wu * wu_hat
    Ws_final = Ws * Ws_hat
    Wu_final = Wu * Wu_hat
    
    # Migration probabilities: pi_nd ∝ Ms_nd * (ws_d / Ws_n)^θ
    attract_s = (ws_final.reshape(1, -1) / Ws_final.reshape(-1, 1)) ** theta_s
    pi_s_numerator = Ms * attract_s
    pi_s = pi_s_numerator / pi_s_numerator.sum(axis=1, keepdims=True)
    
    attract_u = (wu_final.reshape(1, -1) / Wu_final.reshape(-1, 1)) ** theta_u
    pi_u_numerator = Mu * attract_u
    pi_u = pi_u_numerator / pi_u_numerator.sum(axis=1, keepdims=True)
    
    # Save
    pi_s_df = pd.DataFrame(pi_s, columns=[f'dest_{i}' for i in range(N)])
    pi_s_df.insert(0, 'origin', range(N))
    pi_s_df.to_csv(outputs + output_file.replace('.csv', '_migration_skilled.csv'), index=False)
    
    pi_u_df = pd.DataFrame(pi_u, columns=[f'dest_{i}' for i in range(N)])
    pi_u_df.insert(0, 'origin', range(N))
    pi_u_df.to_csv(outputs + output_file.replace('.csv', '_migration_unskilled.csv'), index=False)
    
    print(f"Migration matrices saved")
    