import numpy as np
import pandas as pd

def compute_fundamental_shocks(shock_amen, shock_prod, simYear, 
                               tempAs, tempAu, tempTs, tempTu,
                               climate_data, outputs, baseline_file, id_2, boost_ids=None, boost_factor=1.0):
    """
    Compute fundamental shocks to amenities and productivity from climate change
    
    Returns:
        barAs_hat, barAu_hat, barTs_hat, barTu_hat
    """
    base = pd.read_csv(outputs + baseline_file)
    barAs = base['barAs'].to_numpy().reshape(-1, 1)
    barAu = base['barAu'].to_numpy().reshape(-1, 1)
    barTs = base['barTs'].to_numpy().reshape(-1, 1)
    barTu = base['barTu'].to_numpy().reshape(-1, 1)
    
    temp_baseline = climate_data[:, [1]]
    temp_2050 = climate_data[:, [2]]
    temp_2100 = climate_data[:, [3]]
    pr_baseline = climate_data[:, [7]]
    pr_2050 = climate_data[:, [8]]
    pr_2100 = climate_data[:, [9]]
    
    # Select climate scenario
    if simYear == 2100:
        delta_temp = temp_2100 - temp_baseline
        delta_pr = pr_2100 - pr_baseline
    elif simYear == 2050:
        delta_temp = temp_2050 - temp_baseline
        delta_pr = pr_2050 - pr_baseline
    
    # Initialize new fundamentals
    barAs_new = barAs.copy()
    barAu_new = barAu.copy()
    barTs_new = barTs.copy()
    barTu_new = barTu.copy()
    
    # Apply 1.5x boost to target districts (capped at district 9991)
    if boost_ids is not None:
        # Get ceiling from district 9991
        idx_9991 = np.where(id_2 == 9991)[0][0]
        ceiling_As = 1 * barAs_new[idx_9991][0]
        ceiling_Au = 1 * barAu_new[idx_9991][0]
        ceiling_Ts = 1 * barTs_new[idx_9991][0]
        ceiling_Tu = 1 * barTu_new[idx_9991][0]
        
        # Boost target districts
        for district_id in boost_ids:
            idx = np.where(id_2 == district_id)[0]
            if len(idx) > 0:
                idx = idx[0]
                
                 
                # Store original
                orig_As = barAs_new[idx][0]
                orig_Au = barAu_new[idx][0]
                orig_Ts = barTs_new[idx][0]
                orig_Tu = barTu_new[idx][0]
                
                
                barAs_new[idx] = min(barAs_new[idx] * boost_factor, ceiling_As)
                barAu_new[idx] = min(barAu_new[idx] * boost_factor, ceiling_Au)
                barTs_new[idx] = min(barTs_new[idx] * boost_factor, ceiling_Ts)
                barTu_new[idx] = min(barTu_new[idx] * boost_factor, ceiling_Tu)
                
                # Print transformation
                print(f"\nDistrict {district_id}: (Ceiling: As={ceiling_As:.6f}, Au={ceiling_Au:.6f}, Ts={ceiling_Ts:.6f}, Tu={ceiling_Tu:.6f})")
                print(f"  As: {orig_As:.6f} → {barAs_new[idx][0]:.6f} {'(CAPPED)' if orig_As*boost_factor> ceiling_As else ''}")
                print(f"  Au: {orig_Au:.6f} → {barAu_new[idx][0]:.6f} {'(CAPPED)' if orig_Au*boost_factor > ceiling_Au else ''}")
                print(f"  Ts: {orig_Ts:.6f} → {barTs_new[idx][0]:.6f} {'(CAPPED)' if orig_Ts*boost_factor > ceiling_Ts else ''}")
                print(f"  Tu: {orig_Tu:.6f} → {barTu_new[idx][0]:.6f} {'(CAPPED)' if orig_Tu*boost_factor > ceiling_Tu else ''}")
                
    
    # Apply amenity shocks
    if shock_amen == 1:
        barAs_new = np.exp((tempAs/100) * delta_temp) * barAs_new
        barAu_new = np.exp((tempAu/100) * delta_temp) * barAu_new
    
    # Apply productivity shocks
    if shock_prod == 1:
        barTs_new = np.exp((tempTs/100) * delta_temp) * barTs_new
        barTu_new = np.exp((tempTu/100) * delta_temp) * barTu_new
    
    # Compute hats
    barAs_hat = barAs_new / barAs
    barAu_hat = barAu_new / barAu
    barTs_hat = barTs_new / barTs
    barTu_hat = barTu_new / barTu
    
    return barAs_hat, barAu_hat, barTs_hat, barTu_hat