import numpy as np
import pandas as pd

def compute_fundamental_shocks(shock_amen, shock_prod, simYear, 
                               tempAs, tempAu, tempTs, tempTu,
                               climate_data, outputs, baseline_file):
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
    
    # Apply amenity shocks
    if shock_amen == 1:
        barAs_new = np.exp((tempAs/100) * delta_temp) * barAs
        barAu_new = np.exp((tempAu/100) * delta_temp) * barAu
    
    # Apply productivity shocks
    if shock_prod == 1:
        barTs_new = np.exp((tempTs/100) * delta_temp) * barTs
        barTu_new = np.exp((tempTu/100) * delta_temp) * barTu
    
    # Compute hats
    barAs_hat = barAs_new / barAs
    barAu_hat = barAu_new / barAu
    barTs_hat = barTs_new / barTs
    barTu_hat = barTu_new / barTu
    
    return barAs_hat, barAu_hat, barTs_hat, barTu_hat