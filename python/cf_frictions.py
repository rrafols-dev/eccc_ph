import numpy as np
import pandas as pd

def compute_cost_matrices(data_dist, sameisland, sameprov, diff_lat, diff_long, 
                          beta_hat, betas_hat, betau_hat, islands, islandu, 
                          provs, provu, lats, latu, longs, longu, homes, homeu,
                          sigma, theta, migsimple, I):
    """
    Compute bilateral cost matrices for trade (T) and migration (Ms, Mu)
    """

    
    # Rescale coefficients by their exponentiated parameters
    beta_hat_scaled = beta_hat / (sigma - 1)
    betas_hat_scaled = betas_hat / theta
    betau_hat_scaled = betau_hat / theta
    provs_scaled = provs / theta
    provu_scaled = provu / theta
    lats_scaled = lats / theta
    latu_scaled = latu / theta
    longs_scaled = longs / theta
    longu_scaled = longu / theta
    homes_scaled = homes / theta
    homeu_scaled = homeu / theta
    islands_scaled = islands / theta
    islandu_scaled = islandu / theta
    
    # Trade cost matrix
    T = np.exp(beta_hat_scaled * np.arcsinh(data_dist))
    
    # Migration cost matrices
    if migsimple == 0:
        # Skilled migration with full controls
        Ms = np.exp(
            betas_hat_scaled * np.arcsinh(data_dist)
            + islands_scaled * sameisland
            + provs_scaled * sameprov
            + lats_scaled * diff_lat
            + longs_scaled * diff_long
            + homes_scaled * I
        )
        
        # Unskilled migration with full controls
        Mu = np.exp(
            betau_hat_scaled * np.arcsinh(data_dist)
            + islandu_scaled * sameisland
            + provu_scaled * sameprov
            + latu_scaled * diff_lat
            + longu_scaled * diff_long            
            + homeu_scaled * I
        )
        
    elif migsimple == 1:
        # Simple specification with distance only
        Ms = np.exp(betas_hat_scaled * np.arcsinh(data_dist))
        Mu = np.exp(betau_hat_scaled * np.arcsinh(data_dist))
    
    # Rescale matrices to have 1 on diagonals for faster convergence
    Ms = Ms / Ms.diagonal()[:, None]
    Mu = Mu / Mu.diagonal()[:, None]
    T = T / T.diagonal()[:, None]
    
    return T, Ms, Mu