import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from root import ROOT

##  Set parameters
sigma = 4   # goods substitution
theta = 1.5 # frechet parameter
kappa = 1.5 # skill substitutability
#trade beta for distance
beta_hat = - 1.199
#migration beta for distance (asinh specification)
betas_hat = -1.171
betau_hat = -1.315
islands = 0.274
islandu = 0.145
provs = 1.37
provu = 1.146
homes = 2.911 
homeu = 2.579
longs = 0.164
longu = 0.291
lats = -.02
latu = -.098
#congestion parameter:
eta = -.05
#chauvin et al india estimate for agglomeration parameter (alpha in the paper, not differenatiated by skill):
mus = .02
muu = .02
N = 1600
vgamma=(theta-1)/theta
I = np.eye(N)


def norm01(data):
    return (data - np.min(data)) / (np.max(data) - np.min(data))

def norm(data):
    return np.divide(data,np.mean(data))

def noForeign(path):
    df = pd.read_csv(path)
    df = df[df.iloc[:, 0] != 8887]     # only filter by first column
    return df.to_numpy()



#SET INVERSION PROCEDURE YEAR
year = 2010 # run this twice, once for 2000, another for 2010

# ===========================================================================
#  Step #1: Loading data
# ===========================================================================
OUTPUTS = ROOT + 'outputs_csv/mig_full/'
CSVFILES = ROOT + 'inputs_csv/'

data_dist = np.array(pd.read_csv(CSVFILES + 'distance.csv'))
data_dist = np.delete(data_dist, 1597, axis=0)   # remove the row (foreign)
data_dist = np.delete(data_dist, 1598, axis=1)   # remove the column (foreign)
data_dist = data_dist[:, 1:]

# diff_long
diff_long = np.array(pd.read_csv(CSVFILES + 'diff_long.csv'))
diff_long = np.delete(diff_long, 1597, axis=0)   # remove the row (foreign)
diff_long = np.delete(diff_long, 1598, axis=1)   # remove the column (foreign)
diff_long = diff_long[:, 1:]

# diff_lat
diff_lat = np.array(pd.read_csv(CSVFILES + 'diff_lat.csv'))
diff_lat = np.delete(diff_lat, 1597, axis=0)
diff_lat = np.delete(diff_lat, 1598, axis=1)
diff_lat = diff_lat[:, 1:]

# sameprov
sameprov = np.array(pd.read_csv(CSVFILES + 'sameprov.csv'))
sameprov = np.delete(sameprov, 1597, axis=0)
sameprov = np.delete(sameprov, 1598, axis=1)
sameprov = sameprov[:, 1:]

# sameisland
sameisland = np.array(pd.read_csv(CSVFILES + 'sameisland.csv'))
sameisland = np.delete(sameisland, 1597, axis=0)
sameisland = np.delete(sameisland, 1598, axis=1)
sameisland = sameisland[:, 1:]


data_Y2010      = noForeign(CSVFILES + 'localGDP.csv')
data_L2010      = noForeign(CSVFILES + 'Ld_2010.csv')
data_Lo2010     = noForeign(CSVFILES + 'Lo_2010.csv')
data_wage2010   = noForeign(CSVFILES + 'wagesresidSector_2010.csv')
data_L2000      = noForeign(CSVFILES + 'Ld_2000.csv')
data_Lo2000     = noForeign(CSVFILES + 'Lo_2000.csv')
data_wage2000   = noForeign(CSVFILES + 'wagesresidSector_2000.csv')
data_slr = np.array(pd.read_csv(CSVFILES + 'slr_area.csv'))
area = data_slr[:,[1]]
GDP = data_Y2010[:,[4]] #current PPP
popcell = data_Y2010[:,[5]]

if year==2010:
    resultsCSV= OUTPUTS+'resultspy_new_2010.csv'
    resultsCSV2= OUTPUTS+'u_resultspy_new_2010.csv'
    id_2 = data_L2010[:,[0]]    
    Ls = data_L2010[:,[2]]
    Lu = data_L2010[:,[3]]
    L = data_L2010[:,[1]]
    Lso = data_Lo2010[:,[2]]
    Luo = data_Lo2010[:,[3]]
    Lo = data_Lo2010[:,[1]]
    wages = data_wage2010[:,[2]]
    wageu = data_wage2010[:,[3]]



    
if year==2000:
    resultsCSV= OUTPUTS+'resultspy_new_2000.csv'
    resultsCSV2= OUTPUTS+'u_resultspy_new_2000.csv'
    id_2 = data_L2000[:,[0]]    
    Ls = data_L2000[:,[2]]
    Lu = data_L2000[:,[3]]
    L = data_L2000[:,[1]]
    Lso = data_Lo2000[:,[2]]
    Luo = data_Lo2000[:,[3]]
    Lo = data_Lo2000[:,[1]]
    wages = data_wage2000[:,[2]]
    wageu = data_wage2000[:,[3]]
    tfp_scale_factor = 1.1874    


    


# ===========================================================================
#  Step #2: Preparing matrices
# ===========================================================================
#extract L, w from dataframe
#calculate skilled-labor ratio for rescaling
Lsshare=np.divide(np.sum(Ls),np.sum(L))
Lushare=np.divide(np.sum(Lu),np.sum(L))

#GDP (rescale to match total L-count in each municipality)
GDP = ((GDP*1000000000)/popcell)*L
GDP = norm(GDP)

#income
Y = (np.multiply(Ls,wages)) + (np.multiply(Lu,wageu))

Ys = (np.multiply(Ls,wages))
Yu = (np.multiply(Lu,wageu))
dens = np.divide(L,area)

#Normalize to mean
dens=norm(dens)
Y=norm(Y)
Ys = norm(Ys)
Yu = norm(Yu)


# Normalize both skill groups by the same common mean
wage_mean = np.mean(np.concatenate([wages, wageu]))
labor_mean = np.mean(np.concatenate([Ls, Lu]))
mean_Lu_o = np.mean(np.concatenate([Luo, Lso]))
wages = wages / wage_mean
wageu = wageu / wage_mean
Ls = Ls / labor_mean
Lu = Lu / labor_mean
Luo = Luo / mean_Lu_o
Lso = Lso / mean_Lu_o


# ===========================================================================
#  Step #2.1: Matrics of distance costs
# ===========================================================================
#rescale coefficients by their exponentiatied parameters
beta_hat = beta_hat / (sigma - 1)
betas_hat, betau_hat, provs, provu, lats, latu, longs, longu, homeu, homes, islands, islandu = [x / theta for x in [betas_hat, betau_hat, provs, provu, lats, latu, longs, longu, homeu, homes, islands, islandu]]


T = data_dist
T = np.exp(beta_hat * np.arcsinh(data_dist))
Ms = data_dist
Mu = data_dist

# constructing distance matrix
Ms = (betas_hat * np.arcsinh(data_dist)
+ islands * sameisland
+ provs * sameprov
+ lats * diff_lat
+ longs * diff_long
+ homes * I
)

Mu = (betau_hat * np.arcsinh(data_dist)
+ islandu * sameisland
+ provu * sameprov
+ latu * diff_lat
+ longu * diff_long
+ homeu * I
)


#recale matrix to have 1 diagonals for slightly faster convergence (not by much)
Mu = np.exp(Mu)
Ms = np.exp(Ms)
Ms = Ms / Ms.diagonal()[:, None] 
Mu = Mu / Mu.diagonal()[:, None] 
T = T / T.diagonal()[:, None] 


# ===========================================================================
#  Step #3: Setting-up guesses
# ===========================================================================
Tu = np.ones((N,1))
Ts = np.ones((N,1))
Ws = np.ones((N,1))
ws = np.ones((N,1))
Wu = np.ones((N,1))
wu = np.ones((N,1))
P = np.ones((N,1))
p = np.ones((N,1))
diff1 = 1
diff2 = 1
tol = 1e-09
alpha_update = 1
iter=1

# ===========================================================================
#  Welfare and Prices
# ===========================================================================

while diff1 > tol:
    print(iter, '--', diff1)
    
    #Ws_new = W^\theta
    Ws_new = Ms @ (1 / ws)                 
    Wu_new = Mu @ (1 / wu)  
    
    #ws_new = \omega^-\theta
    ws_new = np.sum([np.multiply(Ms, np.matmul(1/Ls, np.transpose(np.multiply(Lso,np.power(Ws,-1)))))], axis=2)   
    ws_new = ws_new.T #take reciprocal

    wu_new = np.sum([np.multiply(Mu, np.matmul(1/Lu, np.transpose(np.multiply(Luo,np.power(Wu,-1)))))], axis=2)
    wu_new = wu_new.T

    #storing unnormalized welfare
    u_ws = ws_new
    u_wu = wu_new
    
    #normalizing ...
    ws_new = norm(ws_new)
    wu_new = norm(wu_new)
    
    # this is ps_new = p_s^(\sigma-1)
    p_new = np.sum([
        np.multiply(T, np.matmul(1/Y, np.transpose(Y * (P**(sigma-1)))))], axis=2)
    p_new = (p_new.T) 
    
    
    #normalizing 
    u_p = p_new
    p_new = norm(p_new)
    P_new = T @ (p_new**(-1))         # this is P^(1−\sigma)



    # ---------------------------------------------------------------------------
    #  Calculating Tolerance 
    # ---------------------------------------------------------------------------
    dA = np.linalg.norm((ws - ws_new),2)
    dB = np.linalg.norm((Ws - Ws_new),2)
    dC = np.linalg.norm((Wu - Wu_new),2)
    dD = np.linalg.norm((wu - wu_new),2)
    eA = np.linalg.norm((p - p_new),2)
    eB = np.linalg.norm((P - P_new),2)
 
    diff1 = dA + dB + dC + dD +  eA + eB 
            
            
    # ---------------------------------------------------------------------------
    #  Updating welfare and prices
    # ---------------------------------------------------------------------------        
    
    iter = iter + 1       
    ws = alpha_update * ws_new + (1 - alpha_update) * ws
    Ws = alpha_update * Ws_new + (1 - alpha_update) * Ws
    wu = alpha_update * wu_new + (1 - alpha_update) * wu
    Wu = alpha_update * Wu_new + (1 - alpha_update) * Wu

    p = alpha_update * p_new + (1 - alpha_update) * p
    P = alpha_update * P_new + (1 - alpha_update) * P
    

# ---------------------------------------------------------------------------
#  Recovering Amenities 
# (it's wages that give us the source of variation in these results)
# ---------------------------------------------------------------------------
  
Ws = Ws**(1/theta)   
Wu = Wu**(1/theta)   
ws = ws**(-1/theta)   
wu = wu**(-1/theta) 
p = p**(1/(sigma-1))
P = P**(1/(1 - sigma))


#Exogenous amenities from indirect utility fn
As = (P*ws) / wages
Au = (P*wu) / wageu
barAs = As / (dens**eta) 
barAu = Au / (dens**eta)


# ---------------------------------------------------------------------------
#  Recovering Productivities + calibrated CES-shares
# ---------------------------------------------------------------------------  
if year == 2000:
    data_2010 = np.array(pd.read_csv(OUTPUTS + 'resultspy_new_2010.csv'))
    xi = data_2010[:,[39]]
    Ts = data_2010[:,[12]] 
    Tu = data_2010[:,[13]] 
    
    # Scale 2010 productivities backward to 2000
    Ts = Ts / tfp_scale_factor
    Tu = Tu / tfp_scale_factor

    
if year == 2010:
    
    
    def solve_production_params_bounded(wages, wageu, Ls, Lu, GDP, p, kappa):
        """
        Solve with bounds to keep parameters positive
        """
        def equations(params):
            xi, Ts, Tu = params
            
            # Equation 1: FOC for high-skill
            eq1 = wages - p * xi * (Ts**((kappa-1)/kappa)) * (Ls**(-1/kappa)) * (GDP**(1/kappa))
            
            # Equation 2: FOC for low-skill  
            eq2 = wageu - p * (1-xi) * (Tu**((kappa-1)/kappa)) * (Lu**(-1/kappa)) * (GDP**(1/kappa))
            
            # Equation 3: Production function
            term_H = xi * (Ts * Ls)**((kappa-1)/kappa)
            term_L = (1-xi) * (Tu * Lu)**((kappa-1)/kappa)
            GDP_pred = (term_H + term_L)**(kappa/(kappa-1))
            eq3 = GDP - GDP_pred
            
            return [eq1, eq2, eq3]
        
        # Initial guess
        x0 = [0.5, 1.0, 1.0]
        
        # Bounds: xi in (0.01, 0.99), As and Au > 0.001
        bounds = ([0.01, 0.001, 0.001], [0.99, 100, 100])
        
        # Solve with bounds
        result = least_squares(equations, x0, bounds=bounds)
        
        return result.x, result.success
    
    # Loop with error handling
    results = []
    failed_count = 0
    
    
    # CREATE DATAFRAME FROM EXISTING ARRAYS
    df = pd.DataFrame({
        'wages': wages.flatten(),
        'wageu': wageu.flatten(),
        'Ls': Ls.flatten(),
        'Lu': Lu.flatten(),
        'output': GDP.flatten(),
        'prices': p.flatten()
    })
    
    
    for idx, row in df.iterrows():
        try:
            params, success = solve_production_params_bounded(
                wages=row['wages'],
                wageu=row['wageu'],
                Ls=row['Ls'],
                Lu=row['Lu'],
                GDP=row['output'],
                p=row['prices'],
                kappa=kappa
            )
            
            if success:
                xi, Ts, Tu = params
                results.append({
                    'municipality': idx,
                    'xi': xi,
                    'Ts': Ts,
                    'Tu': Tu,
                    'converged': True
                })
            else:
                failed_count += 1
                results.append({
                    'municipality': idx,
                    'xi': np.nan,
                    'Ts': np.nan,
                    'Tu': np.nan,
                    'converged': False
                })
        except Exception as e:
            failed_count += 1
            results.append({
                'municipality': idx,
                'xi': np.nan,
                'Ts': np.nan,
                'Tu': np.nan,
                'converged': False
            })
    
    
    print(f"\nSuccessfully solved: {len(results) - failed_count}/{len(results)} municipalities")
    print(f"Failed: {failed_count}")
    params_df = pd.DataFrame(results)
    
    # After the main loop, retry failed cases with different initial guesses
    failed_indices = params_df[~params_df['converged']]['municipality'].values
    
    if len(failed_indices) > 0:
        print(f"\n=== RETRYING {len(failed_indices)} FAILED MUNICIPALITIES ===")
        
        for failed_idx in failed_indices:
            row = df.iloc[failed_idx]
            
            # Try multiple different initial guesses
            initial_guesses = [
                [0.3, 1.0, 1.0],   # Lower xi
                [0.7, 1.0, 1.0],   # Higher xi
                [0.5, 0.5, 0.5],   # Lower productivities
                [0.5, 2.0, 2.0],   # Higher productivities
            ]
            
            converged = False
            for x0_try in initial_guesses:
                try:
                    result = least_squares(
                        lambda params: [
                            row['wages'] - row['prices'] * params[0] * (params[1]**((kappa-1)/kappa)) * (row['Ls']**(-1/kappa)) * (row['output']**(1/kappa)),
                            row['wageu'] - row['prices'] * (1-params[0]) * (params[2]**((kappa-1)/kappa)) * (row['Lu']**(-1/kappa)) * (row['output']**(1/kappa)),
                            row['output'] - (params[0]*(params[1]*row['Ls'])**((kappa-1)/kappa) + (1-params[0])*(params[2]*row['Lu'])**((kappa-1)/kappa))**(kappa/(kappa-1))
                        ],
                        x0_try,
                        bounds=([0.01, 0.001, 0.001], [0.99, 100, 100])
                    )
                    
                    if result.success:
                        xi_val, Ts_val, Tu_val = result.x
                        params_df.loc[failed_idx, 'xi'] = xi_val
                        params_df.loc[failed_idx, 'Ts'] = Ts_val
                        params_df.loc[failed_idx, 'Tu'] = Tu_val
                        params_df.loc[failed_idx, 'converged'] = True
                        print(f"Municipality {failed_idx} solved with x0={x0_try}")
                        converged = True
                        break
                except:
                    continue
            
            if not converged:
                print(f"Municipality {failed_idx} still failed after retries")
    
    
    # Extract calibrated values from params_df and convert back to numpy arrays
    xi = params_df['xi'].values.reshape(-1, 1)
    Ts = params_df['Ts'].values.reshape(-1, 1)  # Using As as skilled productivity
    Tu = params_df['Tu'].values.reshape(-1, 1)  # Using Au as unskilled productivity
    
    
    
    
    # === DIAGNOSTIC CODES HERE ===
    df['wage_bill_skilled'] = df['wages'] * df['Ls']
    df['wage_bill_unskilled'] = df['wageu'] * df['Lu']
    df['total_wage_bill'] = df['wage_bill_skilled'] + df['wage_bill_unskilled']
    df['skilled_share'] = df['wage_bill_skilled'] / df['total_wage_bill']
    
    print("=== WAGE BILL SHARE STATISTICS ===")
    print(df['skilled_share'].describe())
    
    # === CALIBRATED PARAMETERS DIAGNOSTIC ===
    print("\n" + "="*60)
    print("CALIBRATED PARAMETERS DIAGNOSTICS")
    print("="*60)
    
    # Basic statistics
    print("\n=== XI (CES SHARE PARAMETER) ===")
    print(params_df['xi'].describe())
    print(f"Min: {params_df['xi'].min():.6f}")
    print(f"Max: {params_df['xi'].max():.6f}")
    print(f"Range: {params_df['xi'].max() - params_df['xi'].min():.6f}")
    print(f"Hitting upper bound (>0.98): {(params_df['xi'] > 0.98).sum()}")
    print(f"Hitting lower bound (<0.02): {(params_df['xi'] < 0.02).sum()}")
    
    print("\n=== Ts (SKILLED PRODUCTIVITY) ===")
    print(params_df['Ts'].describe())
    print(f"Min: {params_df['Ts'].min():.6f}")
    print(f"Max: {params_df['Ts'].max():.6f}")
    print(f"Range: {params_df['Ts'].max() - params_df['Ts'].min():.6f}")
    print(f"Coefficient of variation: {params_df['Ts'].std() / params_df['Ts'].mean():.3f}")
    
    print("\n=== Tu (UNSKILLED PRODUCTIVITY) ===")
    print(params_df['Tu'].describe())
    print(f"Min: {params_df['Tu'].min():.6f}")
    print(f"Max: {params_df['Tu'].max():.6f}")
    print(f"Range: {params_df['Tu'].max() - params_df['Tu'].min():.6f}")
    print(f"Coefficient of variation: {params_df['Tu'].std() / params_df['Tu'].mean():.3f}")
    
    print("\n=== Ts/Tu RATIO ===")
    params_df['TsTu_ratio'] = params_df['Ts'] / params_df['Tu']
    print(params_df['TsTu_ratio'].describe())
    print(f"Min ratio: {params_df['TsTu_ratio'].min():.3f}")
    print(f"Max ratio: {params_df['TsTu_ratio'].max():.3f}")
    
    print("\n=== CONVERGENCE ===")
    print(f"Successfully converged: {params_df['converged'].sum()}/{len(params_df)}")
    print(f"Failed: {(~params_df['converged']).sum()}")
    
    print("\n=== TOP 5 MOST EXTREME Ts/Tu RATIOS ===")
    extreme = params_df.nlargest(5, 'TsTu_ratio')
    for idx in extreme['municipality'].values:
        if pd.notna(params_df.loc[idx, 'xi']):  # Skip NaN
            print(f"\nMunicipality {idx}:")
            print(f"  xi={params_df.loc[idx, 'xi']:.6f}, Ts={params_df.loc[idx, 'Ts']:.6f}, Tu={params_df.loc[idx, 'Tu']:.6f}")
            print(f"  Ts/Tu ratio: {params_df.loc[idx, 'Ts']/params_df.loc[idx, 'Tu']:.1f}")
            print(f"  Input wages: {df.iloc[idx]['wages']:.3f}, wageu: {df.iloc[idx]['wageu']:.3f}")
            print(f"  Input labor: Ls={df.iloc[idx]['Ls']:.3f}, Lu={df.iloc[idx]['Lu']:.3f}")
            print(f"  Wage bill share (skilled): {(df.iloc[idx]['wages']*df.iloc[idx]['Ls'])/(df.iloc[idx]['wages']*df.iloc[idx]['Ls']+df.iloc[idx]['wageu']*df.iloc[idx]['Lu']):.3f}")
    
    print("\n=== BOTTOM 5 (LOWEST Ts/Tu RATIOS) ===")
    low_extreme = params_df.nsmallest(5, 'TsTu_ratio')
    for idx in low_extreme['municipality'].values:
        if pd.notna(params_df.loc[idx, 'xi']):
            print(f"\nMunicipality {idx}:")
            print(f"  xi={params_df.loc[idx, 'xi']:.6f}, Ts={params_df.loc[idx, 'Ts']:.6f}, Tu={params_df.loc[idx, 'Tu']:.6f}")
            print(f"  Ts/Tu ratio: {params_df.loc[idx, 'Ts']/params_df.loc[idx, 'Tu']:.6f}")
            print(f"  Input wages: {df.iloc[idx]['wages']:.3f}, wageu: {df.iloc[idx]['wageu']:.3f}")
            print(f"  Input labor: Ls={df.iloc[idx]['Ls']:.3f}, Lu={df.iloc[idx]['Lu']:.3f}")
    # Convert to DataFrame
    params_df = pd.DataFrame(results)

# Exogenous productivities 
barTs = Ts / (dens**mus)
barTu = Tu / (dens**muu)

#normalize these objects to their mean to highlight inequalities (for visualization)
normAs=norm(As)
normTs=norm(Ts)
normTu=norm(Tu)
normAu=norm(Au)
normbarAs=norm(barAs)
normbarTs=norm(barTs)
normbarTu=norm(barTu)
normbarAu=norm(barAu)






# ---------------------------------------------------------------------------
#  Model predicted wages, output and Ls (should be rescaled for output)
# ---------------------------------------------------------------------------  
wages_cfx= (P*ws) / As
wageu_cfx= (P*wu) / Au

Ls_cfx = np.multiply(np.power(ws,theta),np.matmul(Ms, np.divide(Lso,np.power(Ws,theta))))
Lu_cfx = np.multiply(np.power(wu,theta),np.matmul(Mu, np.divide(Luo,np.power(Wu,theta))))




# ---------------------------------------------------------------------------
#  Manual calculation to check if  the same (robustness)
# ---------------------------------------------------------------------------  
manualWs = (Ms @ ((As * wages_cfx / P)**theta))**(1/theta)
manualWu = (Mu @ ((Au * wageu_cfx / P)**theta))**(1/theta)


# ---------------------------------------------------------------------------
#  Calculate unnormalized Prices, Welfare, Productivities, Amenities
# ---------------------------------------------------------------------------  
u_Ws = np.matmul(Ms, np.power(u_ws,theta))
u_Ws = np.power(u_Ws, 1/theta)
u_Wu = np.matmul(Mu, np.power(u_wu,theta))
u_Wu = np.power(u_Wu, 1/theta)

u_P = np.matmul(T, np.power(u_p,1-sigma))
u_P = np.power(u_P,(1/(1-sigma)))


#unnormalized labor given unnormalized welfare (verified match)
u_Lso = data_Lo2010[:,[2]]
u_Luo = data_Lo2010[:,[3]] 
u_Ls = data_L2010[:,[2]]
u_Lu = data_L2010[:,[3]]


#at baseline, wages are taken as given, we import them to obtain exogenous amenities and productivities:
if year==2000:
    u_wages = data_wage2000[:,[2]]
    u_wageu = data_wage2000[:,[3]]
if year==2010:
    u_wages = data_wage2010[:,[2]]
    u_wageu = data_wage2010[:,[3]]

u_As=np.divide(np.multiply(u_P,u_ws), u_wages)
u_Au=np.divide(np.multiply(u_P,u_wu), u_wageu)
u_barAs=np.divide(u_As, np.power(dens,eta))
u_barAu=np.divide(u_Au, np.power(dens,eta))


u_Ts= np.multiply(np.divide(u_wages,u_p), u_Lso)
u_Tu= np.multiply(np.divide(u_wageu,u_p), u_Luo)
u_barTs=np.divide(u_Ts, np.power(dens,mus))
u_barTu=np.divide(u_Tu, np.power(dens,muu))

#distribution of As/Au matters more
u_manualWs = np.power(np.matmul(Ms,np.power(np.divide(np.multiply(As,u_wages),u_P),theta)),1/theta)
u_manualWu = np.power(np.matmul(Mu,np.power(np.divide(np.multiply(Au,u_wageu),u_P),theta)),1/theta)





# ---------------------------------------------------------------------------
#  Exporting
# ---------------------------------------------------------------------------   
variables = [id_2, ws, Ws, Wu, wu, p, P, As, barAs, Au, barAu, Ts, \
             barTs, Tu, barTu, wages, wageu, Ls, Lu, manualWs, manualWu, \
            normAs, normAu, normTs, normTu, normbarAs, normbarAu, normbarTs, \
             normbarTu]
for x in variables:
    x=x[:,0]

id_2=id_2[:,0]
ws=ws[:,0]
Ws=Ws[:,0]
Wu=Wu[:,0]
wu=wu[:,0]
p=p[:,0]
P=P[:,0]
As=As[:,0]
barAs=barAs[:,0]
Au=Au[:,0]
barAu=barAu[:,0]
Ts=Ts[:,0]
barTs=barTs[:,0]
Tu=Tu[:,0]
barTu=barTu[:,0]
wages=wages[:,0]
wageu=wageu[:,0]
Ls=Ls[:,0]
Lu=Lu[:,0]
Lso=Lso[:,0]
Luo=Luo[:,0]
manualWs=manualWs[:,0]
manualWu=manualWu[:,0]
normAs=normAs[:,0]
normAu=normAu[:,0]
normTs=normTs[:,0]
normTu=normTu[:,0]
normbarAs=normbarAs[:,0]
normbarAu=normbarAu[:,0]
normbarTs=normbarTs[:,0]
normbarTu=normbarTu[:,0]
wages_pred = wages_cfx[:,0]
wageu_pred = wageu_cfx[:,0]
Ls_pred = Ls_cfx[:,0]
Lu_pred = Lu_cfx[:,0] 
Y=Y[:,0]
Yu=Yu[:,0]
Ys=Ys[:,0]
xi=xi[:,0]


export = pd.DataFrame({'id': id_2, 
                       'ws': ws, 'wu': wu, 'Ws': Ws, 'Wu': Wu,  
                       'price': p,  'Price': P, 'As':As,'Au':Au, 
                       'barAs':barAs, 'barAu':barAu, 
                       'Ts':Ts, 'Tu':Tu,
                       'barTs':barTs, 'barTu':barTu, 
                       'wages':wages, 'wageu':wageu,
                       'Ls':Ls, 'Lu':Lu, 'Lso':Lso, 'Luo':Luo, 
                       'manualWs':manualWs, 'manualWu':manualWu, 
                       'normAs':normAs, 'normAu':normAu, 
                       'normTs':normTs, 'normTu':normTu,  
                       'normbarAs':normbarAs, 'normbarAu':normbarAu, 
                       'normbarTs':normbarTs, 'normbarTu':normbarTu, 
                       'wages_pred':wages_pred, 'wageu_pred':wageu_pred, 
                       'Ls_pred':Ls_pred, 'Lu_pred':Lu_pred,  
                       'Y':Y, 'Ys':Ys, 'Yu':Yu, 'xi':xi})
export.to_csv(resultsCSV)



# ---------------------------------------------------------------------------
#  Exporting unnormalized version
# ---------------------------------------------------------------------------   
# selected unnormalized metrics you chose to include
u_As=u_As[:,0]
u_barAs=u_barAs[:,0]
u_Au=u_Au[:,0]
u_barAu=u_barAu[:,0]
u_Ts=u_Ts[:,0]
u_barTs=u_barTs[:,0]
u_Tu=u_Tu[:,0]
u_barTu=u_barTu[:,0]
u_ws=u_ws[:,0]
u_Ws=u_Ws[:,0]
u_Wu=u_Wu[:,0]
u_wu=u_wu[:,0]
u_P=u_P[:,0]
u_p=u_p[:,0]
u_wages=u_wages[:,0]
u_wageu=u_wageu[:,0]
u_Lso=u_Lso[:,0] 
u_Luo=u_Luo[:,0] 
u_Ls=u_Ls[:,0] 
u_Lu=u_Lu[:,0] 
u_manualWs=u_manualWs[:,0]
u_manualWu=u_manualWu[:,0]


export = pd.DataFrame({'id': id_2, 'u_ws': u_ws, 'u_wu': u_wu,  
                       'u_Ws': u_Ws, 'u_Wu': u_Wu, 
                       'u_P': u_P, 'u_p': u_p, 
                       'u_As':u_As, 'u_Au':u_Au, 
                       'u_barAs':u_barAs, 'u_barAu':u_barAu, 
                       'u_Ts':u_Ts, 'u_Tu':u_Tu,
                       'u_barTs':u_barTs, 'u_barTu':u_barTu, 
                       'u_wages':u_wages, 'u_wageu':u_wageu,
                       'u_Ls':u_Ls, 'u_Lu':u_Lu, 'u_Lso':u_Lso, 'u_Luo':u_Luo, 
                       'u_manualWs':u_manualWs, 'u_manualWu':u_manualWu })
export.to_csv(resultsCSV2)



