import numpy as np
import pandas as pd
from scipy.optimize import least_squares
from root import ROOT

##  Set parameters
sigma = 4   # goods substitution
theta = 1.5 # frechet parameter
kappa = 1.5 # skill substitutability
eta = -.05
mus = .02
muu = .02
N = 1600
N_foreign = 1  # ROW
N_total = N + N_foreign
vgamma=(theta-1)/theta
I = np.eye(N)
I_full = np.eye(N_total)

beta_hat = - 1.199
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

def norm(data):
    return np.divide(data,np.mean(data))

#SET INVERSION PROCEDURE YEAR
year = 2000  # RUN 2010 FIRST!
migsimple =0

# ===========================================================================
#  Step #1: Loading data
# ===========================================================================
CSVFILES = ROOT + 'inputs_csv/'
OUTPUTS = ROOT + 'outputs_csv/mig_full/'

# Load FULL matrices (including foreign row/column 8887)
data_dist_full = np.array(pd.read_csv(CSVFILES + 'distance.csv'))[:, 1:]
diff_long_full = np.array(pd.read_csv(CSVFILES + 'diff_long.csv'))[:, 1:]
diff_lat_full = np.array(pd.read_csv(CSVFILES + 'diff_lat.csv'))[:, 1:]
sameprov_full = np.array(pd.read_csv(CSVFILES + 'sameprov.csv'))[:, 1:]
sameisland_full = np.array(pd.read_csv(CSVFILES + 'sameisland.csv'))[:, 1:]

# Domestic only for trade
data_dist = data_dist_full[:N, :N]
diff_long = diff_long_full[:N, :N]
diff_lat = diff_lat_full[:N, :N]
sameprov = sameprov_full[:N, :N]
sameisland = sameisland_full[:N, :N]

# Load data (keep row 8887 for migration)
data_L2010_full = pd.read_csv(CSVFILES + 'Ld_2010.csv').to_numpy()
data_Lo2010_full = pd.read_csv(CSVFILES + 'Lo_2010.csv').to_numpy()
data_L2000_full = pd.read_csv(CSVFILES + 'Ld_2000.csv').to_numpy()
data_Lo2000_full = pd.read_csv(CSVFILES + 'Lo_2000.csv').to_numpy()

# Domestic only (filter out row 8887)
data_Y2010 = pd.read_csv(CSVFILES + 'localGDP.csv')
data_Y2010 = data_Y2010[data_Y2010.iloc[:, 0] != 8887].to_numpy()

data_L2010 = data_L2010_full[data_L2010_full[:, 0] != 8887]
data_Lo2010 = data_Lo2010_full[data_Lo2010_full[:, 0] != 8887]
data_wage2010 = pd.read_csv(CSVFILES + 'wagesresidSector_2010.csv')
data_wage2010 = data_wage2010[data_wage2010.iloc[:, 0] != 8887].to_numpy()

data_L2000 = data_L2000_full[data_L2000_full[:, 0] != 8887]
data_Lo2000 = data_Lo2000_full[data_Lo2000_full[:, 0] != 8887]
data_wage2000 = pd.read_csv(CSVFILES + 'wagesresidSector_2000.csv')
data_wage2000 = data_wage2000[data_wage2000.iloc[:, 0] != 8887].to_numpy()

data_slr = np.array(pd.read_csv(CSVFILES + 'slr_area.csv'))
area = data_slr[:,[1]]
GDP = data_Y2010[:,[4]]
popcell = data_Y2010[:,[5]]

if year==2010:
    resultsCSV= OUTPUTS+'resultspy_new_2010_ROW.csv'
    resultsCSV2= OUTPUTS+'u_resultspy_new_2010_ROW.csv'
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
    resultsCSV= OUTPUTS+'resultspy_new_2000_ROW.csv'
    resultsCSV2= OUTPUTS+'u_resultspy_new_2000_ROW.csv'
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
Lsshare=np.divide(np.sum(Ls),np.sum(L))
Lushare=np.divide(np.sum(Lu),np.sum(L))

GDP = ((GDP*1000000000)/popcell)*L
GDP = norm(GDP)

Y = (np.multiply(Ls,wages)) + (np.multiply(Lu,wageu))
Ys = (np.multiply(Ls,wages))
Yu = (np.multiply(Lu,wageu))
dens = np.divide(L,area)

dens=norm(dens)
Y=norm(Y)
Ys = norm(Ys)
Yu = norm(Yu)

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
#  Step #2.1: Distance costs
# ===========================================================================
beta_hat = beta_hat / (sigma - 1)
betas_hat, betau_hat, provs, provu, lats, latu, longs, longu, homeu, homes, islands, islandu = [x / theta for x in [betas_hat, betau_hat, provs, provu, lats, latu, longs, longu, homeu, homes, islands, islandu]]

# FIND DOMESTIC VS FOREIGN INDICES FIRST
domestic_indices = np.where(data_Lo2010_full[:, 0] != 8887)[0]
foreign_indices = np.where(data_Lo2010_full[:, 0] == 8887)[0]

# TRADE MATRICES - domestic only, properly filtered
data_dist = data_dist_full[np.ix_(domestic_indices, domestic_indices)]
diff_long = diff_long_full[np.ix_(domestic_indices, domestic_indices)]
diff_lat = diff_lat_full[np.ix_(domestic_indices, domestic_indices)]
sameprov = sameprov_full[np.ix_(domestic_indices, domestic_indices)]
sameisland = sameisland_full[np.ix_(domestic_indices, domestic_indices)]

T = np.exp(beta_hat * np.arcsinh(data_dist))
T = T / T.diagonal()[:, None]

# MIGRATION MATRICES - FULL (including foreign)
Ms_full = (betas_hat * np.arcsinh(data_dist_full)
    + islands * sameisland_full
    + provs * sameprov_full
    + lats * diff_lat_full
    + longs * diff_long_full
    + homes * I_full
)

Mu_full = (betau_hat * np.arcsinh(data_dist_full)
    + islandu * sameisland_full
    + provu * sameprov_full
    + latu * diff_lat_full
    + longu * diff_long_full
    + homeu * I_full
)

Mu_full = np.exp(Mu_full)
Ms_full = np.exp(Ms_full)
Ms_full = Ms_full / Ms_full.diagonal()[:, None] 
Mu_full = Mu_full / Mu_full.diagonal()[:, None]

# Extract domestic and foreign parts
Ms = Ms_full[np.ix_(domestic_indices, domestic_indices)]
Mu = Mu_full[np.ix_(domestic_indices, domestic_indices)]
Ms_to_foreign = Ms_full[np.ix_(domestic_indices, foreign_indices)]
Mu_to_foreign = Mu_full[np.ix_(domestic_indices, foreign_indices)]

# ===========================================================================
#  Step #2.2: Invert for Omega
# ===========================================================================
print("\n" + "="*60)
print("INVERTING FOR FOREIGN ATTRACTIVENESS (OMEGA)")
print("="*60)

foreign_row = data_Lo2010_full[data_Lo2010_full[:, 0] == 8887]
foreign_Ls = foreign_row[0, 2]
foreign_Lu = foreign_row[0, 3]
domestic_Lso = np.sum(data_Lo2010_full[data_Lo2010_full[:, 0] != 8887, 2])
domestic_Luo = np.sum(data_Lo2010_full[data_Lo2010_full[:, 0] != 8887, 3])

pi_foreign_s = foreign_Ls / domestic_Lso
pi_foreign_u = foreign_Lu / domestic_Luo

print(f"Emigration rate (skilled): {100*pi_foreign_s:.2f}%")
print(f"Emigration rate (unskilled): {100*pi_foreign_u:.2f}%")

avg_mu_inv_s = np.mean(Ms_to_foreign)
avg_mu_inv_u = np.mean(Mu_to_foreign)

ratio_s = pi_foreign_s / (1 - pi_foreign_s)
ratio_u = pi_foreign_u / (1 - pi_foreign_u)

Omega_s = (ratio_s / avg_mu_inv_s) ** (1/theta)
Omega_u = (ratio_u / avg_mu_inv_u) ** (1/theta)

print(f"\nInverted Omega_s: {Omega_s:.4f}")
print(f"Inverted Omega_u: {Omega_u:.4f}")
print("="*60 + "\n")

# ===========================================================================
#  Step #3: Setting-up guesses
# ===========================================================================
Ws = np.ones((N,1))
ws = np.ones((N,1))
Wu = np.ones((N,1))
wu = np.ones((N,1))
P = np.ones((N,1))
p = np.ones((N,1))
diff1 = 1
tol = 1e-09
alpha_update = 1
iter=1

domestic_part = Ms @ (1/ws)
foreign_part = (Omega_s ** theta) * Ms_to_foreign
total = domestic_part + foreign_part


# ===========================================================================
#  Welfare and Prices
# ===========================================================================
while diff1 > tol:
    print(iter, '--', diff1)
    
    # ADD FOREIGN TERM HERE - ONLY CHANGE TO ORIGINAL CODE
    Ws_new = Ms @ (1 / ws) + (Omega_s ** theta) * Ms_to_foreign
    Wu_new = Mu @ (1 / wu) + (Omega_u ** theta) * Mu_to_foreign
    
    # Everything else EXACTLY THE SAME as original
    ws_new = np.sum([np.multiply(Ms, np.matmul(1/Ls, np.transpose(np.multiply(Lso,np.power(Ws,-1)))))], axis=2)   
    ws_new = ws_new.T

    wu_new = np.sum([np.multiply(Mu, np.matmul(1/Lu, np.transpose(np.multiply(Luo,np.power(Wu,-1)))))], axis=2)
    wu_new = wu_new.T

    u_ws = ws_new
    u_wu = wu_new
    
    ws_new = norm(ws_new)
    wu_new = norm(wu_new)
    
    
    p_new = np.sum([
        np.multiply(T, np.matmul(1/Y, np.transpose(Y * (P**(sigma-1)))))], axis=2)
    p_new = (p_new.T) 
    
    u_p = p_new
    p_new = norm(p_new)
    P_new = T @ (p_new**(-1))
    
    
    print(f"\nDEBUG iteration {iter}:")
    print(f"  ws_new: min={ws_new.min():.6f}, max={ws_new.max():.6f}, has NaN={np.isnan(ws_new).any()}")
    print(f"  wu_new: min={wu_new.min():.6f}, max={wu_new.max():.6f}, has NaN={np.isnan(wu_new).any()}")
    print(f"  p_new: min={p_new.min():.6f}, max={p_new.max():.6f}, has NaN={np.isnan(p_new).any()}")
    print(f"  P_new: min={P_new.min():.6f}, max={P_new.max():.6f}, has NaN={np.isnan(P_new).any()}")
    
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
    
    iter = iter + 1       
    ws = alpha_update * ws_new + (1 - alpha_update) * ws
    Ws = alpha_update * Ws_new + (1 - alpha_update) * Ws
    wu = alpha_update * wu_new + (1 - alpha_update) * wu
    Wu = alpha_update * Wu_new + (1 - alpha_update) * Wu
    p = alpha_update * p_new + (1 - alpha_update) * p
    P = alpha_update * P_new + (1 - alpha_update) * P

# Recover fundamentals
Ws = Ws**(1/theta)   
Wu = Wu**(1/theta)   
ws = ws**(-1/theta)   
wu = wu**(-1/theta) 
p = p**(1/(sigma-1))
P = P**(1/(1 - sigma))

As = (P*ws) / wages
Au = (P*wu) / wageu
barAs = As / (dens**eta) 
barAu = Au / (dens**eta)

# Productivities
if year == 2000:
    data_2010 = np.array(pd.read_csv(OUTPUTS + 'resultspy_new_2010_ROW.csv'))
    xi = data_2010[:,[39]]
    Ts = data_2010[:,[12]] 
    Tu = data_2010[:,[13]] 
    Ts = Ts / tfp_scale_factor
    Tu = Tu / tfp_scale_factor

if year == 2010:
    def solve_production_params_bounded(wages, wageu, Ls, Lu, GDP, p, kappa):
        def equations(params):
            xi, Ts, Tu = params
            eq1 = wages - p * xi * (Ts**((kappa-1)/kappa)) * (Ls**(-1/kappa)) * (GDP**(1/kappa))
            eq2 = wageu - p * (1-xi) * (Tu**((kappa-1)/kappa)) * (Lu**(-1/kappa)) * (GDP**(1/kappa))
            term_H = xi * (Ts * Ls)**((kappa-1)/kappa)
            term_L = (1-xi) * (Tu * Lu)**((kappa-1)/kappa)
            GDP_pred = (term_H + term_L)**(kappa/(kappa-1))
            eq3 = GDP - GDP_pred
            return [eq1, eq2, eq3]
        
        x0 = [0.5, 1.0, 1.0]
        bounds = ([0.01, 0.001, 0.001], [0.99, 100, 100])
        result = least_squares(equations, x0, bounds=bounds)
        return result.x, result.success
    
    results = []
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
                results.append({
                    'municipality': idx,
                    'xi': np.nan,
                    'Ts': np.nan,
                    'Tu': np.nan,
                    'converged': False
                })
        except Exception as e:
            results.append({
                'municipality': idx,
                'xi': np.nan,
                'Ts': np.nan,
                'Tu': np.nan,
                'converged': False
            })
    
    params_df = pd.DataFrame(results)
    xi = params_df['xi'].values.reshape(-1, 1)
    Ts = params_df['Ts'].values.reshape(-1, 1)
    Tu = params_df['Tu'].values.reshape(-1, 1)

barTs = Ts / (dens**mus)
barTu = Tu / (dens**muu)

normAs=norm(As)
normTs=norm(Ts)
normTu=norm(Tu)
normAu=norm(Au)
normbarAs=norm(barAs)
normbarTs=norm(barTs)
normbarTu=norm(barTu)
normbarAu=norm(barAu)

# Model predictions
wages_cfx= (P*ws) / As
wageu_cfx= (P*wu) / Au

Ls_cfx = np.multiply(np.power(ws,theta),np.matmul(Ms, np.divide(Lso,np.power(Ws,theta))))
Lu_cfx = np.multiply(np.power(wu,theta),np.matmul(Mu, np.divide(Luo,np.power(Wu,theta))))



manualWs = (Ms @ ((As * wages_cfx / P)**theta))**(1/theta)
manualWu = (Mu @ ((Au * wageu_cfx / P)**theta))**(1/theta)



# Unnormalized
u_Ws = np.matmul(Ms, np.power(u_ws,theta))
u_Ws = np.power(u_Ws, 1/theta)
u_Wu = np.matmul(Mu, np.power(u_wu,theta))
u_Wu = np.power(u_Wu, 1/theta)

u_P = np.matmul(T, np.power(u_p,1-sigma))
u_P = np.power(u_P,(1/(1-sigma)))

u_Lso = data_Lo2010[:,[2]]
u_Luo = data_Lo2010[:,[3]] 
u_Ls = data_L2010[:,[2]]
u_Lu = data_L2010[:,[3]]

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

u_manualWs = np.power(np.matmul(Ms,np.power(np.divide(np.multiply(As,u_wages),u_P),theta)),1/theta)
u_manualWu = np.power(np.matmul(Mu,np.power(np.divide(np.multiply(Au,u_wageu),u_P),theta)),1/theta)



# Flatten for export
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
                       'Y':Y, 'Ys':Ys, 'Yu':Yu, 'xi':xi,
                       'Omega_s':Omega_s, 'Omega_u':Omega_u})
export.to_csv(resultsCSV)

# Unnormalized
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


print(f"\n{'='*60}")
print("INVERSION COMPLETE - ROW VERSION")
print(f"Omega_s: {Omega_s:.4f}, Omega_u: {Omega_u:.4f}")
print(f"{'='*60}\n")