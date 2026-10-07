import numpy as np
import pandas as pd

# ===========================================================================
# Your climate elasticity estimates from the table
# ===========================================================================
eta_As = -0.0348  # Amenity, Skilled
se_As = 0.0114

eta_Au = -0.0030   # Amenity, Low-Skilled
se_Au = 0.0091

eta_Ts = -0.0436   # Productivity, Skilled
se_Ts = 0.0201

eta_Tu = -0.0930   # Productivity, Low-Skilled
se_Tu = 0.0369

print("Climate Elasticity Parameters:")
print("="*70)
print(f"Amenity - Skilled:      {eta_As:7.4f} (SE: {se_As:.4f})")
print(f"Amenity - Low-Skilled:  {eta_Au:7.4f} (SE: {se_Au:.4f})")
print(f"Productivity - Skilled: {eta_Ts:7.4f} (SE: {se_Ts:.4f})")
print(f"Productivity - Low-Skilled: {eta_Tu:7.4f} (SE: {se_Tu:.4f})")
print()

# Check criteria
print("Checking criteria:")
print(f"✓ Both productivity negative: Ts={eta_Ts:.4f} < 0, Tu={eta_Tu:.4f} < 0")
print(f"✓ |Tu| > |Ts|: {abs(eta_Tu):.4f} > {abs(eta_Ts):.4f}")
print(f"✓ Amenity sum: {eta_As + eta_Au:.4f} (NEGATIVE!)")
print()

# ===========================================================================
# Generate bootstrap draws
# ===========================================================================
B = 300  # Number of bootstrap iterations
np.random.seed(42)

print(f"Generating {B} bootstrap draws...")

# Draw from sampling distributions (independent draws)
draws_As = np.random.normal(eta_As, se_As, B)
draws_Au = np.random.normal(eta_Au, se_Au, B)
draws_Ts = np.random.normal(eta_Ts, se_Ts, B)
draws_Tu = np.random.normal(eta_Tu, se_Tu, B)

# Create DataFrame
bootstrap_draws = pd.DataFrame({
    'draw': range(B),
    'eta_As': draws_As,
    'eta_Au': draws_Au,
    'eta_Ts': draws_Ts,
    'eta_Tu': draws_Tu
})

# Save to CSV
bootstrap_draws.to_csv('bootstrap_draws.csv', index=False)
print(f"✓ Saved {B} parameter draws to 'bootstrap_draws.csv'")
print()

# ===========================================================================
# Summary statistics
# ===========================================================================
print("Distribution of bootstrap draws:")
print("="*70)
print(f"{'Parameter':<12} {'Mean':<10} {'2.5%':<10} {'97.5%':<10}")
print("-"*70)

for param_name, draws in [('eta_As', draws_As), ('eta_Au', draws_Au), 
                          ('eta_Ts', draws_Ts), ('eta_Tu', draws_Tu)]:
    ci = np.percentile(draws, [2.5, 97.5])
    print(f"{param_name:<12} {draws.mean():>9.4f} {ci[0]:>9.4f} {ci[1]:>9.4f}")

print()
print("Ready to run bootstrap counterfactuals!")
print(f"Next step: Run cf_main_bootstrap.py")
