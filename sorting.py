import pandas as pd
import numpy as np

# Load the CSV
df = pd.read_csv('boltz/list_with_affinities_boltz.csv')

# Select columns to keep
columns_to_keep = ['smiles', 'id-num', 'vina_affinity', 'boltz_affinity_kcalmol', 
                   'avg_affinity_pred_value', 'avg_affinity_probability_binary']
df_filtered = df[columns_to_keep].copy()

print(f"Total molecules: {len(df)}")
print(f"Molecules with Vina affinity: {df['vina_affinity'].notna().sum()}")
print(f"Molecules with Boltz affinity: {df['boltz_affinity_kcalmol'].notna().sum()}")
print(f"Molecules with both affinities: {(df['vina_affinity'].notna() & df['boltz_affinity_kcalmol'].notna()).sum()}")

# Filter molecules that have both affinities
df_both = df_filtered[(df_filtered['vina_affinity'].notna()) & 
                      (df_filtered['boltz_affinity_kcalmol'].notna())].copy()

# Calculate combined score by averaging
if len(df_both) > 0:
    vina_min = df_both['vina_affinity'].min()
    vina_max = df_both['vina_affinity'].max()
    boltz_min = df_both['boltz_affinity_kcalmol'].min()
    boltz_max = df_both['boltz_affinity_kcalmol'].max()
    
    # For Vina: more negative is better, so we invert it
    # For Boltz: higher is better, so we keep it positive
    # Average them: (vina_affinity + (-boltz_affinity)) / 2
    # This means: lower (more negative) combined score is better
    df_both['combined_score'] = (df_both['vina_affinity'] - df_both['boltz_affinity_kcalmol']) / 2
    
    # Sort by combined score (lowest/most negative = best)
    df_both_sorted = df_both.sort_values('combined_score')
    
    # Get molecules without both affinities
    df_missing = df_filtered[~((df_filtered['vina_affinity'].notna()) & 
                               (df_filtered['boltz_affinity_kcalmol'].notna()))].copy()
    
    # Combine sorted and missing
    df_sorted = pd.concat([df_both_sorted, df_missing], ignore_index=True)
else:
    df_sorted = df_filtered.copy()

# Get top 10
df_best10 = df_sorted.head(10)[columns_to_keep].copy()

# Save results
df_sorted.to_csv('list-sorted.csv', index=False)
df_best10.to_csv('list-best10.csv', index=False)

print(f"\n✓ Saved sorted list with {len(df_sorted)} molecules to 'list-sorted.csv'")
print(f"✓ Saved top 10 molecules to 'list-best10.csv'")

# Print data ranges
if len(df_both) > 0:
    print("\n=== DATA RANGES ===")
    print(f"Vina affinity range: {vina_min:.3f} to {vina_max:.3f}")
    print(f"Boltz affinity range: {boltz_min:.3f} to {boltz_max:.3f}")

print("\n=== SORTING FORMULA EXPLANATION ===")
print("SIMPLE AVERAGING METHOD:\n")
print("1. Vina affinity: More negative = better binding (e.g., -9.5 is better than -8.0)")
print("2. Boltz affinity: Higher values = better binding (e.g., 8.0 is better than 6.5)")
print()
print("3. Combined score formula:")
print("   combined_score = (vina_affinity - boltz_affinity_kcalmol) / 2")
print()
print("   Why subtract Boltz?")
print("   - We want to reward LOW vina (e.g., -9.5)")
print("   - We want to reward HIGH boltz (e.g., 8.0)")
print("   - Subtracting high boltz makes the score MORE negative (better)")
print()
print("   Example:")
print("   - Molecule A: Vina=-9.5, Boltz=8.0 → Score=(-9.5-8.0)/2 = -8.75")
print("   - Molecule B: Vina=-8.0, Boltz=6.5 → Score=(-8.0-6.5)/2 = -7.25")
print("   - Molecule A wins (more negative score)")
print()
print("4. Molecules are sorted by combined score (most negative = best)")
print("5. Molecules missing either value are placed at the end")

print("\n=== TOP 10 MOLECULES ===")
pd.set_option('display.max_colwidth', None)
pd.set_option('display.width', None)
print(df_best10[columns_to_keep].to_string(index=True))

# Show statistics of the top 10
if len(df_both) > 0:
    print("\n=== STATISTICS OF TOP 10 ===")
    molecules_with_both = df_both_sorted.head(10)
    print(f"Total molecules with both affinities: {len(df_both)}")
    print(f"Best combined score: {molecules_with_both['combined_score'].iloc[0]:.3f}")
    if len(molecules_with_both) >= 10:
        print(f"10th best combined score: {molecules_with_both['combined_score'].iloc[9]:.3f}")
    print(f"Average Vina in top 10: {molecules_with_both['vina_affinity'].mean():.3f}")
    print(f"Average Boltz in top 10: {molecules_with_both['boltz_affinity_kcalmol'].mean():.3f}")
    print(f"Average pred_value in top 10: {molecules_with_both['avg_affinity_pred_value'].mean():.3f}")
    print(f"Average probability in top 10: {molecules_with_both['avg_affinity_probability_binary'].mean():.3f}")
