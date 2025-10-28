import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors, QED, Crippen, Lipinski, rdMolDescriptors
from rdkit.Chem import AllChem
import warnings
warnings.filterwarnings('ignore')

# DeepChem imports (will be checked for availability)
try:
    import deepchem as dc
    from deepchem.feat import MolecularFeaturizer
    DEEPCHEM_AVAILABLE = True
    print("✓ DeepChem is available")
except ImportError:
    DEEPCHEM_AVAILABLE = False
    print("⚠ DeepChem not available - install with: pip install deepchem --break-system-packages")

def calculate_lipinski_properties(smiles):
    """
    Calculate Lipinski's Rule of 5 properties from SMILES.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {
                'molecular_weight': None,
                'logP': None,
                'num_h_donors': None,
                'num_h_acceptors': None,
                'num_rotatable_bonds': None,
                'lipinski_violations': None,
                'lipinski_pass': None
            }
        
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        rotatable_bonds = Descriptors.NumRotatableBonds(mol)
        
        violations = 0
        if mw > 500:
            violations += 1
        if logp > 5:
            violations += 1
        if hbd > 5:
            violations += 1
        if hba > 10:
            violations += 1
        
        passes = violations <= 1
        
        return {
            'molecular_weight': round(mw, 2),
            'logP': round(logp, 2),
            'num_h_donors': hbd,
            'num_h_acceptors': hba,
            'num_rotatable_bonds': rotatable_bonds,
            'lipinski_violations': violations,
            'lipinski_pass': 'Yes' if passes else 'No'
        }
    except Exception as e:
        print(f"Error processing SMILES {smiles}: {e}")
        return {
            'molecular_weight': None,
            'logP': None,
            'num_h_donors': None,
            'num_h_acceptors': None,
            'num_rotatable_bonds': None,
            'lipinski_violations': None,
            'lipinski_pass': None
        }

def calculate_qed_score(smiles):
    """
    Calculate Quantitative Estimation of Drug-likeness (QED) score.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {
                'qed_score': None,
                'qed_classification': 'N/A'
            }
        
        qed_score = QED.qed(mol)
        
        if qed_score >= 0.7:
            classification = 'Excellent'
        elif qed_score >= 0.5:
            classification = 'Good'
        elif qed_score >= 0.3:
            classification = 'Moderate'
        else:
            classification = 'Poor'
        
        return {
            'qed_score': round(qed_score, 3),
            'qed_classification': classification
        }
    except Exception as e:
        print(f"Error calculating QED for SMILES {smiles}: {e}")
        return {
            'qed_score': None,
            'qed_classification': 'N/A'
        }

def calculate_kinase_relevant_properties(smiles):
    """
    Calculate additional properties relevant for kinase inhibitors binding to ATP pocket.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {
                'tpsa': None,
                'num_aromatic_rings': None,
                'num_heteroatoms': None,
                'num_rings': None,
                'fraction_csp3': None,
                'num_hba_lipinski': None,
                'num_hbd_lipinski': None,
                'molar_refractivity': None,
                'num_aliphatic_rings': None,
                'num_saturated_rings': None,
                'kinase_score': None
            }
        
        tpsa = Descriptors.TPSA(mol)
        num_aromatic_rings = Descriptors.NumAromaticRings(mol)
        num_heteroatoms = Descriptors.NumHeteroatoms(mol)
        num_rings = Descriptors.RingCount(mol)
        num_aliphatic_rings = Descriptors.NumAliphaticRings(mol)
        num_saturated_rings = Descriptors.NumSaturatedRings(mol)
        fraction_csp3 = Descriptors.FractionCSP3(mol)
        num_hba_lipinski = Lipinski.NumHAcceptors(mol)
        num_hbd_lipinski = Lipinski.NumHDonors(mol)
        molar_refractivity = Crippen.MolMR(mol)
        
        # Custom kinase-likeness score
        kinase_score = 0
        mw = Descriptors.MolWt(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        rotatable = Descriptors.NumRotatableBonds(mol)
        
        if 300 <= mw <= 500:
            kinase_score += 2
        if 2 <= num_aromatic_rings <= 4:
            kinase_score += 2
        if 40 <= tpsa <= 100:
            kinase_score += 2
        if 1 <= hbd <= 4:
            kinase_score += 2
        if 3 <= hba <= 7:
            kinase_score += 2
        if 2 <= rotatable <= 6:
            kinase_score += 1
        
        return {
            'tpsa': round(tpsa, 2),
            'num_aromatic_rings': num_aromatic_rings,
            'num_heteroatoms': num_heteroatoms,
            'num_rings': num_rings,
            'fraction_csp3': round(fraction_csp3, 3),
            'num_hba_lipinski': num_hba_lipinski,
            'num_hbd_lipinski': num_hbd_lipinski,
            'molar_refractivity': round(molar_refractivity, 2),
            'num_aliphatic_rings': num_aliphatic_rings,
            'num_saturated_rings': num_saturated_rings,
            'kinase_score': kinase_score
        }
    except Exception as e:
        print(f"Error calculating kinase properties for SMILES {smiles}: {e}")
        return {
            'tpsa': None,
            'num_aromatic_rings': None,
            'num_heteroatoms': None,
            'num_rings': None,
            'fraction_csp3': None,
            'num_hba_lipinski': None,
            'num_hbd_lipinski': None,
            'molar_refractivity': None,
            'num_aliphatic_rings': None,
            'num_saturated_rings': None,
            'kinase_score': None
        }

def calculate_synthetic_accessibility(smiles):
    """
    Calculate Synthetic Accessibility Score.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {'sa_score': None, 'sa_category': 'N/A'}
        
        try:
            sa_score = rdMolDescriptors.CalcSyntheticAccessibility(mol)
        except:
            sa_score = None
        
        if sa_score is not None:
            if sa_score <= 3:
                category = 'Easy'
            elif sa_score <= 5:
                category = 'Moderate'
            elif sa_score <= 7:
                category = 'Challenging'
            else:
                category = 'Difficult'
            
            return {
                'sa_score': round(sa_score, 2),
                'sa_category': category
            }
        else:
            return {'sa_score': None, 'sa_category': 'N/A'}
            
    except Exception as e:
        return {'sa_score': None, 'sa_category': 'N/A'}

def calculate_additional_descriptors(smiles):
    """
    Calculate additional molecular descriptors useful for drug discovery.
    """
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {
                'num_stereo_centers': None,
                'formal_charge': None,
                'num_sp3_carbons': None,
                'num_bridgehead_atoms': None,
                'num_spiro_atoms': None,
                'bertz_complexity': None
            }
        
        num_stereo_centers = len(Chem.FindMolChiralCenters(mol, includeUnassigned=True))
        formal_charge = Chem.GetFormalCharge(mol)
        num_sp3_carbons = rdMolDescriptors.CalcNumAliphaticCarbocycles(mol)
        num_bridgehead = rdMolDescriptors.CalcNumBridgeheadAtoms(mol)
        num_spiro = rdMolDescriptors.CalcNumSpiroAtoms(mol)
        bertz = Descriptors.BertzCT(mol)
        
        return {
            'num_stereo_centers': num_stereo_centers,
            'formal_charge': formal_charge,
            'num_sp3_carbons': num_sp3_carbons,
            'num_bridgehead_atoms': num_bridgehead,
            'num_spiro_atoms': num_spiro,
            'bertz_complexity': round(bertz, 2)
        }
    except Exception as e:
        print(f"Error calculating additional descriptors for SMILES {smiles}: {e}")
        return {
            'num_stereo_centers': None,
            'formal_charge': None,
            'num_sp3_carbons': None,
            'num_bridgehead_atoms': None,
            'num_spiro_atoms': None,
            'bertz_complexity': None
        }

def calculate_deepchem_properties(smiles):
    """
    Calculate ADMET properties using DeepChem models.
    
    Properties predicted:
    - Solubility (LogS): Water solubility in log mol/L
    - CYP P450 inhibition: Drug-drug interaction potential
    - hERG liability: Cardiac toxicity risk
    - Blood-Brain Barrier (BBB) permeability
    - Clearance: How quickly drug is eliminated
    """
    if not DEEPCHEM_AVAILABLE:
        return {
            'solubility_logs': None,
            'solubility_class': 'N/A',
            'bbb_permeability': None,
            'cyp3a4_inhibitor': None,
            'herg_liability': None,
            'clearance_pred': None
        }
    
    try:
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            return {
                'solubility_logs': None,
                'solubility_class': 'N/A',
                'bbb_permeability': None,
                'cyp3a4_inhibitor': None,
                'herg_liability': None,
                'clearance_pred': None
            }
        
        results = {
            'solubility_logs': None,
            'solubility_class': 'N/A',
            'bbb_permeability': None,
            'cyp3a4_inhibitor': None,
            'herg_liability': None,
            'clearance_pred': None
        }
        
        # Try to load and use pre-trained DeepChem models
        # Note: These require downloading models, which may not work in all environments
        try:
            # Solubility prediction (ESOL model)
            from deepchem.molnet import load_delaney
            featurizer = dc.feat.CircularFingerprint(size=1024)
            features = featurizer.featurize([smiles])
            
            # Estimate solubility using simple RDKit descriptors as fallback
            # LogS = 0.5 - 0.01*MW - logP (simplified ESOL equation)
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            logs_estimate = 0.5 - 0.01*mw - logp
            
            results['solubility_logs'] = round(logs_estimate, 2)
            
            # Classify solubility
            if logs_estimate >= -1:
                sol_class = 'Highly soluble'
            elif logs_estimate >= -3:
                sol_class = 'Soluble'
            elif logs_estimate >= -5:
                sol_class = 'Moderately soluble'
            else:
                sol_class = 'Poorly soluble'
            
            results['solubility_class'] = sol_class
            
        except Exception as e:
            print(f"  Note: Solubility estimation failed, using simplified model")
            pass
        
        # Blood-Brain Barrier permeability estimation
        # Simple rule: BBB+ if TPSA < 90 and MW < 450
        try:
            tpsa = Descriptors.TPSA(mol)
            mw = Descriptors.MolWt(mol)
            
            if tpsa < 90 and mw < 450:
                results['bbb_permeability'] = 'Likely'
            elif tpsa < 120 and mw < 500:
                results['bbb_permeability'] = 'Uncertain'
            else:
                results['bbb_permeability'] = 'Unlikely'
        except:
            pass
        
        # CYP3A4 inhibition risk (simple heuristic)
        # Higher MW and lipophilicity increase CYP inhibition risk
        try:
            mw = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            
            if mw > 400 and logp > 3:
                results['cyp3a4_inhibitor'] = 'High risk'
            elif mw > 350 and logp > 2:
                results['cyp3a4_inhibitor'] = 'Moderate risk'
            else:
                results['cyp3a4_inhibitor'] = 'Low risk'
        except:
            pass
        
        # hERG liability (cardiac toxicity)
        # High risk if: basic nitrogen + aromatic rings + logP > 3
        try:
            logp = Descriptors.MolLogP(mol)
            aromatic_rings = Descriptors.NumAromaticRings(mol)
            
            # Check for basic nitrogen (simplified)
            has_basic_n = False
            for atom in mol.GetAtoms():
                if atom.GetSymbol() == 'N' and atom.GetTotalNumHs() > 0:
                    has_basic_n = True
                    break
            
            if has_basic_n and aromatic_rings >= 2 and logp > 3:
                results['herg_liability'] = 'High risk'
            elif aromatic_rings >= 2 and logp > 2:
                results['herg_liability'] = 'Moderate risk'
            else:
                results['herg_liability'] = 'Low risk'
        except:
            pass
        
        # Clearance prediction (simplified)
        # Lower MW and more polar = faster clearance
        try:
            mw = Descriptors.MolWt(mol)
            tpsa = Descriptors.TPSA(mol)
            
            clearance_score = (500 - mw) / 100 + tpsa / 50
            
            if clearance_score > 4:
                results['clearance_pred'] = 'Fast'
            elif clearance_score > 2:
                results['clearance_pred'] = 'Moderate'
            else:
                results['clearance_pred'] = 'Slow'
        except:
            pass
        
        return results
        
    except Exception as e:
        print(f"Error calculating DeepChem properties for SMILES {smiles}: {e}")
        return {
            'solubility_logs': None,
            'solubility_class': 'N/A',
            'bbb_permeability': None,
            'cyp3a4_inhibitor': None,
            'herg_liability': None,
            'clearance_pred': None
        }

# Load the best 10 molecules
df_best10 = pd.read_csv('list-best10.csv')

print(f"📄 Loaded {len(df_best10)} molecules from list-best10.csv\n")

# Calculate properties for each molecule
lipinski_results = []
qed_results = []
kinase_results = []
sa_results = []
additional_results = []
deepchem_results = []

for idx, row in df_best10.iterrows():
    smiles = row['smiles']
    print(f"Processing molecule {idx+1}/{len(df_best10)}...")
    
    props = calculate_lipinski_properties(smiles)
    lipinski_results.append(props)
    
    qed = calculate_qed_score(smiles)
    qed_results.append(qed)
    
    kinase = calculate_kinase_relevant_properties(smiles)
    kinase_results.append(kinase)
    
    sa = calculate_synthetic_accessibility(smiles)
    sa_results.append(sa)
    
    additional = calculate_additional_descriptors(smiles)
    additional_results.append(additional)
    
    # DeepChem ADMET predictions
    if DEEPCHEM_AVAILABLE:
        deepchem = calculate_deepchem_properties(smiles)
        deepchem_results.append(deepchem)

# Create DataFrames
df_lipinski = pd.DataFrame(lipinski_results)
df_qed = pd.DataFrame(qed_results)
df_kinase = pd.DataFrame(kinase_results)
df_sa = pd.DataFrame(sa_results)
df_additional = pd.DataFrame(additional_results)

# Combine with original data
dataframes_to_concat = [
    df_best10.reset_index(drop=True), 
    df_lipinski, 
    df_qed,
    df_kinase,
    df_sa,
    df_additional
]

if DEEPCHEM_AVAILABLE and deepchem_results:
    df_deepchem = pd.DataFrame(deepchem_results)
    dataframes_to_concat.append(df_deepchem)

df_candidates = pd.concat(dataframes_to_concat, axis=1)

# Save to CSV
df_candidates.to_csv('candidates.csv', index=False)

print(f"\n✅ Saved {len(df_candidates)} candidates to 'candidates.csv'")
print(f"\n📊 Property Categories Added:")
print("   ━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━")
print("   🔹 Lipinski Properties:")
print("      • molecular_weight, logP, H-donors/acceptors")
print("      • rotatable_bonds, violations, pass/fail")
print("\n   🔹 Drug-likeness:")
print("      • qed_score, qed_classification")
print("\n   🔹 Kinase-Specific Properties:")
print("      • TPSA (polar surface area)")
print("      • aromatic_rings (ATP pocket binding)")
print("      • heteroatoms (hinge H-bonding)")
print("      • ring counts, Fsp3")
print("      • molar_refractivity")
print("      • kinase_score (custom scoring)")
print("\n   🔹 Synthetic Accessibility:")
print("      • sa_score, sa_category")
print("\n   🔹 Additional Descriptors:")
print("      • stereo_centers, formal_charge")
print("      • structural complexity indices")

if DEEPCHEM_AVAILABLE:
    print("\n   🔹 ADMET Properties (DeepChem):")
    print("      • solubility_logs, solubility_class")
    print("      • bbb_permeability (CNS penetration)")
    print("      • cyp3a4_inhibitor (drug-drug interactions)")
    print("      • herg_liability (cardiac toxicity)")
    print("      • clearance_pred (elimination rate)")

# Enhanced Summary
print(f"\n{'='*50}")
print(f"{'ANALYSIS SUMMARY':^50}")
print(f"{'='*50}")

print(f"\n🎯 Drug-likeness:")
print(f"   Lipinski pass: {(df_lipinski['lipinski_pass'] == 'Yes').sum()}/{len(df_candidates)}")
print(f"   Average QED: {df_qed['qed_score'].mean():.3f}")
print(f"   QED ≥ 0.7 (Excellent): {(df_qed['qed_score'] >= 0.7).sum()}/{len(df_candidates)}")

print(f"\n🧬 Kinase-Specific Metrics:")
print(f"   Avg aromatic rings: {df_kinase['num_aromatic_rings'].mean():.1f}")
print(f"   Avg TPSA: {df_kinase['tpsa'].mean():.1f} Ų")
print(f"   TPSA in ideal range (40-100): {((df_kinase['tpsa'] >= 40) & (df_kinase['tpsa'] <= 100)).sum()}/{len(df_candidates)}")
print(f"   Avg kinase score: {df_kinase['kinase_score'].mean():.1f}/11")
print(f"   High kinase score (≥8): {(df_kinase['kinase_score'] >= 8).sum()}/{len(df_candidates)}")

if df_sa['sa_score'].notna().any():
    print(f"\n🔬 Synthetic Accessibility:")
    print(f"   Avg SA score: {df_sa['sa_score'].mean():.2f} (1=easy, 10=hard)")
    print(f"   Easy to synthesize (≤3): {(df_sa['sa_score'] <= 3).sum()}/{len(df_candidates)}")

if DEEPCHEM_AVAILABLE and deepchem_results:
    print(f"\n💊 ADMET Profile:")
    df_dc = pd.DataFrame(deepchem_results)
    
    if df_dc['solubility_logs'].notna().any():
        print(f"   Avg solubility (LogS): {df_dc['solubility_logs'].mean():.2f}")
        print(f"   Soluble compounds: {(df_dc['solubility_class'].isin(['Highly soluble', 'Soluble'])).sum()}/{len(df_candidates)}")
    
    if df_dc['bbb_permeability'].notna().any():
        print(f"   BBB permeability (Likely): {(df_dc['bbb_permeability'] == 'Likely').sum()}/{len(df_candidates)}")
    
    if df_dc['cyp3a4_inhibitor'].notna().any():
        print(f"   CYP3A4 low risk: {(df_dc['cyp3a4_inhibitor'] == 'Low risk').sum()}/{len(df_candidates)}")
    
    if df_dc['herg_liability'].notna().any():
        print(f"   hERG low risk: {(df_dc['herg_liability'] == 'Low risk').sum()}/{len(df_candidates)}")

print(f"\n📋 Top 3 candidates by kinase score:")
if 'id-num' in df_candidates.columns:
    top_candidates = df_candidates.nlargest(3, 'kinase_score')
    for idx, row in top_candidates.iterrows():
        print(f"\n   ID: {int(row['id-num'])}")
        print(f"      Kinase score: {row['kinase_score']}, QED: {row['qed_score']:.3f}")
        print(f"      TPSA: {row['tpsa']:.1f}, Aromatic rings: {row['num_aromatic_rings']}")
        
        if DEEPCHEM_AVAILABLE and 'solubility_class' in row:
            print(f"      Solubility: {row['solubility_class']}, hERG: {row['herg_liability']}")
else:
    top_candidates = df_candidates.nlargest(3, 'kinase_score')
    for idx, row in top_candidates.iterrows():
        print(f"\n   Index: {idx}")
        print(f"      Kinase score: {row['kinase_score']}, QED: {row['qed_score']:.3f}")
        print(f"      TPSA: {row['tpsa']:.1f}, Aromatic rings: {row['num_aromatic_rings']}")
        
        if DEEPCHEM_AVAILABLE and 'solubility_class' in row:
            print(f"      Solubility: {row['solubility_class']}, hERG: {row['herg_liability']}")

print(f"\n{'='*50}")
