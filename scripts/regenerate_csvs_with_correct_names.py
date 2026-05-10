#!/usr/bin/env python3
"""
Regenerate CSV files with correct lipid names by fixing corrupted data
This script reads the corrupted CSV files and fixes them using the correct mapping
"""

import pandas as pd
import numpy as np
import os
import re
import glob

# Set paths
base_dir = "/Users/senol/Desktop/GBM/Lipidomics_Metabolomics_Integration"
results_dir = f"{base_dir}/results"

print("="*80)
print("REGENERATING CSV FILES WITH CORRECT LIPID NAMES")
print("="*80)

# Create correct mapping (same as before)
print("\n[1/3] Creating correct lipid name mapping...")

lipidome_pos = pd.read_csv("/Users/senol/Desktop/GBM/lipidome_positive_quantification.txt", sep="\t")
lipidome_neg = pd.read_csv("/Users/senol/Desktop/GBM/lipidome_negative_quantification.txt", sep="\t")
metabolome = pd.read_csv("/Users/senol/Desktop/GBM/metabolome_quantification.txt", sep="\t")
clinical = pd.read_csv("/Users/senol/Desktop/GBM/clinical_patient.txt", sep="\t", skiprows=4)
clinical.columns = clinical.columns.str.strip()
clinical = clinical[clinical['PATIENT_ID'].notna()]

def extract_sample_ids(df):
    return [col for col in df.columns if col.startswith('C3')]

def prepare_data_matrix(df, sample_cols, id_col='ENTITY_STABLE_ID', name_col='NAME'):
    if id_col in df.columns:
        df = df.set_index(id_col)
    elif name_col in df.columns:
        df = df.set_index(name_col)
    else:
        df = df.set_index(df.columns[0])
    available_samples = [s for s in sample_cols if s in df.columns]
    df_subset = df[available_samples].T
    df_subset = df_subset.apply(pd.to_numeric, errors='coerce')
    missing_threshold = 0.5
    df_subset = df_subset.loc[:, df_subset.isna().mean() < missing_threshold]
    df_subset = df_subset.fillna(df_subset.median())
    return df_subset

metabolome_samples = extract_sample_ids(metabolome)
lipidome_pos_samples = extract_sample_ids(lipidome_pos)
lipidome_neg_samples = extract_sample_ids(lipidome_neg)
clinical_patient_ids = clinical['PATIENT_ID'].dropna().unique().tolist()

common_samples = sorted(list(set(metabolome_samples) & set(lipidome_pos_samples) & 
                     set(lipidome_neg_samples) & set(clinical_patient_ids)))

lipidome_pos_matrix = prepare_data_matrix(lipidome_pos, common_samples)
lipidome_neg_matrix = prepare_data_matrix(lipidome_neg, common_samples)

# Create correct mapping
lipid_mapping = {}
lipid_idx = 0

for col in lipidome_pos_matrix.columns:
    if col in lipidome_pos['ENTITY_STABLE_ID'].values:
        name = lipidome_pos.loc[lipidome_pos['ENTITY_STABLE_ID'] == col, 'NAME'].values[0]
    else:
        name = col
    lipid_mapping[f'LIPID_{lipid_idx}'] = name
    lipid_idx += 1

for col in lipidome_neg_matrix.columns:
    if col in lipidome_neg['ENTITY_STABLE_ID'].values:
        name = lipidome_neg.loc[lipidome_neg['ENTITY_STABLE_ID'] == col, 'NAME'].values[0]
    else:
        name = col
    lipid_mapping[f'LIPID_{lipid_idx}'] = name
    lipid_idx += 1

# Create reverse mapping for corrupted patterns
# Pattern: "carnitine(14:1)33" means LIPID_33 was replaced, leaving "3" -> should be LIPID_333
reverse_corrupted_mapping = {}
for old_name, new_name in lipid_mapping.items():
    # Check if this name could have been corrupted
    # For example, if LIPID_333 maps to "Cer(d16:1/18:0);Cer(d18:1/16:0)"
    # and we see "carnitine(14:1)33", it means LIPID_33 was replaced with "carnitine(14:1)"
    # and "3" was left, so the original was LIPID_333
    pass  # This is complex, we'll use a different approach

print(f"Created mapping for {len(lipid_mapping)} lipids")

# ============================================================================
# 2. FIX CSV FILES BY DETECTING CORRUPTED PATTERNS
# ============================================================================
print("\n[2/3] Fixing corrupted CSV files...")

# Strategy: For each corrupted pattern like "carnitine(14:1)33",
# we need to figure out what LIPID_XXX it should be
# The pattern suggests: LIPID_33 was replaced, leaving "3" -> original was LIPID_333

# Create a mapping of what each LIPID_XX maps to
lipid_name_to_id = {name: lid for lid, name in lipid_mapping.items()}

csv_files = glob.glob(f"{results_dir}/*.csv")
print(f"Found {len(csv_files)} CSV files")

for csv_file in csv_files:
    print(f"  Processing: {os.path.basename(csv_file)}")
    try:
        df = pd.read_csv(csv_file)
        original_df = df.copy()
        
        # Check each column for corrupted patterns
        for col in df.columns:
            if df[col].dtype == 'object':
                col_data = df[col].astype(str)
                
                # First, try to find any remaining LIPID_XXX patterns and replace them
                for old_name, new_name in sorted(lipid_mapping.items(), key=lambda x: len(x[0]), reverse=True):
                    pattern = r'\b' + re.escape(old_name) + r'\b'
                    col_data = col_data.str.replace(pattern, new_name, regex=True)
                
                # Now try to fix corrupted patterns
                # Pattern: "name)XX" where XX are digits at the end
                # This suggests "LIPID_XX" was replaced, leaving digits
                corrupted_pattern = r'([^)]+)\)(\d+)$'
                
                def fix_corrupted(match):
                    name_part = match.group(1) + ")"
                    digits = match.group(2)
                    
                    # Try to find which LIPID_XXX this should be
                    # The digits suggest the original was LIPID_XXX where XXX ends with these digits
                    # For example, "carnitine(14:1)33" -> should be LIPID_333
                    # We need to find what LIPID_XX maps to "name_part", then append digits
                    
                    # Find all LIPID_XX that map to something starting with name_part
                    candidates = []
                    for lid, name in lipid_mapping.items():
                        if name.startswith(name_part):
                            # Extract the number from LIPID_XXX
                            lid_num = lid.replace('LIPID_', '')
                            # Check if lid_num ends with digits
                            if lid_num.endswith(digits):
                                candidates.append((lid, name))
                    
                    if len(candidates) == 1:
                        return candidates[0][1]  # Return the correct name
                    elif len(candidates) > 1:
                        # Multiple candidates, use the one that matches best
                        # Prefer the one where the full LIPID_XXX number matches
                        for lid, name in candidates:
                            lid_num = lid.replace('LIPID_', '')
                            if lid_num == digits or lid_num.endswith(digits):
                                return name
                        return candidates[0][1]  # Fallback to first
                    else:
                        # No match found, return original
                        return match.group(0)
                
                # Apply fix to corrupted patterns
                col_data = col_data.str.replace(corrupted_pattern, fix_corrupted, regex=True)
                
                df[col] = col_data
        
        # Also fix column names
        new_columns = {}
        for col in df.columns:
            new_col = col
            for old_name, new_name in sorted(lipid_mapping.items(), key=lambda x: len(x[0]), reverse=True):
                if old_name in new_col:
                    new_col = new_col.replace(old_name, new_name)
            if new_col != col:
                new_columns[col] = new_col
        if new_columns:
            df = df.rename(columns=new_columns)
        
        # Save if changed
        if not df.equals(original_df):
            df.to_csv(csv_file, index=False)
            print(f"    ✓ Fixed and saved")
        else:
            print(f"    - No changes needed")
            
    except Exception as e:
        print(f"    ✗ Error: {e}")
        import traceback
        traceback.print_exc()

print("\n[3/3] Note: Some files may still need manual review or regeneration")
print("="*80)
print("FIXING COMPLETE")
print("="*80)







