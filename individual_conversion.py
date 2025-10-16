#!/usr/bin/env python3
"""
Individual Conversion Script
Converts PLINK recoded data to format suitable for ancestry analysis
"""

import pandas as pd
import numpy as np
import sys
import os

def convert_genotypes():
    """
    Convert PLINK recoded .raw file to individual_genotypes2.input format
    """
    
    # File paths - will be dynamically replaced
    raw_file = "individual_1.raw"
    bim_file = "HT.bim"  # This will be replaced dynamically with actual sample_id.bim
    output_file = "individual_genotypes2.input"
    
    print(f"Reading raw file: {raw_file}")
    print(f"Reading BIM file: {bim_file}")
    
    try:
        # Read the raw file (PLINK recoded format)
        if not os.path.exists(raw_file):
            raise FileNotFoundError(f"Raw file not found: {raw_file}")
            
        raw_data = pd.read_csv(raw_file, sep=r'\s+', engine='python')
        print(f"Raw data shape: {raw_data.shape}")
        
        # Get sample information
        sample_cols = ['FID', 'IID', 'PAT', 'MAT', 'SEX', 'PHENOTYPE']
        sample_info = raw_data[sample_cols]
        
        # Get genotype columns (everything after the sample info columns)
        genotype_cols = [col for col in raw_data.columns if col not in sample_cols]
        genotype_data = raw_data[genotype_cols]
        
        print(f"Found {len(genotype_cols)} variants")
        
        # Read BIM file for variant information
        if os.path.exists(bim_file):
            bim_data = pd.read_csv(bim_file, sep='\t', header=None, 
                                   names=['chr', 'rsid', 'cm', 'pos', 'a1', 'a2'])
            print(f"BIM data shape: {bim_data.shape}")
            
            # Create variant IDs matching the raw file column names
            bim_data['var_id'] = bim_data['rsid'].astype(str) + '_' + bim_data['a1'].astype(str)
        else:
            print(f"Warning: BIM file not found: {bim_file}")
            bim_data = None
        
        # Prepare output format
        output_lines = []
        
        # Header line with variant information
        if bim_data is not None:
            # Match variants between raw file and BIM file
            header_parts = ['IID']
            for col in genotype_cols:
                # Extract chromosome and position from column name if possible
                if '_' in col:
                    parts = col.split('_')
                    if len(parts) >= 2:
                        header_parts.append(col)
                else:
                    header_parts.append(col)
            
            output_lines.append('\t'.join(header_parts))
        else:
            # Use column names as-is
            header_parts = ['IID'] + genotype_cols
            output_lines.append('\t'.join(header_parts))
        
        # Add genotype data for each sample
        for idx, row in raw_data.iterrows():
            sample_id = row['IID']
            genotypes = []
            
            for col in genotype_cols:
                # Convert NA/missing to appropriate format
                if pd.isna(row[col]):
                    genotypes.append('NA')
                else:
                    # PLINK codes: 0=homozygous A1, 1=heterozygous, 2=homozygous A2
                    genotypes.append(str(int(row[col])))
            
            output_lines.append(f"{sample_id}\t" + '\t'.join(genotypes))
        
        # Write output file
        with open(output_file, 'w') as f:
            f.write('\n'.join(output_lines))
        
        print(f"Conversion completed successfully!")
        print(f"Output written to: {output_file}")
        print(f"Total samples: {len(raw_data)}")
        print(f"Total variants: {len(genotype_cols)}")
        
        # Create a summary file
        summary_file = "conversion_summary.txt"
        with open(summary_file, 'w') as f:
            f.write(f"Conversion Summary\n")
            f.write(f"==================\n")
            f.write(f"Input raw file: {raw_file}\n")
            f.write(f"Input BIM file: {bim_file}\n")
            f.write(f"Output file: {output_file}\n")
            f.write(f"Total samples: {len(raw_data)}\n")
            f.write(f"Total variants: {len(genotype_cols)}\n")
            f.write(f"Sample IDs: {', '.join(raw_data['IID'].astype(str).tolist())}\n")
        
        return True
        
    except Exception as e:
        print(f"Error during conversion: {e}")
        import traceback
        traceback.print_exc()
        return False

def main():
    """Main function"""
    print("Starting individual genotype conversion...")
    print(f"Working directory: {os.getcwd()}")
    
    # List files in current directory for debugging
    print("\nFiles in current directory:")
    for f in os.listdir('.'):
        if f.endswith(('.raw', '.bim', '.bed', '.fam')):
            print(f"  - {f}")
    
    success = convert_genotypes()
    
    if success:
        print("\n✓ Conversion completed successfully")
        sys.exit(0)
    else:
        print("\n✗ Conversion failed")
        sys.exit(1)

if __name__ == "__main__":
    main()