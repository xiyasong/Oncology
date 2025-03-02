# cmd + / 
# Read basic database
import json
import os
import pandas as pd
import argparse
import gzip
import shutil


# Define argument parser
# parser = argparse.ArgumentParser(description="Pipeline for processing FDA and gene panel data with JSON input.")
# parser.add_argument("--json", required=True, help="Path to the input JSON file.")
# parser.add_argument("--path", choices=["local", "sysmed"], default="local", help="Choose the file system to use.")
# parser.add_argument("--output", required=True, help="Path for the output simplified FDA file.")

# args = parser.parse_args()
# json_path = args.json
# output_path = args.output


input_folder = "/mnt/storage_pool/xiya/Oncology/input_pool"
finished_folder = "/mnt/storage_pool/xiya/Oncology/input_pool_finished"
output_folder = "/mnt/storage_pool/xiya/Oncology/output_pool"
# Ensure output and finished folders exist
os.makedirs(output_folder, exist_ok=True)
os.makedirs(finished_folder, exist_ok=True)



############ ipynb file 
fda_path = '/mnt/storage_pool/xiya/Oncology/Newest-List of Cleared or Approved Companion Diagnostic Devices (In Vitro and Imaging Tools)  FDA.xlsx'
TSO_523_list_file = '/mnt/storage_pool/xiya/Oncology/gene_list_trusight_oncology_500.xlsx'
gene_panel_file_path = '/mnt/storage_pool/xiya/Oncology/Updated_Cancer_Gene_List_with_Pan-cancer.csv'
# Read the Excel file
fda_df = pd.read_excel(fda_path,header = 1)
#print("fda_df.shape",fda_df.shape)
# Read the Excel file
TSO_523_list = pd.read_excel(TSO_523_list_file,header = 0)
TSO_523_genes = set(TSO_523_list.loc[TSO_523_list['Small variants'] == 'ü', 'Gene symbol'])
#print("TSO_523_list.shape",TSO_523_list.shape)

gene_panel_df = pd.read_csv(gene_panel_file_path, header = 0)
#gene_panel_df.to_csv('/Users/xiyas/Oncology/Results-unknown-cancer-solid/P2_FreshFrozenTissue_DNA/P2_FreshFrozenTissue_DNA/gene_panel_df.txt', sep='\t', index=False)
gene_panel_expanded = gene_panel_df.copy()
gene_panel_expanded['Genes'] = gene_panel_expanded['Genes'].str.split(', ')
gene_panel_expanded = gene_panel_expanded.explode('Genes').reset_index(drop=True)
#print("gene_panel_expanded.shape ",gene_panel_expanded.shape)

#(54, 9)
# Define functions 
def amino_acid_converter(sequence, to_three_letter=True):
    # Mapping dictionaries for one-letter to three-letter and vice versa
    one_to_three = {
        'A': 'Ala', 'R': 'Arg', 'N': 'Asn', 'D': 'Asp', 'C': 'Cys',
        'Q': 'Gln', 'E': 'Glu', 'G': 'Gly', 'H': 'His', 'I': 'Ile',
        'L': 'Leu', 'K': 'Lys', 'M': 'Met', 'F': 'Phe', 'P': 'Pro',
        'S': 'Ser', 'T': 'Thr', 'W': 'Trp', 'Y': 'Tyr', 'V': 'Val'
    }
    three_to_one = {v: k for k, v in one_to_three.items()}

    if not to_three_letter:
        # Convert from three-letter to one-letter code, ignore anything else
        for three_letter, one_letter in three_to_one.items():
            sequence = sequence.replace(three_letter, one_letter)
        return sequence
    return sequence

# Example usage
one_letter_sequence = "ARND"
three_letter_sequence = "NP_006258.3:p.(Leu739Met)"
mutation = "V600E"

# print(amino_acid_converter(one_letter_sequence, to_three_letter=True))  # Output: Ala-Arg-Asn-Asp
# print(amino_acid_converter(three_letter_sequence, to_three_letter=False))  # Output: ARND
# print(amino_acid_converter(mutation, to_three_letter=True))  # Output: Val600Glu

def filter_consequences(row):
    # Check if "pathogenic" is in ClinVar Significance; treat empty as non-pathogenic
    has_pathogenic = any("pathogenic" in s for s in row['ClinVar Significance']) if row['ClinVar Significance'] else False
    # If "pathogenic" exists, keep the row
    if has_pathogenic:
        return True  # Keep the row if it contains "pathogenic" 
    # Define the skip terms for skippable consequences
    skip_terms = ["synonymous_variant", "3_prime_UTR_variant", "5_prime_UTR_variant", "intron_variant"]
    # Flatten the list of lists in consequence
    flattened_consequences = [item for sublist in row['consequence'] for item in sublist]
    # Check if all consequences are in skip_terms
    all_consequences_skippable = all(consequence in skip_terms for consequence in flattened_consequences)
    # Only keep the row if not all consequences are skippable
    return not all_consequences_skippable  # Keep the row if not all consequences are skippable

## define a function to matching the FDA table and returns positive/negative 
def map_biomarkers(row, oncology_df_filtered):
    marker = row['rsID genomics positions']
    biomarker_gene = row['Biomarker(s)']
    print("screening marker:", marker + " which is inside gene:", biomarker_gene)
    matching_indices = []  # To store indices of matching rows in oncology_df
    match_details = []     # To store formatted match details (Chromosome, Position, Reference Allele, Alternative Alleles)

    # Skip markers related to RNA, TMB, or MSI ---- these are handled separately
    if marker.startswith('RNA') or 'TMB' in marker or 'MSI' in marker:
        print("Skip RNA fusion target, TMB and MSI value")
        return 'skip', "", matching_indices

    # Gene-specific markers
    elif marker.startswith('Gene'):
        gene_names = marker.split(':')[1].split(',')
        print("testing target gene_names:", gene_names)

        # Find matching rows in oncology_df based on gene names
        gene_matches = oncology_df_filtered[oncology_df_filtered['HGNC'].apply(lambda patient_genes: any(gene in patient_genes for gene in gene_names))]
        if not gene_matches.empty:
            matching_indices = gene_matches.index.tolist()
            match_details = [
                f"{row['Chromosome']}, {row['Position']}, {row['Reference Allele']}, {row['Alternative Alleles']}, {row['Variant Frequencies']},{row['Unique_HGNC']}, {row['HGVSc']}, {row['HGVSp']}"
                for _, row in gene_matches.iterrows()
            ]
            return 'positive', "; ".join(match_details), matching_indices
        else:
            return 'negative', "", matching_indices

    # Amino acid-specific changes
    elif marker.startswith('AA'):
        aa_change = marker.split(':')[1]
        print("testing target aa_change:", aa_change)

        # First filter by gene match, then by amino acid change
        gene_matches = oncology_df_filtered[oncology_df_filtered['HGNC'].apply(lambda patient_genes: biomarker_gene in patient_genes)]
        if not gene_matches.empty:
            aa_matches = gene_matches[gene_matches['HGVSp_Converted'].apply(lambda changes: any(aa_change in change for change in changes))]
            if not aa_matches.empty:
                matching_indices = aa_matches.index.tolist()
                match_details = [
                    f"{row['Chromosome']}, {row['Position']}, {row['Reference Allele']}, {row['Alternative Alleles']}, {row['Variant Frequencies']}, {row['HGNC']}, {row['HGVSc']}, {row['HGVSp']}"
                    for _, row in aa_matches.iterrows()
                ]
                return 'positive', "; ".join(match_details), matching_indices
        return 'negative', "", matching_indices

    # Exon 19 or specific mutation (L858R) conditions
    elif marker.startswith('Exon 19 deletion or exon 21 L858R substitution mutation'):
        exon_19_gene_matches = oncology_df_filtered[oncology_df_filtered['HGNC'].apply(lambda patient_genes: biomarker_gene in patient_genes)]
        if not exon_19_gene_matches.empty:
            exon_19_matches = exon_19_gene_matches[exon_19_gene_matches['exons'].apply(lambda exons: '19' in exons if isinstance(exons, list) else False)]
            l858r_matches = exon_19_gene_matches[exon_19_gene_matches['HGVSp_Converted'].apply(lambda changes: any('L858R' in change for change in changes))]
        
            if not exon_19_matches.empty or not l858r_matches.empty:
                matching_indices = list(set(exon_19_matches.index.tolist() + l858r_matches.index.tolist()))
                match_details = [
                    f"{row['Chromosome']}, {row['Position']}, {row['Reference Allele']}, {row['Alternative Alleles']}, {row['Variant Frequencies']}, {row['HGNC']},{row['HGVSc']}, {row['HGVSp']}"
                    for _, row in exon_19_gene_matches.iterrows()
                ]
                return 'positive', "; ".join(match_details), matching_indices
            else:
                return 'negative', "", matching_indices

    # Handling for KRAS wild-type (absence of mutations in codons 12 and 13)
    elif marker.startswith('WildType:KRAS(12,13)'):
        kras_mutant = oncology_df_filtered[
            (oncology_df_filtered['HGNC'].apply(lambda genes: 'KRAS' in genes)) &
            (oncology_df_filtered['HGVSp_Converted'].apply(lambda mutations: any('12' in mut or '13' in mut for mut in mutations)))
        ]
        if not kras_mutant.empty:
            matching_indices = kras_mutant.index.tolist()
            match_details = [
                f"{row['Chromosome']}, {row['Position']}, {row['Reference Allele']}, {row['Alternative Alleles']}, {row['Variant Frequencies']}, {row['HGNC']},{row['HGVSc']}, {row['HGVSp']}"
                for _, row in kras_mutant.iterrows()
            ]
            return 'negative', "; ".join(match_details), matching_indices
        else:
            return 'positive', "", []

    # Handling for KRAS and NRAS wild-type (absence of mutations in exons 2, 3, and 4)
    elif marker.startswith('WildType:KRAS(2,3,4,NRAS)'):
        kras_exons_wildtype = oncology_df_filtered[
            (oncology_df_filtered['HGNC'].apply(lambda genes: 'KRAS' in genes)) &
            (~oncology_df_filtered['exons'].apply(lambda exons: any(exon in ['2', '3', '4'] for exon in exons if isinstance(exons, list))))
        ]
        nras_exons_wildtype = oncology_df_filtered[
            (oncology_df_filtered['HGNC'].apply(lambda genes: 'NRAS' in genes)) &
            (~oncology_df_filtered['exons'].apply(lambda exons: any(exon in ['2', '3', '4'] for exon in exons if isinstance(exons, list))))
        ]
        if not kras_exons_wildtype.empty or not nras_exons_wildtype.empty:
            matching_indices = list(set(kras_exons_wildtype.index.tolist() + nras_exons_wildtype.index.tolist()))
            match_details = [
                f"{row['Chromosome']}, {row['Position']}, {row['Reference Allele']}, {row['Alternative Alleles']}, {row['Variant Frequencies']}, {row['HGNC']},{row['HGVSc']}, {row['HGVSp']}"
                for _, row in pd.concat([kras_exons_wildtype, nras_exons_wildtype]).iterrows()
            ]
            return 'negative', "; ".join(match_details), matching_indices
        else:
            return 'positive', "", []

    return 'negative', "; ".join(match_details), matching_indices

def map_oncology_to_fda(fda_df_simplified, oncology_df_filtered, cancer_type_specific_gene_panel):
    # Step 1: Filter oncology_df based on genes of interest, ensuring unique HGNC values
    oncology_filtered_sec3 = oncology_df_filtered[oncology_df_filtered['Unique_HGNC'].apply(lambda gene_list: any(gene in cancer_type_specific_gene_panel for gene in gene_list))]

    # Step 2: Prepare new rows with required details for Section 3
    section_3_rows = []
    for idx, row in oncology_filtered_sec3.iterrows():
        # Get unique genes that match in the HGNC list
        matching_genes = [gene for gene in row['Unique_HGNC'] if gene in cancer_type_specific_gene_panel]
        # Construct Biomarker(s) Details
        biomarker_details = ', '.join(matching_genes)
        # Set Mapping Result as "positive"
        mapping_result = "positive" 
        # Generate Mapping Details string
        mapping_details = f"{row['Chromosome']}, {row['Position']}, {row['Reference Allele']}, " \
                          f"{row['Alternative Alleles']}, {row['Variant Frequencies']}, " \
                          f"{row['HGNC']}, {row['HGVSc']}, {row['HGVSp']}"
        
        # Get matching indices
        matching_indices = [idx]  # Using the current index as a placeholder, adjust if needed
        # Add a new row to the list for Section 3
        section_3_rows.append({
            "Biomarker(s) (Details)": biomarker_details,
            "Mapping Result": mapping_result,
            "Mapping Details": mapping_details,
            "Mapping Indices": matching_indices,
            "Section": "3"
        })
    # Step 3: Convert section 3 rows to DataFrame
    section_3_df = pd.DataFrame(section_3_rows)
    # Step 4: Prepare rows for Section 4 with unmatched variants
    oncology_unmatched = oncology_df_filtered[~oncology_df_filtered.index.isin(oncology_filtered_sec3.index)]
    section_4_rows = []
    for idx, row in oncology_unmatched.iterrows():
        # Biomarker(s) Details will include all genes in Unique_HGNC, as none matched
        biomarker_details = ', '.join(row['Unique_HGNC'])
        # Set Mapping Result as "negative" for unmatched
        mapping_result = "positive"
        # Generate Mapping Details string
        mapping_details = f"{row['Chromosome']}, {row['Position']}, {row['Reference Allele']}, " \
                          f"{row['Alternative Alleles']}, {row['Variant Frequencies']}, " \
                          f"{row['HGNC']}, {row['HGVSc']}, {row['HGVSp']}"
        
        # Use the current index for unmatched rows
        matching_indices = [idx]
        # Add a new row to the list for Section 4
        section_4_rows.append({
            "Biomarker(s) (Details)": biomarker_details,
            "Mapping Result": mapping_result,
            "Mapping Details": mapping_details,
            "Mapping Indices": matching_indices,
            "Section": "4"
        })

    # Step 5: Convert section 4 rows to DataFrame
    section_4_df = pd.DataFrame(section_4_rows)
    # Step 6: Append both sections to fda_df
    fda_df_simplified = pd.concat([fda_df_simplified, section_3_df, section_4_df], ignore_index=True)
    return fda_df_simplified

# 3. Load JSON data into a df that contains the same content as CombinedVariantOutput.tsv 

def open_json_file(json_path):
    if json_path.endswith(".gz"):
        with gzip.open(json_path, 'rt', encoding='utf-8') as f:
            return json.load(f)
    else:
        with open(json_path, 'r', encoding='utf-8') as f:
            return json.load(f)

def process_json_file(input_path, output_path):
    data = open_json_file(input_path)

    # Initialize records list
    records = []
    # Loop through each position and extract relevant information
    for position in data['positions']:
        # Check if 'filters' field contains "PASS" to keep only the PASSED variants
        if "PASS" in position.get('filters', []):
            for variant in position['variants']:
                # Check if ClinVar Significance contains any "benign" or "likely benign" to skip
                clinvar_significance = [c.get('significance', []) for c in variant.get('clinvar', [])]
                clinvar_review_status = [c.get('reviewStatus', []) for c in variant.get('clinvar', [])]
                hgvsc_values = []
                hgvsp_values = []
                hgvsp_converted_values = []
                hgnc_values = []
                exons_values = []
                consequence_values = []
                is_canonical = []
                # Extracting the transcript-level details for hgvsc, hgvsp, and consequence
                for transcript in variant.get('transcripts', []):
                    if transcript.get('isCanonical', False) and transcript.get('hgnc') in TSO_523_genes:# this is to make sure take the correct transcripts
                        # If not skipped, collect the transcript data
                        hgvsc_values.append(transcript.get('hgvsc', ""))
                        hgvsp_values.append(transcript.get('hgvsp', ""))
                        hgvsp_converted_values.append([amino_acid_converter(hgvsp.split('.')[-1].replace('(', '').replace(')', ''), to_three_letter=False) for hgvsp in hgvsp_values])
                        hgnc_values.append(transcript.get('hgnc', ""))
                        exons_values.append(transcript.get('exons', ""))
                        consequence_values.append(transcript.get('consequence', []))
                        is_canonical.append(transcript.get('isCanonical', False))

                # Creating a record for each variant
                record = [
                    position['chromosome'],
                    position['position'],
                    position['refAllele'],
                    position['altAlleles'],
                    position['filters'],
                    position['cytogeneticBand'],
                    position['samples'][0].get('genotype'),
                    position['samples'][0].get('variantFrequencies'),
                    position['samples'][0].get('totalDepth'),
                    position['samples'][0].get('alleleDepths'),
                    variant.get('vid'),
                    variant.get('variantType'),
                    variant.get('hgvsg'),
                    hgvsc_values,
                    hgvsp_values,
                    hgvsp_converted_values,
                    hgnc_values,
                    exons_values,
                    consequence_values,
                    is_canonical,
                    variant.get('ancestralAllele'),
                    [c.get('id') for c in variant.get('clinvar', [])],
                    clinvar_significance,
                    clinvar_review_status,
                    [c.get('id') for c in variant.get('cosmic', [])],
                    [c.get('cancerTypesAndCounts') for c in variant.get('cosmic', [])],
                    variant.get('dbsnp')
                ]
                records.append(record)

    # Define column names for the DataFrame
    columns = [
        "Chromosome", "Position", "Reference Allele", "Alternative Alleles", "Filter",
        "Cytogenetic Band", "Genotype", "Variant Frequencies", "Total Depth", "Allele Depths",
        "Variant ID", "Variant Type", "HGVSg", "HGVSc", "HGVSp", "HGVSp_Converted", "HGNC", "exons", "consequence", "is_canonical", "Ancestral Allele", "ClinVar IDs",
        "ClinVar Significance", "ClinVar Review Status","COSMIC IDs", "COSMIC Cancer Types", "dbSNP IDs"
    ]

    # Create DataFrame from records
    oncology_df = pd.DataFrame(records, columns=columns)
    #oncology_df.shape
    # oncology_df contains the same amount of records with CombinedVariantsOutputs (but with annotations from NIRVANA).
    # 4. start the filter applied to oncology_df
    # Remove duplicates within each HGNC list
    oncology_df['Unique_HGNC'] = oncology_df['HGNC'].apply(lambda x: list(set(x))) 

    # Filter the Unique_HGNC list in each row to keep only genes name in TSO_523_genes
    oncology_df['Unique_HGNC'] = oncology_df['Unique_HGNC'].apply(lambda genes: [gene for gene in genes if gene in TSO_523_genes])
    # Step 2: Remove rows where Unique_HGNC is empty (no genes matched with TSO_523_genes)
    oncology_df_filtered = oncology_df[oncology_df['Unique_HGNC'].str.len() > 0]

    print("oncology_df.shape ",oncology_df.shape)
    #sample 257: 667,26
    print("oncology_df_filtered.shape (removed the variants with empty gene name) ",oncology_df_filtered.shape)

    ## Filter 1: Remove known benign variants 
    oncology_df_filtered = oncology_df_filtered[
        oncology_df_filtered['ClinVar Significance'].apply(lambda significance: all("benign" not in s for s in significance))
    ]
    print("Shape after filtering for benign significance: ", oncology_df_filtered.shape)

    ## Filter 2: Remove variants based on molecular consequences, but only if it is not known pathogenic
    skip_terms = ["synonymous_variant", "3_prime_UTR_variant", "5_prime_UTR_variant", "intron_variant"]
    # Identify rows with "pathogenic" in ClinVar Significance
    # Apply the filtering function
    oncology_df_filtered = oncology_df_filtered[oncology_df_filtered.apply(filter_consequences, axis=1)]
    # Display the final shape after filtering
    print("Final filtered shape:", oncology_df_filtered.shape)

    # 6. Apply the function to each row in the FDA DataFrame and unpack the result, mapping details, and matching indices

    fda_df[['Mapping Result', 'Mapping Details', 'Mapping Indices']] = fda_df.apply(
        lambda row: pd.Series(map_biomarkers(row, oncology_df_filtered)), axis=1
    )
    #simplifying the fda_df 
    # only keep Biomarker(s) (Details),Mapping Result,Mapping Details,Mapping Indices
    # Keep only the specified columns in fda_df
    fda_df_simplified = fda_df[['Biomarker(s) (Details)','Mapping Result', 'Mapping Details', 'Mapping Indices']]
    # Filter to keep only rows where 'Mapping Result' is 'positive'
    fda_df_simplified = fda_df_simplified[fda_df_simplified['Mapping Result'] == 'positive']
    # Keep only unique rows in the DataFrame
    # Convert any unhashable columns (e.g., lists) to strings
    fda_df_simplified = fda_df_simplified.applymap(lambda x: str(x) if isinstance(x, list) else x)
    fda_df_simplified = fda_df_simplified.drop_duplicates()

    fda_df_simplified['Section'] = "1"
    # Display the simplified DataFrame
    #print(fda_df_simplified.head())
    # 7. Adding section 3 and 4
    # Section 3: cancer-type specific 
    # Section 4: All other variants

    cancer_type_specific_gene_panel = set(gene_panel_expanded['Genes'])
    print(len(cancer_type_specific_gene_panel))
    fda_df_simplified_new = map_oncology_to_fda(fda_df_simplified, oncology_df_filtered, cancer_type_specific_gene_panel)

    fda_df_simplified_new.shape
    section_counts = fda_df_simplified_new['Section'].value_counts()

    # Retrieve counts for sections "1" and "3" specifically
    section_1_count = section_counts.get("1", 0)  # Default to 0 if section "1" is not present
    section_3_count = section_counts.get("3", 0)  # Default to 0 if section "3" is not present
    section_4_count = section_counts.get("4", 0)  # Default to 0 if section "3" is not present

    # Display results
    print(f"Shape of fda_df_simplified_new: {fda_df_simplified_new.shape}")
    print(f"Count of Section 1: {section_1_count}")
    print(f"Count of Section 3: {section_3_count}")
    print(f"Count of Section 4: {section_4_count}")

    # final output
    fda_df_simplified_new.to_csv(output_path, sep='\t', index=False)


    # Process all json.gz files in the input folder
for filename in os.listdir(input_folder):
    if filename.endswith(".json.gz") or filename.endswith(".json"):
        input_path = os.path.join(input_folder, filename)
        output_path = os.path.join(output_folder, f"{os.path.splitext(filename)[0]}.txt")

        try:
            # Process the JSON.gz file
            process_json_file(input_path, output_path)
            print(f"Processed: {filename} -> Output: {output_path}")

            # Move the input file to the finished folder
            shutil.move(input_path, os.path.join(finished_folder, filename))
            print(f"Moved {filename} to finished folder.")
        except Exception as e:
            print(f"Error processing {filename}: {e}")
