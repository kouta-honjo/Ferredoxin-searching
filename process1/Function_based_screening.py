import pandas as pd
import re
import os

# --- Define the categorization function (provided by user) ---
def categorize_product_full_redefinition(product_description):
    product_description = str(product_description).lower()
    if "hydabc" in product_description or "[fefe]-hydrogenase" in product_description or "electron-bifurcating hydrogenase" in product_description:
        return "Electron Bifurcation - HydABC"
    elif "nfn" in product_description and ("oxidoreductase" in product_description or "transhydrogenase" in product_description) and ("ferredoxin" in product_description or "nadp" in product_description):
        return "Electron Bifurcation - NfnAB"
    elif "ech complex" in product_description or "[nife]-hydrogenase" in product_description or "energy-converting hydrogenase" in product_description:
        return "Electron Bifurcation - Ech Complex"
    elif "co dehydrogenase" in product_description or "acetyl-coa synthase" in product_description or "carbon monoxide dehydrogenase" in product_description:
        return "Wood-Ljungdahl Pathway - Carbonyl Branch (CODH/ACS)"
    elif "formate dehydrogenase" in product_description:
        return "Wood-Ljungdahl Pathway - Carbonyl Branch (Formate Dehydrogenase)"
    elif "methyltransferase" in product_description or "corrinoid iron-sulfur protein" in product_description or "methylenetetrahydrofolate reductase" in product_description or "formyltetrahydrofolate synthetase" in product_description or "methylene-h4f dehydrogenase" in product_description or "folate-dependent" in product_description:
        return "Wood-Ljungdahl Pathway - Methyl Branch"
    elif "pyruvate ferredoxin oxidoreductase" in product_description or "pfor" in product_description:
        return "Key Electron Donor - PFOR"
    elif "ferredoxin" in product_description:
        return "Ferredoxin"
    elif any(keyword in product_description for keyword in ['carbohydrate', 'glycolysis', 'pentose phosphate', 'sugar', 'fructose', 'glucose', 'starch', 'cellulose', 'sucrose', 'mannose', 'galactose', 'trehalose', 'disaccharide', 'monosaccharide', 'polysaccharide', 'pyruvate', 'lactate', 'acetate', 'butanol', 'ethanol', 'glyceraldehyde', 'phosphogluconate', 'aldolase', 'isomerase', 'mutase', 'decarboxylase', 'synthase', 'kinase', 'phosphatase']):
        return "Carbohydrate Metabolism"
    elif any(keyword in product_description for keyword in ['amino acid', 'alanine', 'arginine', 'asparagine', 'aspartate', 'cysteine', 'glutamate', 'glutamine', 'glycine', 'histidine', 'isoleucine', 'leucine', 'lysine', 'methionine', 'phenylalanine', 'proline', 'serine', 'threonine', 'tryptophan', 'tyrosine', 'valine', 'biosynthesis', 'degradation', 'racemase', 'aminotransferase', 'deaminase']):
        return "Amino Acid Metabolism"
    elif any(keyword in product_description for keyword in ['nucleotide', 'purine', 'pyrimidine', 'adenine', 'guanine', 'cytosine', 'thymine', 'uracil', 'atp', 'gtp', 'ctp', 'ttp', 'utp', 'datp', 'dgtp', 'dctp', 'dttp', 'dna', 'rna', 'ribonucleotide', 'deoxyribonucleotide']):
        return "Nucleotide Metabolism"
    elif any(keyword in product_description for keyword in ['lipid', 'fatty acid', 'acylglycerol', 'phospholipid', 'glycerol', 'biosynthesis', 'degradation', 'lipase', 'acyltransferase']):
        return "Lipid Metabolism"
    elif any(keyword in product_description for keyword in ['cofactor', 'vitamin', 'nad', 'fad', 'coenzyme a', 'biotin', 'folate', 'thiamine', 'riboflavin', 'pyridoxal', 'cobalamin', 'lipoate', 'heme', 'iron-sulfur', 'molybdopterin', 'tetrahydrofolate']):
        return "Cofactor Metabolism"
    elif any(keyword in product_description for keyword in ['atp synthase', 'electron transfer', 'oxidoreductase', 'dehydrogenase', 'cytochrome', 'nadh', 'fadh2', 'hydrogenase', 'methanogenesis', 'tca cycle', 'krebs cycle', 'anaerobic respiration', 'flavodoxin']):
        return "Energy Metabolism/Electron Transport"
    elif any(keyword in product_description for keyword in ['abc transporter', 'permease', 'antiporter', 'symporter', 'efflux pump', 'channel', 'transport system', 'porin']):
        return "Transport"
    elif any(keyword in product_description for keyword in ['kinase', 'phosphatase', 'regulator', 'sensor', 'response regulator', 'transcription factor', 'signaling', 'two-component', 'gtpase', 'atpase', 'cyclase']):
        return "Regulatory/Signaling"
    elif any(keyword in product_description for keyword in ['cell wall', 'membrane', 'peptidoglycan', 'flagella', 'pilus', 'fimbrial', 's-layer']):
        return "Cell Structure/Motility"
    elif any(keyword in product_description for keyword in ['dna replication', 'dna repair', 'rna polymerase', 'ribosomal protein', 'trna synthetase', 'translation elongation factor', 'transcription', 'replication', 'repair', 'recombinase', 'helicase', 'topoisomerase', 'integrase', 'ligase']):
        return "Replication/Repair/Transcription/Translation"
    elif any(keyword in product_description for keyword in ['chaperone', 'heat shock', 'cold shock', 'oxidative stress', 'antibiotic resistance', 'stress response', 'detoxification', 'peroxidase', 'catalase', 'superoxide dismutase']):
        return "Stress Response"
    elif any(keyword in product_description for keyword in ['hydrolase', 'transferase', 'lyase', 'ligase', 'synthase', 'mutase', 'epimerase', 'racemase', 'protease', 'esterase', 'peptidase', 'isomerase', 'nuclease', 'glycosylase']):
        return "Other Enzymes"
    elif any(keyword in product_description for keyword in ['hypothetical protein', 'uncharacterized protein', 'putative protein', 'domain of unknown function', 'orf']):
        return "Hypothetical/Unknown"
    else:
        return "Miscellaneous"

# --- Core analysis function ---
def analyze_ferredoxin_neighborhoods(target_functional_category: str):
    print(f"\n--- Analyzing ferredoxins related to: '{target_functional_category}' ---")

    # Load necessary files (assuming they are in the same directory as the script)
    genbank_file_path = 'Moorella_thermoacetica_ATCC_39073_complete_genome.gb'
    protein_features_path = 'Moorella_protein_features.tsv'
    ferredoxin_fes_path = 'Moorella_ferredoxin_FeS_annotated.tsv'

    # --- 1. Robust GenBank file parsing ---
    genome_sequence = ""
    gene_features_parsed_list = [] # List to store parsed feature dictionaries
    current_feature_data = None
    in_features_block = False
    in_origin_block = False

    try:
        with open(genbank_file_path, 'r') as f:
            lines = f.readlines()

        line_idx = 0
        while line_idx < len(lines):
            line = lines[line_idx]

            if line.startswith("FEATURES"):
                in_features_block = True
                in_origin_block = False
                line_idx += 1
                continue
            elif line.strip().startswith("ORIGIN"):
                if current_feature_data and current_feature_data.get('type') == 'CDS':
                    gene_features_parsed_list.append(current_feature_data)
                current_feature_data = None
                in_features_block = False
                in_origin_block = True
                line_idx += 1
                continue
            elif line.startswith("//"):
                if current_feature_data and current_feature_data.get('type') == 'CDS':
                    gene_features_parsed_list.append(current_feature_data)
                break # End of record
            
            if in_origin_block:
                genome_sequence += "".join(re.findall(r'[a-zA-Z]', line)).upper()
                line_idx += 1
                continue
            
            if in_features_block:
                feature_line_match = re.match(r'^\s{5,}(\S+)\s+(.+)', line)
                
                if feature_line_match: # New feature line
                    if current_feature_data and current_feature_data.get('type') == 'CDS':
                        gene_features_parsed_list.append(current_feature_data)
                    
                    feature_type = feature_line_match.group(1)
                    location_raw = feature_line_match.group(2).strip()

                    current_feature_data = {'type': feature_type, 'qualifiers': {}, 'location_raw': location_raw}
                    
                    if feature_type == 'CDS':
                        coords = re.findall(r'(\d+)\.\.(\d+)', location_raw)
                        if coords:
                            current_feature_data['start'] = int(coords[0][0])
                            current_feature_data['end'] = int(coords[0][1])
                            current_feature_data['strand'] = -1 if "complement" in location_raw else 1
                        else:
                            current_feature_data['start'] = 0
                            current_feature_data['end'] = 0
                            current_feature_data['strand'] = 1
                    else: # Non-CDS feature types, initialize their location data to defaults
                        current_feature_data['start'] = 0
                        current_feature_data['end'] = 0
                        current_feature_data['strand'] = 1

                    line_idx += 1
                    continue
                    
                elif re.match(r'^\s{21,}/', line): # Qualifier line
                    if current_feature_data is None:
                        line_idx += 1
                        continue

                    qualifier_match = re.match(r'\s+/(.+?)=(.+)', line)
                    if qualifier_match:
                        key = qualifier_match.group(1)
                        value = qualifier_match.group(2).strip().strip('"')
                        
                        temp_line_idx = line_idx + 1
                        while temp_line_idx < len(lines) and \
                              not re.match(r'^\s{5,}\S', lines[temp_line_idx]) and \
                              not re.match(r'^\s{21,}/', lines[temp_line_idx]) and \
                              not lines[temp_line_idx].strip().startswith("ORIGIN") and \
                              not lines[temp_line_idx].strip().startswith("//"):
                            
                            value += lines[temp_line_idx].strip().strip('"')
                            temp_line_idx += 1
                        
                        current_feature_data['qualifiers'][key] = value
                        line_idx = temp_line_idx
                        continue
                    
                    line_idx += 1
                    continue
            
            line_idx += 1
            continue # Ensure line_idx always advances

        # After loop, add the very last feature if it was being processed
        if current_feature_data and current_feature_data.get('type') == 'CDS':
            gene_features_parsed_list.append(current_feature_data)

    except Exception as e:
        print(f"Error parsing GenBank file: {e}")
        return # Exit if parsing fails

    # Convert parsed features to DataFrame
    df_genome_features = pd.DataFrame([
        {
            'locus_tag': f['qualifiers'].get('locus_tag', ''),
            'start': f.get('start', 0),
            'end': f.get('end', 0),
            'strand': f.get('strand', 1),
            'product': f['qualifiers'].get('product', ''),
            'protein_id': f['qualifiers'].get('protein_id', ''),
            'translation': f['qualifiers'].get('translation', '')
        } for f in gene_features_parsed_list if f.get('type') == 'CDS'
    ])

    df_genome_features_cleaned = df_genome_features[df_genome_features['locus_tag'] != ''].copy()
    df_genome_features_cleaned['start'] = df_genome_features_cleaned['start'].astype(int)
    df_genome_features_cleaned['end'] = df_genome_features_cleaned['end'].astype(int)

    print(f"GenBank parsing complete. Total genome length: {len(genome_sequence)} bp, CDS features parsed: {len(df_genome_features_cleaned)}")
    print(df_genome_features_cleaned.head())

    # --- 2. Categorize all protein features and identify ferredoxins ---
    df_protein_features = pd.read_csv(protein_features_path, sep='\t')
    df_ferredoxin_features = pd.read_csv(ferredoxin_fes_path, sep='\t')

    df_protein_features['functional_category_refined'] = df_protein_features['product'].apply(categorize_product_full_redefinition)

    ferredoxin_locus_tags_from_file = set(df_ferredoxin_features['locus_tag'])
    df_protein_features['is_ferredoxin'] = df_protein_features['locus_tag'].isin(ferredoxin_locus_tags_from_file)
    
    # Sort proteins by locus_tag to maintain genomic order for neighborhood analysis
    df_protein_features['locus_numeric'] = df_protein_features['locus_tag'].apply(
        lambda x: int(x.split('_')[1]) if '_' in x and x.split('_')[1].isdigit() else 0
    )
    df_protein_features_sorted = df_protein_features.sort_values(by='locus_numeric').reset_index(drop=True)

    # --- 3. Analyze ferredoxin neighborhoods ---
    window_size_genes = 5 
    ferredoxin_neighborhood_data_refined = []

    for index, row in df_protein_features_sorted.iterrows():
        if row['is_ferredoxin']:
            ferredoxin_locus_tag = row['locus_tag']
            ferredoxin_product = row['product']
            ferredoxin_category = row['functional_category_refined']

            start_index = max(0, index - window_size_genes)
            end_index = min(len(df_protein_features_sorted) - 1, index + window_size_genes)

            neighborhood_slice = df_protein_features_sorted.iloc[start_index : end_index + 1]

            neighborhood_info = {
                'ferredoxin_locus_tag': ferredoxin_locus_tag,
                'ferredoxin_product': ferredoxin_product,
                'ferredoxin_category': ferredoxin_category,
                'neighborhood_genes': [n_row['locus_tag'] for _, n_row in neighborhood_slice.iterrows()],
                'neighborhood_functional_categories_refined': [n_row['functional_category_refined'] for _, n_row in neighborhood_slice.iterrows()],
                'neighborhood_products': [n_row['product'] for _, n_row in neighborhood_slice.iterrows()]
            }
            ferredoxin_neighborhood_data_refined.append(neighborhood_info)

    df_ferredoxin_neighborhoods_refined_full_set = pd.DataFrame(ferredoxin_neighborhood_data_refined)

    # --- 4. Filter based on user's target functional category ---
    # Ferredoxins whose neighborhood contains the target_functional_category
    df_filtered_by_target_category = df_ferredoxin_neighborhoods_refined_full_set[
        df_ferredoxin_neighborhoods_refined_full_set['neighborhood_functional_categories_refined'].apply(
            lambda categories: target_functional_category in categories
        )
    ].copy()

    if df_filtered_by_target_category.empty:
        print(f"No ferredoxins found with '{target_functional_category}' in their genomic neighborhood.")
        return # Exit if no relevant Fds found

    # Add 'representative_related_pathway' column (picking the target category if present)
    # The priority order for this column can be simplified as we're already filtering by a specific category.
    # If the target category is found, that's the representative. Otherwise, it would be 'None' from previous logic.
    df_filtered_by_target_category['representative_related_pathway'] = target_functional_category

    # --- 5. Output TSV file for the filtered ferredoxins ---
    sanitized_category_name = re.sub(r"[^a-zA-Z0-9_]", "_", target_functional_category)
    output_tsv_path = f'Moorella_Fd_related_to_{sanitized_category_name}.tsv'

    output_columns_tsv = [
        'ferredoxin_locus_tag',
        'ferredoxin_product',
        'ferredoxin_category',
        'representative_related_pathway' # New column
        # You can add 'neighborhood_genes', 'neighborhood_functional_categories_refined', 'neighborhood_products'
        # if you want the full neighborhood context in the TSV
    ]
    df_filtered_by_target_category[output_columns_tsv].to_csv(output_tsv_path, sep='\t', index=False)
    print(f"\nSuccessfully generated TSV file for ferredoxins related to '{target_functional_category}': {output_tsv_path}")

    # --- 6. Generate GenBank files for the filtered ferredoxins' neighborhoods ---
    output_gb_folder = f"GenBank_Neighborhoods_for_{sanitized_category_name}"
    os.makedirs(output_gb_folder, exist_ok=True)
    
    # Create a mapping from locus_tag to its index in the sorted full genome DataFrame
    df_genome_features_cleaned['locus_numeric'] = df_genome_features_cleaned['locus_tag'].apply(
        lambda x: int(x.split('_')[1]) if '_' in x and x.split('_')[1].isdigit() else 0
    )
    df_genome_features_sorted_by_locus = df_genome_features_cleaned.sort_values(by='locus_numeric').reset_index(drop=True)
    locus_tag_to_sorted_index = {tag: i for i, tag in enumerate(df_genome_features_sorted_by_locus['locus_tag'])}

    num_gb_files_generated = 0
    for fd_locus_tag in df_filtered_by_target_category['ferredoxin_locus_tag'].unique():
        fd_index_in_sorted = locus_tag_to_sorted_index.get(fd_locus_tag)
        
        if fd_index_in_sorted is None:
            print(f"Warning: Ferredoxin {fd_locus_tag} not found in the parsed genome features. Skipping GenBank file generation.")
            continue

        fd_row_from_genome_df = df_genome_features_sorted_by_locus.iloc[fd_index_in_sorted]

        neighborhood_start_idx_linear = max(0, fd_index_in_sorted - window_size_genes)
        neighborhood_end_idx_linear = min(len(df_genome_features_sorted_by_locus) - 1, fd_index_in_sorted + window_size_genes)

        neighborhood_genes_df_for_gb = df_genome_features_sorted_by_locus.iloc[
            neighborhood_start_idx_linear : neighborhood_end_idx_linear + 1
        ].copy()

        fragment_start_abs = neighborhood_genes_df_for_gb['start'].min()
        fragment_end_abs = neighborhood_genes_df_for_gb['end'].max()

        # Extract the DNA sequence fragment (1-based GenBank coords to 0-based Python slicing)
        extracted_seq = genome_sequence[fragment_start_abs - 1 : fragment_end_abs]
        
        features_for_gb = []
        for _, gene_row in neighborhood_genes_df_for_gb.iterrows():
            relative_start = gene_row['start'] - fragment_start_abs + 1
            relative_end = gene_row['end'] - fragment_start_abs + 1

            feature_str = f'     CDS             '
            if gene_row['strand'] == -1:
                feature_str += f'complement({relative_start}..{relative_end})'
            else:
                feature_str += f'{relative_start}..{relative_end}'

            # Add qualifiers for SnapGene visibility
            # Combine locus_tag and product for product qualifier as requested for display in SnapGene
            display_product = f"{gene_row['locus_tag']} ({gene_row['product']})" if gene_row['product'] else gene_row['locus_tag']
            feature_str += f'\n                     /locus_tag="{gene_row["locus_tag"]}"' # Keep standard locus_tag
            feature_str += f'\n                     /product="{display_product}"' # Use product for combined display
            if gene_row['protein_id']:
                feature_str += f'\n                     /protein_id="{gene_row["protein_id"]}"'
            if gene_row['translation']:
                translation_lines = [gene_row['translation'][i:i+60] for i in range(0, len(gene_row['translation']), 60)]
                feature_str += f'\n                     /translation="{translation_lines[0]}'
                for line_part in translation_lines[1:]:
                    feature_str += f'\n                                 {line_part}'
                feature_str += '"'
            features_for_gb.append(feature_str)
            
        genbank_content = f"LOCUS       {fd_locus_tag}_neighborhood {len(extracted_seq):>10} bp    DNA     linear BCT\n"
        genbank_content += f"DEFINITION  Genomic neighborhood of {fd_locus_tag} related to {target_functional_category} in Moorella thermoacetica ATCC 39073.\n"
        genbank_content += "ACCESSION   CUSTOM_NEIGHBORHOOD\n"
        genbank_content += "VERSION     CUSTOM_NEIGHBORHOOD.1\n"
        genbank_content += "KEYWORDS    {target_functional_category}; Ferredoxin; Moorella thermoacetica.\n"
        genbank_content += "SOURCE      Moorella thermoacetica ATCC 39073\n"
        genbank_content += "  ORGANISM  Moorella thermoacetica ATCC 39073\n"
        genbank_content += "            Bacteria; Firmicutes; Clostridia; Thermoanaerobacterales;\n"
        genbank_content += "            Thermoanaerobacteraceae; Moorella group; Moorella.\n"
        genbank_content += "FEATURES             Location/Qualifiers\n"

        for feature_str in features_for_gb:
            genbank_content += feature_str + "\n"

        genbank_content += "ORIGIN\n"
        for i in range(0, len(extracted_seq), 60):
            sub_seq = extracted_seq[i : i+60]
            formatted_line = f"{i+1:>9} " + " ".join([sub_seq[j:j+10] for j in range(0, len(sub_seq), 10)])
            genbank_content += formatted_line + "\n"
        genbank_content += "//\n"

        output_file_name_gb = os.path.join(output_gb_folder, f"{fd_locus_tag}_neighborhood.gb")
        with open(output_file_name_gb, 'w') as out_f:
            out_f.write(genbank_content)
        num_gb_files_generated += 1

    print(f"Successfully generated {num_gb_files_generated} GenBank files in the '{output_gb_folder}' folder.")
    print(f"Each file contains the specified ferredoxin gene and its {window_size_genes} upstream and {window_size_genes} downstream gene neighbors, with annotation information visible in SnapGene (locus_tag and product combined in /product).")


# --- User Input and Execution ---
if __name__ == "__main__":
    # Get all possible refined categories for user reference
    # This requires running a dummy categorization or knowing the categories beforehand
    # For now, let's list some key categories for the user to choose from:
    all_possible_categories = [
        "Electron Bifurcation - HydABC",
        "Electron Bifurcation - NfnAB",
        "ATP_synthesis - Ech Complex",
        "Wood-Ljungdahl Pathway - Carbonyl Branch (CODH/ACS)",
        "Wood-Ljungdahl Pathway - Carbonyl Branch (Formate Dehydrogenase)",
        "Wood-Ljungdahl Pathway - Methyl Branch",
        "Key Electron Donor - PFOR",
        "Ferredoxin", # Ferredoxin itself
        "Carbohydrate Metabolism",
        "Amino Acid Metabolism",
        "Nucleotide Metabolism",
        "Lipid Metabolism",
        "Cofactor Metabolism",
        "Energy Metabolism/Electron Transport",
        "Transport",
        "Regulatory/Signaling",
        "Cell Structure/Motility",
        "Replication/Repair/Transcription/Translation",
        "Stress Response",
        "Other Enzymes",
        "Hypothetical/Unknown",
        "Miscellaneous"
    ]
    
    print("Available functional categories for search:")
    for cat in all_possible_categories:
        print(f"- {cat}")
    
    user_input_category = input("\nEnter the exact functional category you want to search for: ").strip()

    if user_input_category not in all_possible_categories:
        print(f"Error: '{user_input_category}' is not a recognized functional category. Please choose from the list above.")
    else:
        analyze_ferredoxin_neighborhoods(user_input_category)