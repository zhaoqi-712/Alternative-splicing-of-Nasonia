"""
Gene Data Processing Pipeline

This script processes gene data from multiple sources, performs various analyses including:
- FDR threshold checks (FDR < x in n time points)
- Gene Ontology term matching (3 GO terms and structural terms)
- Motif counting (Target sequences of Tra)
- Ortholog identification between Nasonia and Apis (Whether this gene has ortholog in Aphis,
  significantly identified as sexual spliced gene in embryonic stage)
- Expression data processing
- Scoring calculations

The final output is a comprehensive CSV file with all processed data.
"""

import csv
import pandas as pd

# ======================================================================
# CONSTANTS AND CONFIGURATION
# ======================================================================

# File paths for all input and output data
# 1. The file with FDR values and GO terms of each event.
# The result of step 2-4
DATA_PATH = r"C:/Users/15611/Desktop/test_script/2-4_Change formate/Output/Gene_info_result_wide.csv"

# 2. The file of max match (of Tra target motif per exon) of each gene
# The result of step 3-0A + 3-1
COUNTING_DATA_PATH = r"C:/Users/15611/Desktop/test_script/3-1_Tra result modify ID/Output/max_match_per_gene_ID_modified.csv"

# 3. The file containing the maximum matches of Tra target motifs per exon for each gene, specifically for female exons.
# The result of step 3-0B + 3-1, actually we did not use this data in scoring,
# because the female-exons we used as the reference was not correct
FEMALE_COUNTING_DATA_PATH = r"C:/Users/15611/Desktop/test_script/3-1_Tra result modify ID/Output/max_match_per_gene_female_ID_modified.csv"

# 4. Apis FDR document of 25-40h and 55-70h. Gene list with FDR < 0.05.
# The result of step 2-5.
APIS_2540_PATH = r"C:/Users/15611/Desktop/test_script/5_Create scoring matrix/Input/genes_25-40h_FDR_0.05.csv"
APIS_5570_PATH = r"C:/Users/15611/Desktop/test_script/5_Create scoring matrix/Input/genes_55-70h_FDR_0.05.csv"

# 5. The relationship between Apis and Nasonia orthologs
# The result of step 4
ORTHO_PATH = r"C:/Users/15611/Desktop/test_script/4_Apis ortholog gene finding/Output/ortholog_mapping.csv"

# 6. Gene expression data used for validation
# The treatment information of the expression data
TREATMENT_PATH = r"C:/Users/15611/Desktop/test_script/5_Create scoring matrix/Input/expression data_treatments.xlsx"
# The gene expression data of males and females in different time stages
GENE_EXPRESSION_PATH = r"C:/Users/15611/Desktop/test_script/5_Create scoring matrix/Input/expression data_ID_transferred_genome annotation.xlsx"

# Out put path
OUTPUT_PATH = r"C:/Users/15611/Desktop/test_script/5_Create scoring matrix/Output/Scoring_matrix_without cell.csv"

# Time point keys to check for FDR analysis (columns in the input data)
KEYS_TO_CHECK = ['3_8h', '4_12h', '5_17h', '6_22h']

# Fields used for calculating the total score with their weights
COUNTING_SCORE_FIELDS = [
    'number_FDR_<_0.01_(8-22h)',
    "number_FDR_punish",
    'GOTERM_BP_regulation: regulation of transcription, regulation of gene expression, regulation of DNA-templated transcription',
    'GOTERM_MF_DNA/RNA_binding: RNA binding, DNA-binding, DNA binding, mRNA binding',
    'GOTERM_CC_nucleus: GO:0005634~nucleus',
    'INTERPRO_ZnF: Znf_C2H2, DM_DNA-bd, DMRT',
    'GOTERM_BP_punishment: protein phosphorylation, transmembrane transport',
    'GOTERM_CC_punishment: cytoplasm, membrane, plasma membrane, cytosol, mitochondrion',
    'Target Sequence Score (female exons)',
    "Apis_ortholog_24-50h_FDR<0.05",
    "Apis_ortholog_55-70h_FDR<0.05"
]

# Configuration for Gene Ontology term processing
# Each entry specifies:
# - key: The field name in the input data to check
# - keywords: List of terms to search for in that field, check whether it's a part of the term
#   The keywords can be replaced to the terms of interest.
# - new_key: Name of the new field to create
# - func: Type of operation ('inclusion' adds positive score, 'punish' adds negative)
GO_CONFIGS = [
    # Positive scoring terms for Biological Process
    {
        'key': "GOTERM_BP_DIRECT",
        'keywords': ["regulation of transcription", "regulation of gene expression",
                     "regulation of DNA-templated transcription"],
        'new_key': "GOTERM_BP_regulation",
        'func': 'inclusion'
    },
    # Positive scoring terms for Molecular Function
    {
        'key': "GOTERM_MF_DIRECT",
        'keywords': ["RNA binding", "DNA-binding", "DNA binding", "mRNA binding"],
        'new_key': "GOTERM_MF_DNA/RNA_binding",
        'func': 'inclusion'
    },
    # Positive scoring terms for Cellular Component
    {
        'key': "GOTERM_CC_DIRECT",
        'keywords': ["GO:0005634~nucleus"],
        'new_key': "GOTERM_CC_nucleus",
        'func': 'inclusion'
    },
    # Positive scoring terms for InterPro domains
    {
        'key': "INTERPRO",
        'keywords': ["Znf_C2H2", "DM_DNA-bd", "DMRT"],
        'new_key': "INTERPRO_ZnF",
        'func': 'inclusion'
    },
    # Negative scoring terms for Biological Process
    {
        'key': "GOTERM_BP_DIRECT",
        'keywords': ["protein phosphorylation", "transmembrane transport"],
        'new_key': "GOTERM_BP_punishment",
        'func': 'punish'
    },
    # Negative scoring terms for Cellular Component
    {
        'key': "GOTERM_CC_DIRECT",
        'keywords': ["cytoplasm", "membrane", "plasma membrane", "cytosol", "mitochondrion"],
        'new_key': "GOTERM_CC_punishment",
        'func': 'punish'
    }
]


# ======================================================================
# UTILITY FUNCTIONS
# ======================================================================

def read_csv_to_dicts(file_path):
    """
    Read a CSV file and convert each row to a dictionary using headers as keys.

    Args:
        file_path (str): Path to the CSV file to read

    Returns:
        list: List of dictionaries where each dictionary represents a row
    """
    data = []
    try:
        with open(file_path, mode='r', encoding='utf-8') as file:
            reader = csv.DictReader(file)
            data = [row for row in reader]
    except FileNotFoundError:
        print(f"File not found: {file_path}")
    except Exception as e:
        print(f"Error reading file: {e}")
    return data


def excel_to_dict(file_path):
    """
    Read an Excel file and convert to dictionary with first column as keys and the rest as values..

    Args:
        file_path (str): Path to the Excel file

    Returns:
        dict: Dictionary where keys are from first column and values are lists of remaining values
    """
    # Read Excel without headers
    df = pd.read_excel(file_path, header=None)
    # Remove duplicate keys, keeping first occurrence
    df = df.drop_duplicates(subset=[0], keep='first')
    # Convert to dictionary with first column as keys
    return df.set_index(0).T.to_dict(orient='list')


def save_gene_info_to_csv(gene_info, output_file):
    """
    Save processed gene information to a CSV file.

    Args:
        gene_info (list): List of dictionaries containing gene data
        output_file (str): Path to save the output CSV file
    """
    # Get headers from the keys of the first dictionary
    headers = list(gene_info[0].keys()) if gene_info else []

    with open(output_file, mode='w', newline='', encoding='utf-8') as file:
        writer = csv.writer(file)
        # Write header row
        writer.writerow(headers)
        # Write each gene's data
        for gene in gene_info:
            row = [gene.get(key, None) for key in headers]
            writer.writerow(row)
    print(f"Data successfully written to file: {output_file}")


# ======================================================================
# DATA PROCESSING FUNCTIONS
# ======================================================================

def clean_and_count(data, threshold, keys_to_check):
    """
    Clean data and count values under specific keys that meet FDR threshold.

    Args:
        data (list): List of gene data dictionaries
        threshold (float): FDR threshold value to compare against
        keys_to_check (list): List of column names to check for FDR values

    Returns:
        list: Updated list with new count fields added
    """
    result = []
    for row in data:
        # Clean NA values by converting them to None
        cleaned_row = {k: (None if v == 'NA' else v) for k, v in row.items()}
        count = 0

        # Check each specified time point
        for key in keys_to_check:
            fdr_value = cleaned_row.get(key)
            if fdr_value is not None:
                try:
                    # Count if FDR value is below threshold
                    if float(fdr_value) < threshold:
                        count += 1
                except ValueError:
                    # Skip if value can't be converted to float
                    pass

        # Add count field to the row
        cleaned_row[f'count_FDR_<_{threshold}_(8-22h)'] = count
        result.append(cleaned_row)
    return result


def add_number_based_on_count(data, threshold):
    """
    Convert FDR counts (number of time points in 8-12h) to categorical scores (0, 1, or 2).
    Number of significant (FDR < threshold) time points --> scores:
    0 --> 0
    1 --> 1
    more than 1 --> 2

    Args:
        data (list): List of gene data dictionaries
        threshold (float): FDR threshold used to identify the count field

    Returns:
        list: Updated list with new score fields added
    """
    for row in data:
        # Construct field names based on threshold
        count_key = f"count_FDR_<_{threshold}_(8-22h)"
        number_key = f"number_FDR_<_{threshold}_(8-22h)"

        # Get count value, default to 0 if missing
        count_value = row.get(count_key, 0)

        try:
            count_value = int(count_value)
        except ValueError:
            count_value = 0

        # Convert count to categorical score
        if count_value == 0:
            number = 0
        elif count_value == 1:
            number = 1
        else:
            number = 2

        # Add score field to the row
        row[number_key] = number

    return data


def add_key_based_on_inclusion(data, key_to_check, keywords_list, new_key_name):
    """
    Add new field based on presence of keywords in specified field.
    Positive scoring - adds 0.5 if any keyword is found.
    This function is used to detect specific keywords in each kind of GO terms,
    and give 0.5 score if exists.

    Args:
        data (list): List of gene data dictionaries
        key_to_check (str): Field name to search in
        keywords_list (list): Keywords to search for
        new_key_name (str): Name for new field to add

    Returns:
        list: Updated list with new fields added
    """
    # Create descriptive field name listing all keywords
    key_name = f"{new_key_name}: {', '.join(keywords_list)}"

    for dic in data:
        # Get value to check (empty string if field missing)
        target_value = dic.get(key_to_check, "")
        # Set value to 0.5 if any keyword found, else 0
        dic[key_name] = 0.5 if any(keyword in target_value for keyword in keywords_list) else 0
    return data


def add_key_based_punish(data, key_to_check, keywords_list, new_key_name):
    """
    Add new field based on presence of keywords in specified field.
    Negative scoring - adds -100 if any keyword is found.

    Args:
        data (list): List of gene data dictionaries
        key_to_check (str): Field name to search in
        keywords_list (list): Keywords to search for
        new_key_name (str): Name for new field to add

    Returns:
        list: Updated list with new fields added
    """
    # Create descriptive field name listing all keywords
    key_name = f"{new_key_name}: {', '.join(keywords_list)}"

    for dic in data:
        # Get value to check (empty string if field missing)
        target_value = dic.get(key_to_check, "")
        # Set value to -100 if any keyword found, else 0
        dic[key_name] = -100 if any(keyword in target_value for keyword in keywords_list) else 0
    return data

def process_motif_data(event_info):
    """
    Process motif counting data from files (The result of step 3-1).
    The max number of Tra binding sites per exon of each gene.
    Scoring of max(hit/exon):
    0, 1 --> 0
    2, 3, 4 --> 1
    more than 4 --> 2

    Args:
        event_info (list): List of gene data dictionaries

    Returns:
        list: Updated list with motif count fields added
    """
    # Read and process main motif counting data
    counting_data = read_csv_to_dicts(COUNTING_DATA_PATH)
    counting_dict = {data['Gene ID']: data['Target Sequence Count (max matches per exon)']
                     for data in counting_data}

    # Read and process female-specific motif counting data
    female_counting_data = read_csv_to_dicts(FEMALE_COUNTING_DATA_PATH)
    female_counting_dict = {data['Gene ID']: data['Target Sequence Count (max)']
                            for data in female_counting_data}

    for dic in event_info:
        # Add main motif count for all exons
        count_motif = int(counting_dict.get(dic['gene_id'], 0))
        dic['Target Sequence Count (all exons) (max matches per exon)'] = count_motif

        # Process female-specific motif data if available
        if dic['gene_id'] in female_counting_dict:
            count_motif = int(female_counting_dict[dic['gene_id']])
            dic['Target Sequence Count (female exons) (max matches per exon)'] = count_motif

            # Convert count to categorical score
            if count_motif == 0 or count_motif == 1:
                value = 0
            elif 2 <= count_motif <= 4:
                value = 1
            elif count_motif >= 5:
                value = 2
            dic['Target Sequence Score (female exons)'] = value
        else:
            # Handle genes not expressed in females
            dic['Target Sequence Count (female exons) (max matches per exon)'] = "not expressed in female"
            dic['Target Sequence Score (female exons)'] = 0

        # Add FDR punishment score
        # move it to elsewhere
        count_FDR = int(dic.get('count_FDR_<_0.05_(8-22h)', 0))
        dic["number_FDR_punish"] = -100 if count_FDR == 0 else 0

    return event_info


def process_ortholog_data(event_info):
    """
    Process ortholog data from honeybee to identify conserved genes.

    Args:
        event_info (list): List of gene data dictionaries

    Returns:
        list: Updated list with ortholog fields added
    """
    # Read and process Apis (honeybee) data at 25-40h
    apis_2540 = [row['gene_id'] for row in read_csv_to_dicts(APIS_2540_PATH)]
    # Read and process Apis data at 55-70h
    apis_5570 = [row['gene_id'] for row in read_csv_to_dicts(APIS_5570_PATH)]

    # Create ortholog mapping between Apis and Nasonia genes
    ortho_info = {row['Apis mellifera Gene']: row['Nasonia vitripennis Gene']
                  for row in read_csv_to_dicts(ORTHO_PATH)}

    # Find Nasonia orthologs for significant Apis genes
    apis_2540_nasonia = [ortho_info[gene] for gene in apis_2540 if gene in ortho_info]
    apis_5570_nasonia = [ortho_info[gene] for gene in apis_5570 if gene in ortho_info]

    # Add ortholog information to each gene
    for dic in event_info:
        dic["Apis_ortholog_24-50h_FDR<0.05"] = 1 if dic['gene_id'] in apis_2540_nasonia else 0
        dic["Apis_ortholog_55-70h_FDR<0.05"] = 1 if dic['gene_id'] in apis_5570_nasonia else 0

    return event_info


def process_expression_data(event_info):
    """
    Process gene expression data.
    Gene expression of males and females in different stages.
    Using TPM value.

    Args:
        event_info (list): List of gene data dictionaries

    Returns:
        list: Updated list with expression fields added
    """
    # Read treatment metadata
    treatments = excel_to_dict(TREATMENT_PATH)
    # Read expression data
    expression_data = excel_to_dict(GENE_EXPRESSION_PATH)

    # Reorganize treatment data by averaging replicates
    treatments_new = {
        "Name": [],  # Treatment names
        "sexes": [],  # Sex information
        "stage": [],  # Developmental stage
        "RNAi": []  # RNAi treatment type
    }

    total_length = len(treatments['sexes'])
    exp_length = int(total_length / 5)  # 5 replicates per treatment

    # Prepare dictionary for averaged expression data
    expression_data_new = {str(gene): [] for gene in expression_data.keys()}

    # Process each treatment by averaging replicates
    for i in range(exp_length):
        # Create averaged treatment name
        name_new = treatments["Name"][i * 5][:-2] + "_ave"
        treatments_new["Name"].append(name_new)
        treatments_new["sexes"].append(treatments["sexes"][i * 5])
        treatments_new["stage"].append(treatments["stage"][i * 5])
        treatments_new["RNAi"].append(treatments["RNAi"][i * 5])

        # Calculate average expression across replicates for each gene
        for gene in expression_data:
            sum_exp = sum(expression_data[gene][i * 5 + j] for j in range(5))
            expression_data_new[str(gene)].append(sum_exp / 5)

    # Add expression data to each gene record
    for dic in event_info:
        gene_id = dic["gene_id"]
        if gene_id in expression_data_new:
            # Add expression values for GFP-RNAi treatments only
            for i, name_exp in enumerate(treatments_new["Name"]):
                if treatments_new["RNAi"][i] == "GFP-RNAi":
                    dic[name_exp] = expression_data_new[gene_id][i]

            # Calculate adult averages across tissues
            female_ave = (dic['Female_Abdomen_GFP-RNAi_ave'] +
                          dic['Female_Head_GFP-RNAi_ave'] +
                          dic['Female_Thorax_GFP-RNAi_ave']) / 3
            male_ave = (dic['Male_Abdomen_GFP-RNAi_ave'] +
                        dic['Male_Head_GFP-RNAi_ave'] +
                        dic['Male_Thorax_GFP-RNAi_ave']) / 3

            dic["Female_adult_ave"] = female_ave
            dic["Male_adult_ave"] = male_ave

            # Count how many female time points have expression >= 8 TPM
            female_exps = ['Female_1 day_GFP-RNAi_ave', 'Female_7 days_GFP-RNAi_ave',
                           'Female_9 days_GFP-RNAi_ave', 'Female_11 days_GFP-RNAi_ave',
                           'Female_adult_ave']
            dic["Female_exp_num"] = sum(1 for exp in female_exps if dic.get(exp, 0) >= 8)

            # Count how many male time points have expression >= 8 TPM
            male_exps = ['Male_1 day_GFP-RNAi_ave', 'Male_7 days_GFP-RNAi_ave',
                         'Male_9 days_GFP-RNAi_ave', 'Male_11 days_GFP-RNAi_ave',
                         'Male_adult_ave']
            dic["Male_exp_num"] = sum(1 for exp in male_exps if dic.get(exp, 0) >= 8)
        else:
            # Handle genes with no expression data
            for name_exp in treatments_new["Name"]:
                if treatments_new["RNAi"][treatments_new["Name"].index(name_exp)] == "GFP-RNAi":
                    dic[name_exp] = 0
            dic["Female_adult_ave"] = 0
            dic["Male_adult_ave"] = 0
            dic["Female_exp_num"] = 0
            dic["Male_exp_num"] = 0

    return event_info


def calculate_total_scores(event_info):
    """
    Calculate total composite score for each gene by summing all scoring fields.

    Args:
        event_info (list): List of gene data dictionaries

    Returns:
        list: Updated list with total_score field added
    """
    for dic in event_info:
        # Sum all specified scoring fields to get total score
        total_score = sum(float(dic.get(name, 0)) for name in COUNTING_SCORE_FIELDS)
        dic["total_score"] = total_score
    return event_info


# ======================================================================
# MAIN PROCESSING PIPELINE
# ======================================================================

def main():
    """
    Main execution function that runs the entire processing pipeline.
    """
    # 1. Read and process initial gene data
    print("Reading initial gene data...")
    event_info = read_csv_to_dicts(DATA_PATH)

    # 2. Perform FDR threshold analysis
    print("Performing FDR threshold analysis...")
    for threshold in [0.05, 0.01]:
        event_info = clean_and_count(event_info, threshold, KEYS_TO_CHECK)
        event_info = add_number_based_on_count(event_info, threshold)

    # 3. Process Gene Ontology terms
    print("Processing Gene Ontology terms...")
    for config in GO_CONFIGS:
        if config['func'] == 'inclusion':
            event_info = add_key_based_on_inclusion(
                event_info,
                config['key'],
                config['keywords'],
                config['new_key']
            )
        else:
            event_info = add_key_based_punish(
                event_info,
                config['key'],
                config['keywords'],
                config['new_key']
            )

    # 4. Process additional data sources
    print("Processing motif data...")
    event_info = process_motif_data(event_info)

    print("Processing ortholog data...")
    event_info = process_ortholog_data(event_info)

    print("Processing expression data...")
    event_info = process_expression_data(event_info)

    print("Calculating total scores...")
    event_info = calculate_total_scores(event_info)

    # 5. Save final results
    print("Saving results...")
    save_gene_info_to_csv(event_info, OUTPUT_PATH)

    print("Processing complete!")

if __name__ == "__main__":
    # Execute main function when run as a script
    main()