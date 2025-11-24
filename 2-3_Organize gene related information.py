"""
Gene Annotation Script

This script performs the following steps:
1. Reads splicing event data from a CSV file (long format, no FDR filtering).
2. Converts the raw data into a list of dictionaries with basic gene event information.
3. Loads GO and InterPro annotation data from a separate CSV file.
4. Matches and appends the GO annotations to the corresponding genes.
5. Replaces empty lists or placeholder values with None for clarity.
6. Exports the final annotated gene information to a CSV file.
"""

import csv


def read_csv_to_list(file_path):
    """
    Read a CSV file and return the contents as a list of tuples.

    Args:
        file_path (str): Path to the input CSV file.

    Returns:
        list: A list of tuples, each representing a row in the CSV file.
    """
    data = []
    with open(file_path, mode='r', encoding='utf-8') as file:
        reader = csv.reader(file)
        next(reader)  # Skip header
        for row in reader:
            data.append(tuple(row))
    return data


def basic_info(gene_data):
    """
    Convert a list of tuples into a list of dictionaries with basic gene information.

    Args:
        gene_data (list): List of tuples containing gene information.

    Returns:
        list: List of dictionaries with keys: gene_id, time, event, fdr, exons_of_event.
    """
    gene_info = []
    for gene_id, time, event, fdr, exons_of_event in gene_data:
        dic = {
            'gene_id': gene_id,
            'time': time,
            'event': event,
            'fdr': fdr,
            'exons_of_event': exons_of_event
        }
        gene_info.append(dic)
    return gene_info


def get_GO_terms(file_path_GO):
    """
    Parse GO annotation file and return a dictionary mapping gene IDs to GO terms.

    Args:
        file_path_GO (str): Path to the GO annotation CSV file.

    Returns:
        dict: Dictionary mapping gene_id to GO and InterPro terms.
    """
    gene_data_GO = read_csv_to_list(file_path_GO)
    GO_dic = {}
    for gene_id, Gene_Name, GO_BP, GO_CC, GO_MF, InterPro in gene_data_GO:
        GO_BP = ["GO:" + item for item in GO_BP[3:-1].split(",GO:")]
        GO_CC = ["GO:" + item for item in GO_CC[3:-1].split(",GO:")]
        GO_MF = ["GO:" + item for item in GO_MF[3:-1].split(",GO:")]
        InterPro = InterPro[:-1].split(",")

        GO_dic[gene_id] = {
            'Gene_Name': Gene_Name,
            'GOTERM_BP_DIRECT': GO_BP,
            'GOTERM_CC_DIRECT': GO_CC,
            'GOTERM_MF_DIRECT': GO_MF,
            'INTERPRO': InterPro
        }
    return GO_dic


def delete_empty(gene_info):
    """
    Replace empty list or placeholder values in gene info dictionaries with None.

    Args:
        gene_info (list): List of dictionaries containing gene information.

    Returns:
        list: Cleaned list with empty values replaced by None.
    """
    for entry in gene_info:
        for key, value in entry.items():
            if value in ([''], [], ['GO:']):
                entry[key] = None
    return gene_info


def save_gene_info_to_csv(gene_info, output_file):
    """
    Save gene information (as list of dictionaries) to a CSV file.

    Args:
        gene_info (list): List of dictionaries containing gene-related information.
        output_file (str): Path to save the output CSV file.
    """
    headers = list(gene_info[0].keys())
    with open(output_file, mode='w', newline='', encoding='utf-8') as file:
        writer = csv.writer(file)
        writer.writerow(headers)
        for gene in gene_info:
            row = []
            for key in headers:
                value = gene.get(key, None)
                if isinstance(value, list):
                    value = "; ".join(str(v) for v in value)
                row.append(value)
            writer.writerow(row)
    print(f"Data successfully written to: {output_file}")





def main():
    # ========== MAIN EXECUTION SECTION ==========

    # Step 1: Load gene splicing data (the result of step 1, long format, no FDR filtering)
    input_gene_file = r"C:/Users/15611/Desktop/test_script/2-3_Organize gene related information/Input/Nasonia_all_FDR.csv"
    gene_data = read_csv_to_list(input_gene_file)

    # Step 2: Convert to list of dictionaries
    gene_info = basic_info(gene_data)

    # Step 3: Load GO annotations
    input_GO_file = r"C:/Users/15611/Desktop/test_script/2-3_Organize gene related information/Input/GO_terms_annotation.csv"
    GO_dic = get_GO_terms(input_GO_file)

    # Step 4: Annotate each gene_info item with GO terms
    for entry in gene_info:
        gene_id = entry['gene_id']
        if gene_id in GO_dic:
            for key, value in GO_dic[gene_id].items():
                entry[key] = value

    # Step 5: Replace empty values with None
    gene_info = delete_empty(gene_info)

    # Step 6: Save final annotated data to CSV
    output_gene_file = r"C:/Users/15611/Desktop/test_script/2-3_Organize gene related information/Output/Gene_info_result_all.csv"
    save_gene_info_to_csv(gene_info, output_gene_file)


if __name__ == "__main__":
    main()