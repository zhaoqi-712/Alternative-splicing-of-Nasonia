"""
Description:
    This script reads a tab-delimited file containing orthologous gene group information
    for multiple species. It filters the data to find orthologous gene pairs between
    *Nasonia vitripennis* and *Apis mellifera*, constructs a mapping dictionary between
    their gene IDs, and exports the results to a CSV file.

Steps:
    1. Read and filter the ortholog group data.
    2. Identify groups that contain both species.
    3. Create a gene ID mapping from Nasonia to Apis.
    4. Export the mapping as a CSV file.

Author:
    Zhaoqi Leng

Date:
    2025-01-22
"""

import csv

def filter_and_map_orthologs(file_path, output_path):
    """
    Reads a tab-delimited file of ortholog groups, filters to include only pairs
    of genes from Nasonia vitripennis and Apis mellifera, and exports the mappings
    to a CSV file.

    Args:
        file_path (str): Path to the input tab-delimited file.
        output_path (str): Path to the output CSV file.

    Returns:
        dict: A dictionary mapping Nasonia vitripennis gene IDs to Apis mellifera gene IDs.
    """
    ortholog_dict = {}

    # Read the file and collect ortholog groups for the two target species
    ortholog_groups = {}
    with open(file_path, 'r', encoding='utf-8') as file:
        reader = csv.reader(file, delimiter='\t')  # Assuming the file is tab-delimited
        for row in reader:
            if len(row) < 5:
                continue  # Skip malformed rows
            species = row[4]
            gene_id = row[2]
            ortholog_group = row[1]

            if species in ["Nasonia vitripennis", "Apis mellifera"]:
                if ortholog_group not in ortholog_groups:
                    ortholog_groups[ortholog_group] = {}
                ortholog_groups[ortholog_group][species] = gene_id

    # Build a mapping for ortholog groups that contain both species
    for ortholog_group, genes in ortholog_groups.items():
        if "Nasonia vitripennis" in genes and "Apis mellifera" in genes:
            ortholog_dict[genes["Nasonia vitripennis"]] = genes["Apis mellifera"]

    # Write the result to a CSV file
    with open(output_path, 'w', newline='', encoding='utf-8') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(["Nasonia vitripennis Gene", "Apis mellifera Gene"])
        for nv_gene, am_gene in ortholog_dict.items():
            writer.writerow([nv_gene, am_gene])

    return ortholog_dict


def main():
    # ========== MAIN EXECUTION SECTION ==========

    # Replace with your actual file paths
    file_path = r"C:/Users/15611/Desktop/test_script/4_Apis ortholog gene finding/Input/Hymenoptera_HGD_Ortho.tab"
    output_path = r"C:/Users/15611/Desktop/test_script/4_Apis ortholog gene finding/Output/ortholog_mapping.csv"

    # Run the mapping function
    ortholog_mapping = filter_and_map_orthologs(file_path, output_path)

    # Print result confirmation
    print(f"Mapping results have been written to: {output_path}")

if __name__ == "__main__":
    main()