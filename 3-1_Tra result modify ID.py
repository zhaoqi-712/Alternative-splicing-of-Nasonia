"""
Gene ID Processor for NCBI Entrez

This script processes a CSV file containing gene identifiers (e.g., LOC-prefixed or gene symbols).
It performs the following:
1. If the gene ID starts with "LOC", the prefix is removed.
2. Otherwise, the gene name is queried against NCBI Entrez to retrieve the corresponding Entrez Gene ID
   for the species *Nasonia vitripennis*.
3. The processed gene IDs are saved into a new output CSV file.

Dependencies:
- Biopython (for Entrez access)
- Internet connection (to query NCBI)

Author: Zhaoqi Leng
Date: 2025-1-22
"""

import csv
from Bio import Entrez

# Required by NCBI Entrez API
Entrez.email = "zhaoqi.leng@wur.nl"  # Replace with your actual email

def get_gene_id_from_ncbi(gene_name, species_name="Nasonia vitripennis"):
    """
    Query NCBI Entrez to get the Entrez Gene ID for a given gene name and species.

    Args:
        gene_name (str): The name of the gene (e.g., "hexamerin").
        species_name (str): The full name of the species (default is "Nasonia vitripennis").

    Returns:
        str or None: The first matching Entrez Gene ID if found; otherwise, None.
    """
    try:
        query = f"{gene_name}[Gene Name] AND {species_name}[Organism]"
        search_result = Entrez.esearch(db="gene", term=query)
        record = Entrez.read(search_result)
        if int(record["Count"]) > 0:
            return record["IdList"][0]
        else:
            return None
    except Exception as e:
        print(f"Error querying gene {gene_name}: {e}")
        return None


def process_csv(input_file, output_file):
    """
    Read a CSV file, process gene IDs (removing 'LOC' or replacing with Entrez ID),
    and write the updated records to a new output file.

    Args:
        input_file (str): Path to the input CSV file.
        output_file (str): Path to the output CSV file.
    """
    try:
        with open(input_file, mode='r', encoding='utf-8') as infile, \
             open(output_file, mode='w', encoding='utf-8', newline='') as outfile:

            reader = csv.reader(infile)
            writer = csv.writer(outfile)

            # Read and write header row
            headers = next(reader)
            writer.writerow(headers)

            # Process each row
            for row in reader:
                gene_id = row[0]
                if gene_id.startswith("LOC"):
                    row[0] = gene_id.replace("LOC", "")
                if gene_id.startswith("gene-LOC"):
                    row[0] = gene_id.replace("gene-LOC", "")
                else:
                    entrez_id = get_gene_id_from_ncbi(gene_id)
                    row[0] = entrez_id if entrez_id else gene_id  # fallback to original if query fails
                writer.writerow(row)

        print(f"Processing completed. Results saved to {output_file}")

    except Exception as e:
        print(f"Error processing file: {e}")


# === Example Usage ===
input = r"C:/Users/15611/Desktop/test_script/3-1_Tra result modify ID/Input/max_match_per_gene.csv"
output = r"C:/Users/15611/Desktop/test_script/3-1_Tra result modify ID/Output/max_match_per_gene_ID_modified.csv"
process_csv(input, output)

input = r"C:/Users/15611/Desktop/test_script/3-1_Tra result modify ID/Input/exon_counts_with_strand.csv"
output = r"C:/Users/15611/Desktop/test_script/3-1_Tra result modify ID/Output/exon_counts_with_strand_ID_modified.csv"
process_csv(input, output)



