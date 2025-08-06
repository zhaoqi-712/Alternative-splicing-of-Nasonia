"""
Script to extract all the gene ID
Input: A CSV file with all events, the gene ID is on the first column
Output: A CSV gene list with all unique gene
"""
import csv
from openpyxl import Workbook
from collections import Counter

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

def write_to_excel(unique_gene_ids, file_path):
    """
    Write a list of unique gene IDs to an Excel file (one per row).

    Args:
        unique_gene_ids (list): List of gene IDs.
        file_path (str): Path to save the Excel file.
    """
    wb = Workbook()
    ws = wb.active
    ws.title = "Gene IDs"
    for row_index, gene_id in enumerate(unique_gene_ids, start=1):
        ws.cell(row=row_index, column=1, value=gene_id)
    wb.save(file_path)
    print(f"Data successfully written to: {file_path}")

def extract_unique_gene_ids(data):
    """
    Extract unique gene IDs from a list of tuples (assuming gene_id is the 3rd item).

    Args:
        data (list): List of tuples.

    Returns:
        list: Sorted list of unique gene IDs.
    """
    gene_ids = [row[0] for row in data]
    return sorted(set(gene_ids))


# ========== MAIN EXECUTION SECTION ==========
def main():
    # Step 1: Load gene splicing data (long format, no FDR filtering)
    # Input data is the output of script (1)
    input_gene_file = r"C:/Users/15611/Desktop/test_script/2-1_Extract all gene names/Input/Nasonia_all_FDR.csv"
    gene_data = read_csv_to_list(input_gene_file)

    # Export unique gene IDs to Excel
    unique_gene_ids = extract_unique_gene_ids(gene_data)
    write_to_excel(unique_gene_ids,
                   r"C:/Users/15611/Desktop/test_script/2-1_Extract all gene names/Output/gene_unique_list.xlsx")

if __name__ == "__main__":
    main()

