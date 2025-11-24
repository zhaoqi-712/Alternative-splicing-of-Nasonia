"""
Honeybee Differential Splicing Gene Filter

Description:
This script processes honeybee gene expression data to identify genes showing significant
differential splicing at specific developmental time points (25-40h and 55-70h) based on
a specified False Discovery Rate (FDR) threshold. The results are saved in CSV files
with the FDR threshold embedded in the filename for easy reference.

Key Features:
1. Reads wide-format gene expression data with FDR values
2. Takes FDR threshold as a parameter in the main function
3. Filters genes meeting significance criteria for each time point
4. Saves results with FDR value in filename for traceability
5. Handles missing data by defaulting to non-significant values (FDR=1)

Input Requirements:
- CSV file containing gene IDs and FDR values for time points
- First row should contain headers including 'gene_id', '25-40h', and '55-70h'
- Empty FDR values are treated as non-significant (FDR=1)

Output:
- Two CSV files (one per time point) containing filtered gene lists
- Filenames include time point and FDR threshold (e.g., 'genes_25-40h_FDR_0p05.csv')

Usage:
1. Call main() function with desired FDR threshold parameter
2. Script will generate output files in the specified directory
"""

import csv
import pandas as pd
from typing import List, Dict


def read_csv_to_dicts(file_path: str) -> List[Dict]:
    """
    Read a CSV file and store each row as a dictionary with headers as keys.

    Args:
        file_path (str): Path to the CSV file.

    Returns:
        List[Dict]: A list of dictionaries representing each row of data.
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


def filter_genes_by_fdr(gene_data: List[Dict], time_point: str, fdr_threshold: float) -> List[str]:
    """
    Filter genes based on FDR values for a specific time point.

    Args:
        gene_data (List[Dict]): List of gene data dictionaries
        time_point (str): The time point column to check (e.g., '25-40h')
        fdr_threshold (float): FDR threshold value

    Returns:
        List[str]: List of gene IDs that meet the FDR threshold
    """
    filtered_genes = []
    for data in gene_data:
        try:
            fdr_value = float(data.get(time_point, 1))  # Default to 1 if missing
            if fdr_value < fdr_threshold:
                filtered_genes.append(data['gene_id'])
        except ValueError:
            continue
    return list(set(filtered_genes))  # Remove duplicates


def save_filtered_genes(gene_list: List[str], save_path: str, time_point: str, fdr_threshold: float) -> None:
    """
    Save filtered gene list to CSV with FDR threshold in filename.

    Args:
        gene_list (List[str]): List of gene IDs to save
        save_path (str): Directory to save the file
        time_point (str): Time point identifier
        fdr_threshold (float): FDR threshold used for filtering
    """
    # Create filename with FDR threshold
    filename = f"genes_{time_point}_FDR_{fdr_threshold:.2f}.csv"  # Replace . with p for filename
    full_path = f"{save_path}/{filename}" if not save_path.endswith('/') else f"{save_path}{filename}"

    # Save to CSV with proper formatting
    df = pd.DataFrame({"gene_id": gene_list})
    df.to_csv(full_path, index=False, encoding='utf-8')
    print(f"Saved {len(gene_list)} genes to {filename}")


def main(fdr_threshold: float = 0.05) -> None:
    """
    Main function to process and filter gene data.

    Args:
        fdr_threshold (float): The FDR threshold to use for filtering (default: 0.05)
    """
    # Configuration
    INPUT_FILE = r'C:/Users/15611/Desktop/test_script/2-4_Change formate/Output/Honeybee_data_all_FDR_wide.csv'
    SAVE_DIR = "C:/Users/15611/Desktop/test_script/2-5_Apis significant gene/Output"
    #fdr_threshold = 0.05

    # Validate FDR threshold
    if not 0 < fdr_threshold < 1:
        raise ValueError("FDR threshold must be between 0 and 1")

    # Read gene data
    gene_data = read_csv_to_dicts(INPUT_FILE)
    if not gene_data:
        raise ValueError("No data was read from the input file")

    # Process empty values (set to 1 = no significance)
    for data in gene_data:
        for key, value in data.items():
            if value == "":
                data[key] = "1"

    # Time points to analyze
    time_points = ['25-40h', '55-70h']

    # Filter and save genes for each time point
    for time_point in time_points:
        filtered_genes = filter_genes_by_fdr(gene_data, time_point, fdr_threshold)
        save_filtered_genes(filtered_genes, SAVE_DIR, time_point, fdr_threshold)


if __name__ == "__main__":
    # Example usage with different FDR thresholds:
    # main(fdr_threshold=0.05)  # For 5% FDR
    # main(fdr_threshold=0.01)  # For 1% FDR

    # Default run with 0.05 FDR
    main()