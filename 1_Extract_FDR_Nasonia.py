"""
Script to extract splicing events with FDR < 0.05 for multiple time points
for the species Nasonia vitripennis. Gene names are queried via NCBI using Biopython's Entrez module.

Inputs: the JC results of rMATS of Nasonia data
Output: a CSV file, containing all splicing events
The gene names are ALL replaced with Entrez ID
"""

from Bio import Entrez
import csv

def get_gene_id_from_ncbi(gene_name, species_name="Nasonia vitripennis"):
    """
    Query NCBI Entrez to retrieve the gene ID for a given gene name and species.
    Because in the original results, the genes are annotated as either Entrez ID like "LOCxxxxxxx",
    or the gene name like “DSX”. This function is to transfer all the names into Entrez ID.

    Args:
        gene_name (str): The gene name to query.
        species_name (str): The species to filter the gene search (default is Nasonia vitripennis).

    Returns:
        str or None: The first matching gene ID if found, otherwise None.
    """
    try:
        query = f"{gene_name}[Gene Name] AND {species_name}[Organism]"
        search_result = Entrez.esearch(db="gene", term=query)
        record = Entrez.read(search_result)
        if int(record["Count"]) > 0:
            return record["IdList"][0]  # Return the first matching gene ID
        else:
            return None
    except Exception as e:
        print(f"Error querying gene {gene_name}: {e}")
        return None


def update_geneid_fdr(input_path, time, event, geneid_fdr_list):
    """
    Process a file to extract genes with FDR < 0.05 and append their IDs to a list.

    Args:
        time (str): The time point (e.g., '1_1h').
        event (str): The splicing event type (e.g., 'SE', 'RI').
        geneid_fdr_list (list): A list to append results as (time, event, geneid, fdr) tuples.

    Returns:
        list: The updated geneid_fdr_list including new qualifying gene entries.
    """
    file_path = generate_file_path(input_path, time, event)
    with open(file_path, mode='r', encoding='utf-8') as file:
        reader = csv.reader(file, delimiter='\t')
        header = next(reader)

        # Get the index positions of 'GeneID' and 'FDR' columns
        geneid_index = header.index('GeneID')
        fdr_index = header.index('FDR')

        for row in reader:
            geneid = row[geneid_index]
            fdr = row[fdr_index]
            fdr_float = float(fdr)

            # Convert gene name to NCBI gene ID if not already in LOC format
            if not geneid.startswith("LOC"):
                converted_id = get_gene_id_from_ncbi(geneid)
                geneid = converted_id if converted_id else geneid

            # Remove 'LOC' prefix if present
            if geneid.startswith("LOC"):
                geneid = geneid[3:]

            # Append to result list
            geneid_fdr_list.append((time, event, geneid, fdr))

    return geneid_fdr_list


def generate_file_path(input_path, time, event):
    """
    Construct the file path of the certain rMATS result based on the given time and event.

    Args:
        time (str): The time point.
        event (str): The event type.

    Returns:
        str: The constructed file path to the input data file.
    """
    file_path = f"{input_path}/results_{time}/{event}.MATS.JC.txt"
    return file_path


def extract_all(input_path, gene_list):
    """
    Iterate through all time points and event types, collecting genes with FDR < 0.05.

    Args:
        gene_list (list): List to store the (time, event, geneid, fdr) results.

    Returns:
        list: The fully populated gene_list with all qualifying entries.
    """
    time_all = ["1_1h", "2_4h", "3_8h", "4_12h", "5_17h", "6_22h", "7_27h"]
    event_all = ["SE", "A5SS", "A3SS", "MXE", "RI"]
    for time in time_all:
        for event in event_all:
            update_gene_list = update_geneid_fdr(input_path, time, event, gene_list)
    return update_gene_list


def main():
    """Main execution function"""
    # Set email for NCBI Entrez (required)
    Entrez.email = "zhaoqi.leng@wur.nl"  # Replace with your own email address

    # Define input file path, the root path of the JC results
    input_path = r"C:/Users/15611/Desktop/test_script/1_Extract_FDR_Nasonia/Input"

    # Define output file path
    output_file = 'C:/Users/15611/Desktop/test_script/1_Extract_FDR_Nasonia/Output/gene_data.csv'

    # Initialize empty list to store gene IDs and FDR values
    geneid_fdr_list = []
    gene_all = extract_all(input_path, geneid_fdr_list)
    print(gene_all)

    # Write output to CSV
    with open(output_file, mode='w', newline='', encoding='utf-8') as file:
        writer = csv.writer(file)

        # Write header
        writer.writerow(['time', 'event', 'geneid', 'fdr'])

        # Write each row of data
        for row in gene_all:
            row_with_str_geneid = (row[0], row[1], str(row[2]), row[3])
            writer.writerow(row_with_str_geneid)

    print(f"Data successfully written to {output_file}")

if __name__ == "__main__":
    main()