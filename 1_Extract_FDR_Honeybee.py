"""
Script to extract splicing events for multiple time points
for the species Apis mellifera (honeybee). Gene names are queried via NCBI using Biopython's Entrez module.

Inputs: JC result files from rMATS for honeybee data
Output: A CSV file containing all splicing events, along with exon-related information
"""

from Bio import Entrez
import csv

# Set your email address (required by NCBI Entrez)
Entrez.email = "zhaoqi.leng@wur.nl"  # Replace with your email


def get_gene_id_from_ncbi(gene_name, species_name="Apis mellifera"):
    """
    Query NCBI Gene database to retrieve a GeneID for a given gene name.

    Args:
        gene_name (str): The gene name or symbol.
        species_name (str): The scientific name of the organism (default: "Apis mellifera").

    Returns:
        str or None: The NCBI GeneID as a string if found, otherwise None.
    """
    try:
        query = f"{gene_name}[Gene Name] AND {species_name}[Organism]"
        search_result = Entrez.esearch(db="gene", term=query)
        record = Entrez.read(search_result)
        if int(record["Count"]) > 0:
            return record["IdList"][0]  # Return the first GeneID found
        else:
            return None
    except Exception as e:
        print(f"Error querying gene {gene_name}: {e}")
        return None


def update_geneid_fdr(time, event, geneid_fdr_list):
    """
    Process a single rMATS JC result file and extract entries with FDR < 0.05.
    Optionally convert gene names to GeneIDs via NCBI.

    Args:
        time (str): Developmental time point (e.g., "10-15h").
        event (str): Splicing event type (e.g., "SE", "MXE").
        geneid_fdr_list (list): A list to which filtered results will be appended.

    Returns:
        list: Updated geneid_fdr_list containing entries that passed the FDR filter.
    """
    file_path = generate_file_path(time, event)

    with open(file_path, mode='r', encoding='utf-8') as file:
        reader = csv.reader(file, delimiter='\t')
        header = next(reader)

        # Column index lookups
        geneid_index = header.index('GeneID')
        fdr_index = header.index('FDR')
        strand_index = header.index('strand')
        id_indices = [i for i, col in enumerate(header) if col == 'ID']
        if len(id_indices) < 2:
            raise ValueError("Second 'ID' column not found in header.")
        id1_index = id_indices[1]

        # Columns between 'strand' and second 'ID'
        between_columns = header[strand_index + 1: id1_index]

        # Process rows
        for row in reader:
            geneid = row[geneid_index]
            fdr = row[fdr_index]
            fdr_float = float(fdr)

            # Convert gene name to GeneID if necessary
            if not geneid.startswith("LOC"):
                converted_id = get_gene_id_from_ncbi(geneid)
                geneid = converted_id if converted_id else geneid

            if geneid.startswith("LOC"):
                geneid = geneid[3:]  # Remove LOC prefix

            between_data = {col: row[header.index(col)] for col in between_columns}

            geneid_fdr_list.append({
                'time': time,
                'event': event,
                'geneid': geneid,
                'fdr': fdr,
                'exon_data': between_data
            })

    return geneid_fdr_list


def generate_file_path(time, event):
    """
    Generate the file path for a given time point and splicing event.

    Args:
        time (str): Time point string (e.g., "25-40h").
        event (str): Splicing event type (e.g., "SE", "RI").

    Returns:
        str: File path to the corresponding JC result file.
    """
    file_path = f"C:/Users/15611/Desktop/test_script/1-1_Extract_FDR_honeybee/Input/DS_results_{time}/{event}.MATS.JC.txt"
    return file_path


def extract_all(gene_list):
    """
    Loop through all time points and event types to collect FDR.

    Args:
        gene_list (list): A list to collect all filtered gene entries.

    Returns:
        list: Final list of all gene entries that passed FDR filtering.
    """
    time_all = ["10-15h", "25-40h", "55-70h"]
    event_all = ["SE", "A5SS", "A3SS", "MXE", "RI"]

    for time in time_all:
        for event in event_all:
            update_gene_list = update_geneid_fdr(time, event, gene_list)
            print(f"{event} of {time} finished!")

    return update_gene_list


# Initialize and run extraction
geneid_fdr_list = []
gene_all = extract_all(geneid_fdr_list)

######### OUTPUT #########

# Output file path
output_file = r'C:/Users/15611/Desktop/test_script/1-1_Extract_FDR_honeybee/Output/Honeybee_all_FDR.csv'

# Write results to CSV
with open(output_file, mode='w', newline='', encoding='utf-8') as file:
    writer = csv.writer(file)

    writer.writerow(['gene_id', 'time', 'event', 'fdr', 'exon_data'])  # header

    for row in gene_all:
        exon_string = ''
        for key, value in row['exon_data'].items():
            exon_string += f"{key}:\t{value}\t;"

        base_data = [str(row['geneid']), row['time'], row['event'], row['fdr'], exon_string]
        writer.writerow(base_data)

print(f"Data successfully written to {output_file}")