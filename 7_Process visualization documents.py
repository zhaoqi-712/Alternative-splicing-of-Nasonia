"""
This script processes RNA-seq data to:
1. Extract and filter transcript isoforms for target genes
2. Calculate and incorporate TPM expression values
3. Generate multiple output formats (TSV, GTF, FASTA) for visualization
4. Create Geneious-compatible annotation files

Key Features:
- Processes both male and female samples separately
- Maintains gene-transcript relationships from GTF annotations
- Preserves expression levels (TPM) from quantification data
- Generates genome sequences for target genes

Input Files Required:
1. GTF annotation file (Nv_annotations.gtf)
2. Merged transcriptome GTF (merged_transcriptom.gtf)
3. TPM expression tables (CSV format)
4. Transcript-to-gene mapping dictionary

Output Files Generated:
1. Filtered transcript tables (TSV)
2. Standard GTF files
3. Geneious-optimized GTF files (The annotation of the gene start from 1)
4. Gene sequences (FASTA)

Typical Workflow:
1. Preprocess GTF annotations
2. Filter transcripts for target genes
3. Incorporate expression data
4. Generate multiple output formats
5. Extract genomic sequences

Author: Zhaoqi Leng
Date: 2025-2
"""

import csv
import re
import os
from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

# Main Functions
def process_transcripts(processed_data, trans_ID_dict, gene, TPM_file, output_path, sex):
    """
    Main function to process GTF data, filter transcripts, and save outputs.
    Creates TSV, GTF, and FASTA files for the specified gene.
    """
    output_path_gene = f"{output_path}/{gene}"
    Path(output_path_gene).mkdir(parents=True, exist_ok=True)
    output_file = f"{output_path_gene}/{gene}_{sex}_transcripts.tsv"
    gtf_output_file = f"{output_path_gene}/{gene}_{sex}_transcripts.gtf"
    gtf_output_geneious = f"{output_path_gene}/{gene}_{sex}_geneious.gtf"
    fasta_output_file = f"{output_path_gene}/{gene}_{sex}_gene.fasta"

    transcript_list = write_from_gene(gene, TPM_file)
    if transcript_list:
        filtered_data = filter_transcripts(processed_data, transcript_list, trans_ID_dict)
        # Save in R-compatible format
        write_output(output_file, filtered_data)
        print(f"Processing complete. Output saved to {output_file}")

        # Save as GTF format
        write_gtf_output(gtf_output_file, filtered_data, TPM_file, gene)
        print(f"GTF file saved to {gtf_output_file}")

        # Save as Geneious-compatible GTF format
        write_gtf_output_modified(gtf_output_geneious, filtered_data, TPM_file, gene)
        print(f"Geneious GTF file saved to {gtf_output_file}")

        # Save FASTA
        gene_info = find_gene_info(gene)
        gtf_line = f"{gene_info[1]}\t{gene_info[2]}\t{gene_info[3]}\t{gene_info[4]}\t.\t{gene_info[6]}\t."
        extract_sequence_from_genome(gene, gene_info, fasta_output_file)
        print(f"Fasta file saved to {fasta_output_file}")
    else:
        print(f"{gene} has no transcripts")

# Core Processing Functions
def find_gene_info(gene_id):
    """
    Find gene information in GTF file for specified gene ID.
    Returns modified columns list if found, None otherwise.
    """
    gtf_file = "C:/Users/15611/Desktop/test_script/Genome and annotations/Nv_annotations.gtf"

    if not os.path.exists(gtf_file):
        raise FileNotFoundError(f"File not found: {gtf_file}")

    with open(gtf_file, 'r', encoding='utf-8') as file:
        for line in file:
            if line.startswith('#'):
                continue  # Skip comments

            columns = line.strip().split('\t')
            if len(columns) < 9:
                continue  # Ensure complete data

            if columns[2] == "gene" and f'GeneID:{gene_id}' in columns[8]:
                attributes = columns[8]
                new_attributes = []
                for attr in attributes.split('; '):
                    if 'gene_id' in attr or 'gene_biotype' in attr or 'gene' in attr:
                        new_attributes.append(attr.split(' ', 1)[1].replace('"', ''))
                    elif 'db_xref "GeneID:' in attr:
                        # Extract GeneID number
                        gene_id_from_xref = attr.split('GeneID:')[1].split('"')[0]
                        new_attributes.append(gene_id_from_xref)
                columns = columns[:8] + new_attributes
                return columns

    return None

def preprocess_gtf(file_path):
    """
    Preprocess GTF file to keep only exon/transcript data in standardized format.
    Returns list of processed data rows.
    Optional: keep CDS.
    """
    processed_data = []
    with open(file_path, 'r', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):
                continue  # Skip comments
            fields = line.strip().split('\t')
            if len(fields) < 9 or fields[2] not in ["exon", "transcript"]:
                continue  # Keep only exon/transcript types

            attributes = fields[8]
            transcript_match = re.search(r'transcript_id "(.*?)"', attributes)
            transcript_id = transcript_match.group(1) if transcript_match else "unknown"

            processed_data.append([
                fields[0],  # seqnames
                fields[3],  # start
                fields[4],  # end
                fields[6],  # strand
                fields[2],  # type
                "unknown",  # gene_name (placeholder)
                transcript_id  # transcript_name
            ])
    return processed_data

# File I/O Functions
def read_csv_to_dict(input_path):
    """
    Read CSV into dictionary structure with TPM average data.
    Returns dictionary in format: {key: {subkey: value, ...}, ...}
    """
    data_dict = {}
    with open(input_path, 'r', encoding='utf-8') as f:
        reader = csv.DictReader(f)
        for row in reader:
            key = row["Key"]
            data_dict[key] = {
                "gene_ID": row["gene_ID"],
                "transcript_ID": row["transcript_ID"],
                "TPM_ave": float(row["TPM_ave"])  # Store as float
            }
    return data_dict

def write_from_gene(gene, TPM_file):
    """
    Get all transcript isoforms for a given gene from TPM data.
    Returns list of transcript IDs or None if none found.
    """
    transcripts_isoforms = []
    for key, value in TPM_file.items():
        gene_name = value['gene_ID']
        if gene_name == gene:
            transcripts_isoforms.append(value['transcript_ID'])
    return transcripts_isoforms if transcripts_isoforms else None

def filter_transcripts(processed_data, transcript_list, trans_ID_dict):
    """
    Filter data by transcript IDs and replace gene_name using trans_ID_dict.
    Returns filtered data list.
    """
    filtered_data = []
    for row in processed_data:
        transcript_id = row[6]
        if transcript_id in transcript_list:  # Keep only specified transcripts
            gene_name = trans_ID_dict.get(transcript_id, "unknown")
            filtered_data.append(row[:5] + [gene_name] + row[6:])
    return filtered_data

def write_output(output_path, data):
    """
    Write filtered data to output file with header.
    """
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f, delimiter='\t')
        writer.writerow(["seqnames", "start", "end", "strand", "type", "gene_name", "transcript_name"])
        writer.writerows(data)

def csv_to_dict(file_path):
    """
    Read two-column CSV into dictionary (transcript_id -> gene_name).
    """
    data_dict = {}
    with open(file_path, mode='r', encoding='utf-8') as file:
        reader = csv.reader(file)
        next(reader)  # Skip header
        for row in reader:
            if len(row) >= 2:
                data_dict[row[0]] = row[1]  # transcript_id -> gene_name
    return data_dict

def write_gtf_output(output_path, data, TPM_file, gene_ID):
    """
    Write filtered data in standard GTF format.
    Includes TPM values in attributes.
    """
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        f.write("#seqnames\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tattributes\n")
        gene_info = find_gene_info(gene_ID)
        attributes0 = f'gene_id "{gene_info[9]}"; gene "{gene_info[9]}"; gene_biotype "{gene_info[11][:-1]}";'
        f.write(f"{gene_info[0]}\t{gene_info[1]}\t{gene_info[2]}\t{gene_info[3]}\t{gene_info[4]}\t.\t{gene_info[6]}\t.\t{attributes0}\n")

        for row in data:
            TPM_value = TPM_file[row[6]]['TPM_ave']
            attributes = f'gene_id "{row[5]}"; transcript_id "{row[6]}"; gene_name "{row[5]}"; TPM_ave "{TPM_value}"'
            f.write(f"{row[0]}\tStringTie\t{row[4]}\t{row[1]}\t{row[2]}\t.\t{row[3]}\t.\t{attributes}\n")

def write_gtf_output_modified(output_path, data, TPM_file, gene_ID):
    """
    Write filtered data in Geneious-compatible GTF format.
    Adjusts coordinates relative to gene start.
    """
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        f.write("#seqnames\tsource\tfeature\tstart\tend\tscore\tstrand\tframe\tattributes\n")
        gene_info = find_gene_info(gene_ID)
        start = int(gene_info[3])-1
        attributes0 = (f'gene_id "{gene_info[9]}"; gene "{gene_info[9]}"; gene_biotype "{gene_info[11][:-1]}"; '
                       f'real_start "{gene_info[3]}"; real_end "{gene_info[4]}";')
        f.write(f"{gene_info[9]}\t{gene_info[1]}\t{gene_info[2]}\t{int(gene_info[3]) - start}\t{int(gene_info[4]) - start}\t.\t{gene_info[6]}\t.\t{attributes0}\n")

        for row in data:
            TPM_value = TPM_file[row[6]]['TPM_ave']
            attributes = (f'gene_id "{row[5]}"; transcript_id "{row[6]}"; gene_name "{row[5]}"; TPM_ave "{TPM_value}"; '
                          f'real_start "{row[1]}"; real_end "{row[2]}";')
            f.write(f"{row[5]}\tStringTie\t{row[4]}\t{int(row[1]) - start}\t{int(row[2]) - start}\t.\t{row[3]}\t.\t{attributes}\n")

def extract_sequence_from_genome(gene_id, gtf_list, output_dir):
    """
    Extract gene sequence from genome file based on GTF annotation.
    Saves sequence in FASTA format.
    """
    genome_file = r"C:/Users/15611/Desktop/test_script/Genome and annotations/Nv_genome.fna"
    columns = gtf_list
    chrom = columns[0]
    start = int(columns[3])
    end = int(columns[4])
    strand = columns[6]

    genome_records = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
    if chrom not in genome_records:
        raise ValueError(f"Chromosome {chrom} not found in genome file.")

    seq = genome_records[chrom].seq[start - 1:end].upper()
    seq_record = SeqRecord(seq, id=gene_id, description=f"Gene {gene_id} sequence")
    os.makedirs(os.path.dirname(output_dir), exist_ok=True)
    SeqIO.write(seq_record, output_dir, "fasta")
    print(f"Sequence for gene {gene_id} extracted and saved to {output_dir}")

# Execution Code
if __name__ == "__main__":
    # Initialize dictionaries
    trans_ID_path = r"C:/Users/15611/Desktop/test_script/7_Process visualization documents/Input/Dictionary_transcript_entrez_mapping.csv"
    trans_ID_dict = csv_to_dict(trans_ID_path)
    output_path = r"C:/Users/15611/Desktop/test_script/7_Process visualization documents/Output"

    # Target genes to process
    target_gene_list = ['100678718', '100123363', '100122481', '100302336',
                        '100680400', '100123628', '100117702', '100116002',
                        '100123232', '100121285', '100678824', '100117117',
                        '100678462', '103316989', '100117261', '100116264']

    # Process male samples
    input_gtf_male = r"C:/Users/15611/Desktop/test_script/Genome and annotations/Nv_annotations.gtf"
    processed_data = preprocess_gtf(input_gtf_male)
    TPM_csv_path_male = r"C:/Users/15611/Desktop/test_script/7_Process visualization documents/Input/3-8h_Assembly_-e/8h_average_TPM_male.csv"
    TPM_csv_male = read_csv_to_dict(TPM_csv_path_male)
    for gene in target_gene_list:
        process_transcripts(processed_data, trans_ID_dict, gene, TPM_csv_male, output_path, "male")

    # Process female samples
    input_gtf_female = r"C:/Users/15611/Desktop/test_script/Genome and annotations/Nv_annotations.gtf"
    processed_data = preprocess_gtf(input_gtf_female)
    TPM_csv_path_female = r"C:/Users/15611/Desktop/test_script/7_Process visualization documents/Input/3-8h_Assembly_-e/8h_average_TPM_female.csv"
    TPM_csv_female = read_csv_to_dict(TPM_csv_path_female)
    for gene in target_gene_list:
        process_transcripts(processed_data, trans_ID_dict, gene, TPM_csv_female, output_path, "female")