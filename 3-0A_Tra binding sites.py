"""
Exon Target Sequence Counter
This script aims to find the number of Tra binding sites in the Nv genome (per exon).
Also count the max number of Tra binding sites per exon in each gene.

This script performs the following tasks:
1. Parses a GFF annotation file to extract exon locations in the genome and strand information for each gene.
2. For each exon, counts the occurrences of given target sequences based on strand orientation.
   - If the gene is on the "+" strand, sequences are counted as-is.
   - If on the "-" strand, the reverse complement of the sequences is used.
3. Outputs a CSV file listing the number of matches per exon.
4. Generates a summary CSV file with the maximum match count per exon for each gene.

Dependencies:
- Biopython (for parsing FASTA files and working with sequences)

Author: Zhaoqi Leng
Date: 2025-01-22
"""

import csv
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq


def parse_gff(file_path):
    """
    Parse a GFF file to extract exon locations and strand orientation for each gene.

    Args:
        file_path (str): Path to the GFF annotation file.

    Returns:
        dict: A dictionary mapping gene_id to a list of tuples:
              (exon_ID, chromosome, start, end, strand)
    """
    gene_exons = defaultdict(list)
    seen_exons = set()
    with open(file_path) as gff:
        for line in gff:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if fields[2] == "exon":
                chrom = fields[0]
                start = int(fields[3]) - 1  # Convert from 1-based to 0-based
                end = int(fields[4])
                strand = fields[6]
                attributes = fields[8]

                gene_id = next((x.split("=")[1] for x in attributes.split(";") if "gene" in x), None)
                exon_ID = next((x.split("=")[1] for x in attributes.split(";") if "ID" in x), None)

                if gene_id and exon_ID:
                    exon_key = (gene_id, start, end, strand)
                    if exon_key not in seen_exons:
                        seen_exons.add(exon_key)
                        gene_exons[gene_id].append((exon_ID, chrom, start, end, strand))
    return gene_exons


def count_sequence_in_exons(genome_file, gff_file, target_seqs, output_csv):
    """
    Count occurrences of target sequences in each exon, respecting strand direction,
    and write the results to a CSV file.

    Args:
        genome_file (str): Path to the genome FASTA file.
        gff_file (str): Path to the GFF annotation file.
        target_seqs (list of str): List of target sequences to count.
        output_csv (str): Path to save the exon-level output CSV.

    Returns:
        dict: A dictionary mapping gene_id to a list of tuples:
              (chromosome, strand, total_count_in_exons)
    """
    genome = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))
    gene_exons = parse_gff(gff_file)

    exon_match_counts = defaultdict(list)
    with open(output_csv, mode="w", newline="") as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(
            ["Gene ID", "Exon ID", "Chromosome", "Exon Start", "Exon End", "Strand", "Target Sequence Count"])

        for gene_id, exons in gene_exons.items():
            for exon_ID, chrom, start, end, strand in exons:
                exon_seq_raw = str(genome[chrom].seq[start:end]).upper()
                count = 0
                if strand == "+":
                    exon_seq = exon_seq_raw

                elif strand == "-":
                    exon_seq = Seq(exon_seq_raw).reverse_complement().upper()

                for seq in target_seqs:
                    target_seq = seq.upper()
                    count += exon_seq.count(target_seq)

                """for target_seq in target_seqs:
                    #count += exon_seq.count(target_seq.upper())
                    if strand == "+":
                        count += exon_seq.count(target_seq.upper())
                    elif strand == "-":
                        reverse_complement_seq = str(Seq(target_seq).reverse_complement()).upper()
                        count += exon_seq.count(reverse_complement_seq)"""

                csvwriter.writerow([gene_id, exon_ID, chrom, start, end, strand, count,exon_seq])
                exon_match_counts[gene_id].append((chrom, strand, count))
    return exon_match_counts


def count_max_match_per_gene(exon_match_counts, output_csv):
    """
    For each gene, find the exon with the maximum number of target sequence matches
    and save the result to a CSV file.

    Args:
        exon_match_counts (dict): Output from count_sequence_in_exons().
        output_csv (str): Path to save the gene-level summary CSV.
    """
    with open(output_csv, mode="w", newline="") as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(["Gene ID", "Chromosome", "Strand", "Target Sequence Count (max)"])

        for gene_id, matches in exon_match_counts.items():
            max_match = max(matches, key=lambda x: x[2])
            chrom, strand, max_count = max_match
            csvwriter.writerow([gene_id, chrom, strand, max_count])



def main():
    # ====== FILE PATHS AND TARGET SEQUENCES ======
    # Uses all exons as reference
    # Path to the genome fasta file
    genome_file = r"C:/Users/15611/Desktop/test_script/Genome and annotations/Nv_genome.fna"
    # Path to the genome annotation gff file
    gff_file = r"C:/Users/15611/Desktop/test_script/Genome and annotations/Nv_annotations.gff"
    # test fru
    # gff_file = r"C:\Users\15611\Desktop\Nv_annotations.gff"
    target_seqs = ["CGAAGATA", "CCCTGAAGATTTGC", "GGAAGATA", "GGAAGATC",
                   "TGAAGATC", "CGAAGATA", "TGAAGATT"]
    # target_seqs = ["GAAGAT"]

    output_csv_exons = r"C:/Users/15611/Desktop/test_script/3_Tra binding sites/Output/try_exon_counts_with_strand.csv"
    output_csv_max = r"C:/Users/15611/Desktop/test_script/3_Tra binding sites/Output/try_max_match_per_gene.csv"

    # ====== RUN ANALYSIS ======

    exon_match_counts = count_sequence_in_exons(genome_file, gff_file, target_seqs, output_csv_exons)
    count_max_match_per_gene(exon_match_counts, output_csv_max)

    print(f"Results for exons written to {output_csv_exons}")
    print(f"Max match counts per gene written to {output_csv_max}")

if __name__ == "__main__":
    main()