# We could also use female-exons annotations as reference
# Actually this annotation is the merged annotation of female (included the predicted transcripts),
# thus, it is not good enough. We did not use it in the scoring.

import csv
from collections import defaultdict
from Bio import SeqIO
from Bio.Seq import Seq


def parse_gff(file_path):
    """
    Parse a GFF (or GTF) file and extract exon regions for each gene, including strand information.
    Ensures that duplicated exons are not counted more than once.

    Args:
        file_path (str): Path to the GFF or GTF file.

    Returns:
        dict: A dictionary with gene IDs as keys and lists of exon tuples as values.
              Each exon tuple contains (exon_ID, chromosome, start, end, strand).
    """
    gene_exons = defaultdict(list)
    seen_exons = set()  # Track already processed exons to avoid duplicates

    with open(file_path) as gff:
        for line in gff:
            if line.startswith("#"):
                continue  # Skip comment lines

            fields = line.strip().split("\t")
            if fields[2] == "exon":  # Only process exon entries
                if fields[1] == "StringTie":
                    chrom = fields[0]
                    start = int(fields[3]) - 1  # Convert from 1-based to 0-based
                    end = int(fields[4])
                    strand = fields[6]
                    attributes = fields[8]

                    gene_id = None
                    exon_ID = None

                    # Extract gene_id and exon_number from attributes
                    for attribute in attributes.split("; "):
                        key_value = attribute.split(" ")
                        if len(key_value) == 2:
                            key, value = key_value
                            if key == "gene_name":
                                gene_id = value.strip('"')
                            elif key == "exon_number":
                                exon_ID = value.strip('"')

                    # Add exon to dictionary if gene_id and exon_ID are valid and not yet seen
                    if gene_id and exon_ID:
                        exon_key = (gene_id, start, end, strand)
                        if exon_key not in seen_exons:
                            seen_exons.add(exon_key)
                            gene_exons[gene_id].append((exon_ID, chrom, start, end, strand))

    return gene_exons


def count_sequence_in_exons(genome_file, gff_file, target_seqs, output_csv):
    """
    Count occurrences of specific target sequences in each exon, considering gene strand direction.

    Args:
        genome_file (str): Path to the reference genome FASTA file.
        gff_file (str): Path to the GFF or GTF annotation file.
        target_seqs (list): List of target sequences to search for.
        output_csv (str): Path to write the output CSV containing exon-level counts.

    Returns:
        dict: A dictionary of gene ID to list of tuples (chrom, strand, match count per exon).
    """
    # Load genome sequences as a dictionary
    genome = SeqIO.to_dict(SeqIO.parse(genome_file, "fasta"))

    # Parse exon annotations from GFF/GTF
    gene_exons = parse_gff(gff_file)

    exon_match_counts = defaultdict(list)

    with open(output_csv, mode="w", newline="") as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(
            ["Gene ID", "Exon ID", "Chromosome", "Exon Start", "Exon End", "Strand", "Target Sequence Count"]
        )

        for gene_id, exons in gene_exons.items():
            for exon_ID, chrom, start, end, strand in exons:
                exon_seq = str(genome[chrom].seq[start:end]).upper()
                count = 0

                # Count each target sequence depending on strand orientation
                for target_seq in target_seqs:
                    if strand == "+":
                        count += exon_seq.count(target_seq.upper())
                    elif strand == "-":
                        reverse_complement_seq = str(Seq(target_seq).reverse_complement()).upper()
                        count += exon_seq.count(reverse_complement_seq)

                # Write exon-level data (convert start to 1-based for output)
                csvwriter.writerow([gene_id, exon_ID, chrom, start + 1, end, strand, count])
                exon_match_counts[gene_id].append((chrom, strand, count))

    return exon_match_counts


def count_max_match_per_gene(exon_match_counts, output_csv):
    """
    For each gene, determine the exon with the highest number of target sequence matches,
    and write the result to a CSV file.

    Args:
        exon_match_counts (dict): Dictionary from count_sequence_in_exons().
        output_csv (str): Path to write the summary CSV.
    """
    with open(output_csv, mode="w", newline="") as csvfile:
        csvwriter = csv.writer(csvfile)
        csvwriter.writerow(["Gene ID", "Chromosome", "Strand", "Target Sequence Count (max)"])

        for gene_id, matches in exon_match_counts.items():
            # Identify the exon with the maximum match count
            max_match = max(matches, key=lambda x: x[2])
            chrom, strand, max_count = max_match
            csvwriter.writerow([gene_id, chrom, strand, max_count])

def main():
    # ========== MAIN EXECUTION SECTION ==========

    # Input file paths
    genome_file = r"C:/Users/15611/Desktop/test_script/Genome and annotations/Nv_genome.fna"
    gff_file = r"E:/Research_Practice/data/Female_exon_annotation/merged_transcriptom_females_modified.gtf"

    # List of target sequences to count
    target_seqs = [
        "CGAAGATA", "CCCTGAAGATTTGC", "GGAAGATA", "GGAAGATC",
        "TGAAGATC", "CGAAGATA", "TGAAGATT"
    ]

    # Output file paths
    output_csv_exons = r"C:/Users/15611/Desktop/test_script/3_Tra binding sites/Output/exon_counts_with_strand_female.csv"
    output_csv_max = r"C:/Users/15611/Desktop/test_script/3_Tra binding sites/Output/max_match_per_gene_female.csv"

    # Run processing functions
    exon_match_counts = count_sequence_in_exons(genome_file, gff_file, target_seqs, output_csv_exons)
    count_max_match_per_gene(exon_match_counts, output_csv_max)

    # Print completion messages
    print(f"Results for exons written to {output_csv_exons}")
    print(f"Max match counts per gene written to {output_csv_max}")


if __name__ == "__main__":
    main()
