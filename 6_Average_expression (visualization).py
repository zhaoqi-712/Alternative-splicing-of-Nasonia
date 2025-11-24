"""
GTF File Processing Pipeline

Description:
This script processes GTF files containing transcript expression data (TPM values) and performs:
1. GTF file parsing and attribute extraction
2. Transcript-to-gene ID mapping
3. Expression analysis and averaging across replicates
4. Target gene filtering based on minimum TPM thresholds
5. Output of processed data to CSV files

Key Features:
- Processes multiple GTF files from different timepoints and replicates
- Calculates average TPM values across biological replicates
- Filters for specific target genes of interest
- Generates separate output files for male and female samples

Input Requirements:
- GTF files with TPM values in the attributes field
- Transcript-to-gene ID mapping file (CSV format)
- List of target gene IDs to analyze

Output:
- CSV files containing average TPM values for target genes
- Separate files for male and female samples
- Files named by timepoint and sample type

Usage:
1. Configure input paths and target gene list
2. Run script to process all specified timepoints
3. Check output directory for results
"""

import csv
import os
from typing import Dict, List, Optional, Set


def csv_to_dict(file_path: str) -> Dict[str, str]:
    """
    Read a two-column CSV file into a dictionary.

    Args:
        file_path: Path to input CSV file (transcript_id, gene_id)

    Returns:
        Dictionary mapping transcript IDs to gene IDs
    """
    data_dict = {}
    with open(file_path, mode='r', encoding='utf-8') as file:
        reader = csv.reader(file)
        next(reader)  # Skip header
        data_dict = {row[0]: row[1] for row in reader if len(row) >= 2}
    return data_dict


def parse_attributes(attribute_string: str) -> Dict[str, str]:
    """
    Parse GTF attribute field into key-value pairs.

    Args:
        attribute_string: The 9th field from GTF file

    Returns:
        Dictionary of attribute key-value pairs
    """
    attributes = {}
    for attr in attribute_string.strip().strip(";").split('; '):
        key_value = attr.strip().split(' ', 1)
        if len(key_value) == 2:
            key, value = key_value
            attributes[key] = value.strip('"').strip(";")
    return attributes


def parse_gtf_line(line: str, feature: str = "transcript") -> Optional[Dict]:
    """
    Parse a single line from GTF file.

    Args:
        line: Input line from GTF file
        feature: Feature type to filter for (default: "transcript")

    Returns:
        Parsed record as dictionary or None if line doesn't match feature
    """
    if line.startswith("#"):
        return None

    fields = line.strip().split('\t')
    if len(fields) < 9 or fields[2] != feature:
        return None

    attributes = parse_attributes(fields[8])
    return {
        "seqname": fields[0],
        "source": fields[1],
        "feature": fields[2],
        "transcript_id": attributes.get("transcript_id", ""),
        "TPM": float(attributes.get("TPM", "0").strip(";")),
        "start": int(fields[3]),
        "end": int(fields[4]),
        "strand": fields[6],
        "attributes": attributes
    }


def parse_gtf(gtf_file: str, feature: str = "transcript", TPM_min: float = 0.0) -> List[Dict]:
    """
    Parse GTF file and filter records by feature type and minimum TPM.

    Args:
        gtf_file: Path to GTF file
        feature: Feature type to extract (default: "transcript")
        TPM_min: Minimum TPM threshold (default: 0.0)

    Returns:
        List of parsed records meeting criteria
    """
    with open(gtf_file, 'r', encoding='utf-8') as file:
        return [
            record for line in file
            if (record := parse_gtf_line(line, feature)) and record["TPM"] >= TPM_min
        ]


def add_gene_id(gtf_records: List[Dict], trans_dict: Dict[str, str]) -> List[Dict]:
    """
    Add gene ID to GTF records using transcript-to-gene mapping.

    Args:
        gtf_records: List of parsed GTF records
        trans_dict: Transcript-to-gene ID mapping dictionary

    Returns:
        Updated records with gene IDs added
    """
    for record in gtf_records:
        record["gene_id"] = trans_dict.get(record["transcript_id"], record["transcript_id"])
    return gtf_records


def process_gtf(file_path: str, trans_dict: Dict[str, str], TPM_min: float = 0.1) -> List[Dict]:
    """
    Process GTF file with transcript-to-gene mapping and TPM filtering.

    Args:
        file_path: Path to GTF file
        trans_dict: Transcript-to-gene ID mapping
        TPM_min: Minimum TPM threshold (default: 0.1)

    Returns:
        Processed GTF records with gene IDs
    """
    return add_gene_id(parse_gtf(file_path, TPM_min=TPM_min), trans_dict)


def generate_file_paths(prefix: str, timepoint: str, base_path: str) -> Dict[str, str]:
    """
    Generate file paths for all sample replicates.

    Args:
        prefix: Sample prefix (e.g., "3" for timepoint 3)
        timepoint: Timepoint identifier (e.g., "8h")
        base_path: Base directory path

    Returns:
        Dictionary mapping sample names to file paths
    """
    folder = f"{prefix}-{timepoint}_Assembly_-e"
    return {
        f"{sex}{prefix}{rep}": os.path.join(base_path, folder, f"Nv{sex}_{prefix}{rep}-sorted.gtf")
        for sex in "FM" for rep in "ABC"
    }


def analyze_expression(
        gtf_data: Dict[str, Dict[str, Dict]],
        sample_names: List[str],
        trans_dict: Dict[str, str],
        target_gene_list: List[str],
        min_tpm: float = 0.01
) -> Dict[str, Dict]:
    """
    Analyze expression data and calculate average TPM for target genes.

    Args:
        gtf_data: Dictionary of processed GTF data
        sample_names: List of sample names to analyze
        trans_dict: Transcript-to-gene ID mapping
        target_gene_list: List of target gene IDs
        min_tpm: Minimum average TPM threshold (default: 0.01)

    Returns:
        Dictionary of filtered results with average TPM values
    """
    # Collect all unique transcript IDs
    unique_ids: Set[str] = set()
    for sample_data in gtf_data.values():
        if isinstance(sample_data, dict):
            unique_ids.update(sample_data.keys())

    # Calculate expression values
    expression_dict = {}
    for transcript_id in unique_ids:
        tpm_values = [
            gtf_data[sample].get(transcript_id, {}).get("TPM", 0)
            for sample in sample_names
            if sample in gtf_data and isinstance(gtf_data[sample], dict)
        ]
        avg_tpm = sum(tpm_values) / len(sample_names) if sample_names else 0

        # Store if matches target gene and meets TPM threshold
        if transcript_id in trans_dict:
            gene_id = trans_dict[transcript_id]
            if gene_id in target_gene_list and avg_tpm > min_tpm:
                expression_dict[transcript_id] = {
                    "gene_ID": gene_id,
                    "transcript_ID": transcript_id,
                    "TPM_ave": avg_tpm
                }

    return expression_dict


def write_dict_to_csv(output_path: str, data_dict: Dict[str, Dict]) -> None:
    """
    Write dictionary data to CSV file.

    Args:
        output_path: Path to output CSV file
        data_dict: Dictionary containing data to write
    """
    with open(output_path, 'w', newline='', encoding='utf-8') as f:
        writer = csv.writer(f)
        writer.writerow(["Key", "gene_ID", "transcript_ID", "TPM_ave"])
        for transcript_id, data in data_dict.items():
            writer.writerow([
                transcript_id,
                data["gene_ID"],
                transcript_id,
                data["TPM_ave"]
            ])


def process_timepoint(prefix: str, timepoint: str) -> None:
    """
    Process all data for a single timepoint.

    Args:
        prefix: Sample prefix (e.g., "3")
        timepoint: Timepoint identifier (e.g., "8h")
    """
    # Configuration
    # The relationship between the transcript ID and the gene ID
    dic_path = r"C:/Users/15611/Desktop/test_script/6_Average_expression (visualization)/Input/Dictionary_transcript_entrez_mapping.csv"
    trans_dict = csv_to_dict(dic_path)

    # The genes of interest
    target_genes = [
        '100678718', '100123363', '100122481', '100302336',
        '100680400', '100123628', '100117702', '100116002',
        '100123232', '100121285', '100678824', '100117117',
        '100678462', '103316989', '100117261', '100116264'
    ]

    # The path of the expression data
    base_path = r"C:/Users/15611/Desktop/test_script/6_Average_expression (visualization)/Input/8h+12h_assembly with -e option"

    output_path = r"C:/Users/15611/Desktop/test_script/6_Average_expression (visualization)/Output"

    # Generate and process file paths (for the specific time point)
    paths = generate_file_paths(prefix, timepoint, base_path)

    gtf_data = {
        key: {rec["transcript_id"]: rec for rec in process_gtf(path, trans_dict)}
        for key, path in paths.items()
    }

    # Analyze expression
    female_samples = [f"F{prefix}{rep}" for rep in "ABC"]
    male_samples = [f"M{prefix}{rep}" for rep in "ABC"]

    female_results = analyze_expression(gtf_data, female_samples, trans_dict, target_genes)
    male_results = analyze_expression(gtf_data, male_samples, trans_dict, target_genes)

    # Write output files
    output_dir = os.path.join(output_path, f"{prefix}-{timepoint}_Assembly_-e")
    os.makedirs(output_dir, exist_ok=True)

    female_output = os.path.join(output_dir, f"{timepoint}_average_TPM_female.csv")
    male_output = os.path.join(output_dir, f"{timepoint}_average_TPM_male.csv")

    write_dict_to_csv(female_output, female_results)
    write_dict_to_csv(male_output, male_results)

    print(f"Processed {timepoint}:")
    print(f" Female samples: {len(female_results)} transcripts saved to {female_output}")
    print(f" Male samples: {len(male_results)} transcripts saved to {male_output}")


if __name__ == "__main__":
    # Process multiple timepoints
    timepoints = [("3", "8h"), ("4", "12h")]
    for prefix, timepoint in timepoints:
        process_timepoint(prefix, timepoint)