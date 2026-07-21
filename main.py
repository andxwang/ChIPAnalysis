import argparse
import json
import pandas as pd
from gene_functions import annotate_peaks

GFF_HEADERS=['seqid', 'source', 'type', 'start', 'end', 'score', 'direction', 'phase', 'attr']


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description="Main script to annotate ChIP-seq peaks.")
    parser.add_argument("-c", "--config", default="config.json", help="Path to configuration JSON file (default: config.json)")
    parser.add_argument("-g", "--genes", help="Path to gene annotation GFF file")
    parser.add_argument("-p", "--peaks", help="Path to peaks GFF file")
    parser.add_argument("-o", "--output", help="Output CSV filename")

    args = parser.parse_args()
    
    config_path = args.config
    with open(config_path, 'r') as f:
        config = json.load(f)

    mab_gff_path = args.genes or config.get('paths', {}).get('genes')
    if not mab_gff_path:
        print("You didn't pass in a path to the gene Mab gff.")
        exit(1)
    try:
        mab_df = pd.read_table(mab_gff_path, header=None, names=GFF_HEADERS)
    except FileNotFoundError:
        print(f"Error: The file '{mab_gff_path}' was not found. Please check the path in the config.json file.")
        exit(1)
    print("Successfully loaded MAB GFF file with column names:", GFF_HEADERS)

    peaks_gff_path = args.peaks or config.get('paths', {}).get('peaks')
    if not peaks_gff_path:
        print("You didn't pass in a path to the peaks gff.")
        exit(1)
    try:
        peaks_df = pd.read_table(peaks_gff_path, header=None, names=GFF_HEADERS)
    except FileNotFoundError:
        print(f"Error: The file '{peaks_gff_path}' was not found. Please check the path in the config.json file.")
        exit(1)
    print("Successfully loaded peaks GFF file with column names:", GFF_HEADERS)

    print("Running analysis:")
    regulated_peaks = annotate_peaks(peaks_df, mab_df)
    print("Output preview:")
    print(regulated_peaks.head())
    
    out_path = args.output or config.get('paths', {}).get('output', 'regulated_peaks_out.csv')
    regulated_peaks.to_csv(out_path, index=False)
    print("Saved analysis to", out_path)
