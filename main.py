import argparse
import json
import pandas as pd
from gene_functions import annotate_peaks

GFF_HEADERS=['seqid', 'source', 'type', 'start', 'end', 'score', 'direction', 'phase', 'attr']


def load_gff_data(gff_path: str):
    if not gff_path:
        print("You didn't pass in a path to the gff file, or didn't input it in config.json.")
        exit(1)
    if not gff_path.endswith('.gff') and not gff_path.endswith('.txt'):
        print(f"Warning: {gff_path} does not have a .gff or .txt extension. Please check the file format.")
    try:
        df = pd.read_table(gff_path, header=None, names=GFF_HEADERS)
        print(f"Successfully loaded GFF file at '{gff_path}' with column names:", GFF_HEADERS)
        return df
    except FileNotFoundError:
        print(f"Error: The file '{gff_path}' was not found. Please check the path.")
        exit(1)


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
    mab_df = load_gff_data(mab_gff_path)

    peaks_gff_path = args.peaks or config.get('paths', {}).get('peaks')
    peaks_df = load_gff_data(peaks_gff_path)

    print("Running analysis:")
    regulated_peaks = annotate_peaks(peaks_df, mab_df)
    print("Output preview:")
    print(regulated_peaks.head())
    
    out_path = args.output or config.get('paths', {}).get('output_path', 'regulated_peaks_out.csv')
    regulated_peaks.to_csv(out_path, index=False)
    print("Saved analysis to", out_path)
