import json
import pandas as pd
from gene_functions import annotate_peaks

GFF_HEADERS=['seqid', 'source', 'type', 'start', 'end', 'score', 'direction', 'phase', 'attr']


if __name__ == '__main__':
    with open('config.json', 'r') as f:
        config = json.load(f)

    MAB_GFF_PATH = config.get('paths', {}).get('genes', 'MabATCC19977_gff.gff')
    try:
        mab_df = pd.read_table(MAB_GFF_PATH, header=None, names=GFF_HEADERS)
    except FileNotFoundError:
        print(f"Error: The file '{MAB_GFF_PATH}' was not found. Please check the path in the config.json file.")
        exit(1)
    print("Successfully loaded MAB GFF file with column names:", GFF_HEADERS)

    PEAKS_GFF_PATH = config.get('paths', {}).get('peaks', 'SigHP1/SigHP1_FDR0.01combo.gff')
    try:
        peaks_df = pd.read_table(PEAKS_GFF_PATH, header=None, names=GFF_HEADERS)
    except FileNotFoundError:
        print(f"Error: The file '{PEAKS_GFF_PATH}' was not found. Please check the path in the config.json file.")
        exit(1)
    print("Successfully loaded peaks GFF file with column names:", GFF_HEADERS)

    print("Running analysis:")
    regulated_peaks = annotate_peaks(peaks_df, mab_df)
    print("Output preview:")
    print(regulated_peaks.head())
    
    out_path = config.get('paths', {}).get('output', 'regulated_peaks_out.csv')
    regulated_peaks.to_csv(out_path, index=False)
    print("Saved analysis to", out_path)
