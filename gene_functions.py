import re
import pandas as pd
from tqdm.auto import tqdm


def parse_gene_name(attr):
    text = str(attr)
    match = re.search(r"name=([^;]+)", text)
    return match.group(1) if match else None


def build_gene_table(mab_df):
    genes = mab_df[["start", "end", "direction", "attr"]].copy()
    genes = genes.rename(columns={"start": "gene_start", "end": "gene_end"})
    genes["gene_name"] = genes["attr"].apply(parse_gene_name)
    genes = genes.dropna(subset=["gene_name"]).copy()
    genes["gene_start"] = pd.to_numeric(genes["gene_start"], errors="coerce")
    genes["gene_end"] = pd.to_numeric(genes["gene_end"], errors="coerce")
    genes = genes.sort_values(["gene_start", "gene_end"]).reset_index(drop=True)
    return genes


def is_intragenic(peak_start, peak_end, gene_start, gene_end):
    return gene_start <= peak_start and peak_end <= gene_end


def get_peak_location(peak_start, peak_end, genes):
    containing = [
        row
        for _, row in genes.iterrows()
        if is_intragenic(peak_start, peak_end, row["gene_start"], row["gene_end"])
    ]

    if len(containing) == 1:
        return f"{containing[0]['gene_name']} (intragenic)"

    overlapping = [
        row
        for _, row in genes.iterrows()
        if max(peak_start, row["gene_start"]) <= min(peak_end, row["gene_end"])
    ]

    if len(overlapping) == 2:
        return f"{overlapping[0]['gene_name']}/{overlapping[1]['gene_name']} (intergenic)"

    left_candidates = genes[genes["gene_end"] < peak_start]
    right_candidates = genes[genes["gene_start"] > peak_end]

    left_gene = left_candidates.iloc[-1]["gene_name"] if not left_candidates.empty else None
    right_gene = right_candidates.iloc[0]["gene_name"] if not right_candidates.empty else None

    if left_gene is not None and right_gene is not None:
        return f"{left_gene}/{right_gene} (intergenic)"

    if left_gene is not None:
        return f"{left_gene}/ - (intergenic)"

    if right_gene is not None:
        return f"-/ {right_gene} (intergenic)"

    return "-"

def is_in_first_third(peak_start, peak_end, gene_start, gene_end):
    if peak_end <= gene_start:
        return False
    if peak_start >= gene_end:
        return False

    gene_length = gene_end - gene_start
    if gene_length <= 0:
        return False

    first_third_end = gene_start + (gene_length / 3)
    overlap_start = max(peak_start, gene_start)
    overlap_end = min(peak_end, gene_end)

    return overlap_end > overlap_start and overlap_start < first_third_end


def is_in_last_third(peak_start, peak_end, gene_start, gene_end):
    if peak_start >= gene_end:
        return False
    if peak_end <= gene_start:
        return False

    gene_length = gene_end - gene_start
    if gene_length <= 0:
        return False

    last_third_start = gene_end - (gene_length / 3)
    overlap_start = max(peak_start, gene_start)
    overlap_end = min(peak_end, gene_end)

    return overlap_end > overlap_start and overlap_end > last_third_start


def find_regulated_genes(peak_start, peak_end, genes, proximity=600):
    candidates = []

    for _, row in genes.iterrows():  # TODO: use binary search or something more efficient
        if row["direction"] == "+":  # forward
            in_first_third = is_in_first_third(
                peak_start, peak_end, row["gene_start"], row["gene_end"]
            )
            near_left_boundary = abs(peak_start - row["gene_start"]) <= proximity
            # near boundary should only be checked if the peak is in front of the gene
            if in_first_third or (near_left_boundary and peak_end <= row['gene_start']):
                candidates.append(row["gene_name"])

        if row["direction"] == "-":  # backward
            if row['gene_name'] == 'MAB_0495c':
                print(f"Checking gene {row['gene_name']} with peak ({peak_start}, {peak_end})")
            in_last_third = is_in_last_third(
                peak_start, peak_end, row["gene_start"], row["gene_end"]
            )
            near_right_boundary = abs(peak_end - row["gene_end"]) <= proximity
            if in_last_third or (near_right_boundary and peak_start >= row['gene_end']):
                candidates.append(row["gene_name"])

    if not candidates:
        return "-"

    return "/".join(candidates)


def annotate_peaks(peaks_df, mab_df):
    tqdm.pandas()
    genes = build_gene_table(mab_df)

    result = peaks_df[["start", "end", "score"]].copy()
    result.columns = ["P1", "P2", "Score"]
    result["Paverage"] = (result["P1"] + result["P2"]) / 2
    result["Peak Location"] = result.progress_apply(
        lambda row: get_peak_location(row["P1"], row["P2"], genes),
        axis=1,
    )
    result["Gene(s) Regulated"] = result.progress_apply(
        lambda row: find_regulated_genes(row["P1"], row["P2"], genes),
        axis=1,
    )
    return result
