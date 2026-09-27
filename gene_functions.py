import re
import numpy as np
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

    return {
        "genes": genes,
        "starts": genes["gene_start"].to_numpy(),
        "ends": genes["gene_end"].to_numpy(),
    }

def get_genetic_location(
    peak_average,
    genes,
    location_edge_cutoff,
):
    """Classify a peak using the average of its start and end coordinates.

    Return the containing gene name, adding ``start of`` or ``end of`` when
    the midpoint is within the corresponding third and configured cutoff.
    Prefixes follow transcription direction, so they reverse for ``-`` genes.
    Join names with ``/`` for overlapping genes or for the nearest flanking
    genes; use ``-`` when an outside flank is absent.
    """
    containing = genes[
        (genes["gene_start"] <= peak_average)
        & (peak_average <= genes["gene_end"])
    ]

    if len(containing) == 1:
        gene = containing.iloc[0]
        gene_start = gene["gene_start"]
        gene_end = gene["gene_end"]
        gene_length = gene_end - gene_start
        gene_name = gene["gene_name"]
        start_cutoff = min(gene_length / 3, location_edge_cutoff)

        if gene["direction"] == "-":
            if gene_end - peak_average < start_cutoff:
                return f"start of {gene_name}"
            if peak_average - gene_start < start_cutoff:
                return f"end of {gene_name}"
        else:
            if peak_average - gene_start < start_cutoff:
                return f"start of {gene_name}"
            if gene_end - peak_average < start_cutoff:
                return f"end of {gene_name}"
        return gene_name

    if len(containing) > 1:
        return "/".join(containing["gene_name"].astype(str))

    left_candidates = genes[genes["gene_end"] < peak_average]
    right_candidates = genes[genes["gene_start"] > peak_average]

    left_gene = left_candidates.iloc[-1]["gene_name"] if not left_candidates.empty else None
    right_gene = right_candidates.iloc[0]["gene_name"] if not right_candidates.empty else None

    if left_gene is not None and right_gene is not None:
        return f"{left_gene}/{right_gene}"

    if left_gene is not None:
        return f"{left_gene}/-"

    if right_gene is not None:
        return f"-/{right_gene}"

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


def is_regulating_gene(peak_start, peak_end, gene_row, proximity=600):
    direction = gene_row["direction"]
    gene_start = gene_row["gene_start"]
    gene_end = gene_row["gene_end"]

    if direction == "+":
        in_first_third = is_in_first_third(peak_start, peak_end, gene_start, gene_end)
        near_left_boundary = abs(peak_start - gene_start) <= proximity
        return in_first_third or (near_left_boundary and peak_end <= gene_start)

    if direction == "-":
        in_last_third = is_in_last_third(peak_start, peak_end, gene_start, gene_end)
        near_right_boundary = abs(peak_end - gene_end) <= proximity
        return in_last_third or (near_right_boundary and peak_start >= gene_end)

    return False


def find_regulated_genes(peak_start, peak_end, gene_table, proximity=600):
    genes = gene_table["genes"]
    starts = gene_table["starts"]
    ends = gene_table["ends"]

    candidate_indices = set()

    # Genes whose start is close enough to possibly satisfy overlap or + strand proximity.
    left = np.searchsorted(starts, peak_start - proximity, side="left")
    right = np.searchsorted(starts, peak_end + proximity, side="right")
    candidate_indices.update(range(left, right))

    # Genes whose end is close enough to possibly satisfy overlap or - strand proximity.
    left = np.searchsorted(ends, peak_start - proximity, side="left")
    right = np.searchsorted(ends, peak_end + proximity, side="right")
    candidate_indices.update(range(left, right))

    candidates = [
        genes.iloc[i]["gene_name"]
        for i in sorted(candidate_indices)
        if is_regulating_gene(
            peak_start,
            peak_end,
            genes.iloc[i],
            proximity=proximity,
        )
    ]

    return "-" if not candidates else "/".join(candidates)


def annotate_peaks(
    peaks_df,
    mab_df,
    proximity=600,
    *,
    location_edge_cutoff,
):
    gene_table = build_gene_table(mab_df)
    genes = gene_table["genes"]

    result = peaks_df[["start", "end", "score"]].copy()
    result.columns = ["P1", "P2", "Score"]
    result["Paverage"] = (result["P1"] + result["P2"]) / 2
    
    tqdm.pandas(desc="Finding genetic locations")
    result["Genetic Location"] = result.progress_apply(
        lambda row: get_genetic_location(
            row["Paverage"], genes, location_edge_cutoff=location_edge_cutoff
        ),
        axis=1,
    )
    
    tqdm.pandas(desc=f"Finding regulated genes with proximity={proximity}")
    result["Gene(s) Regulated"] = result.progress_apply(
        lambda row: find_regulated_genes(row["P1"], row["P2"], gene_table, proximity=proximity),
        axis=1,
    )
    result['comments'] = pd.qcut(result['Score'], q=5, labels=['no real peak', 'small', 'medium', 'large', 'very large'])
    
    return result.drop(columns=["P1", "P2"])
