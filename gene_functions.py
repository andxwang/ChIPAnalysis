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


def find_operon(peak_start, peak_end, genes, operon_gap=100, proximity=600):
    """
    Return the operon associated with a peak.

    An operon is a chain of >=2 adjacent genes where every neighboring pair:
      - is separated by <= operon_gap coordinates
      - has the same direction

    A peak can identify an operon either by:
      - directly overlapping a gene, or
      - being within `proximity` of a gene boundary in the gene's
        transcriptional direction.

    Returns
    -------
    str
        "<first gene> - <last gene>" if an operon is found,
        otherwise "-".
    """

    def expand_operon(idx):
        """Expand from a seed gene index in both directions."""
        cluster = [idx]
        direction = genes.iloc[idx]["direction"]

        # Expand left
        left = idx
        while left > 0:
            prev = genes.iloc[left - 1]
            curr = genes.iloc[left]

            gap = curr["gene_start"] - prev["gene_end"]

            if gap <= operon_gap and prev["direction"] == direction:
                cluster.insert(0, left - 1)
                left -= 1
            else:
                break

        # Expand right
        right = idx
        while right < len(genes) - 1:
            curr = genes.iloc[right]
            nxt = genes.iloc[right + 1]

            gap = nxt["gene_start"] - curr["gene_end"]

            if gap <= operon_gap and nxt["direction"] == direction:
                cluster.append(right + 1)
                right += 1
            else:
                break

        return cluster if len(cluster) >= 2 else []

    # First: direct overlap
    overlapping = genes[
        (genes["gene_start"] <= peak_end) &
        (genes["gene_end"] >= peak_start)
    ]

    seed_indices = list(overlapping.index)

    # need to do this even if no overlap
    # Candidate genes near the peak boundaries
    candidates = genes[
        (abs(genes["gene_start"] - peak_start) <= proximity) |
        (abs(genes["gene_end"] - peak_end) <= proximity)
    ]

    # Apply strand-specific regulatory boundary logic
    for idx, gene in candidates.iterrows():
        if gene["direction"] == "+":
            if abs(peak_start - gene["gene_start"]) <= proximity:
                seed_indices.append(idx)

        elif gene["direction"] == "-":
            if abs(peak_end - gene["gene_end"]) <= proximity:
                seed_indices.append(idx)

    if not seed_indices:
        return "-"

    operon_indices = set()

    for idx in seed_indices:
        cluster = expand_operon(idx)
        operon_indices.update(cluster)

    if not operon_indices:
        return "-"

    operon_indices = sorted(operon_indices)
    first_gene = genes.iloc[operon_indices[0]]["gene_name"]
    last_gene = genes.iloc[operon_indices[-1]]["gene_name"]

    return f"{first_gene} - {last_gene}"


def annotate_peaks(peaks_df, mab_df, proximity=600, operon_gap=30):
    gene_table = build_gene_table(mab_df)
    genes = gene_table["genes"]

    result = peaks_df[["start", "end", "score"]].copy()
    result.columns = ["P1", "P2", "Score"]
    result["Paverage"] = (result["P1"] + result["P2"]) / 2
    
    tqdm.pandas(desc="Finding peak locations")
    result["Peak Location"] = result.progress_apply(
        lambda row: get_peak_location(row["P1"], row["P2"], genes),
        axis=1,
    )
    
    tqdm.pandas(desc=f"Finding regulated genes with proximity={proximity}")
    result["Gene(s) Regulated"] = result.progress_apply(
        lambda row: find_regulated_genes(row["P1"], row["P2"], gene_table, proximity=proximity),
        axis=1,
    )
    result['comments'] = pd.qcut(result['Score'], q=5, labels=['no real peak', 'small', 'medium', 'large', 'very large'])
    
    tqdm.pandas(desc=f"Finding operons with operon_gap={operon_gap}")
    result["Operon"] = result.progress_apply(
        lambda row: find_operon(
            row["P1"], row["P2"], genes, operon_gap=operon_gap
        ),
        axis=1,
    )
    return result
