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

def find_regulated_genes(
    peak_average,
    gene_table,
    forward_proximity=700,
    antisense_proximity=300,
):
    """
    - If the peak average (Pavg) lies completely within the boundaries of a gene X, it does regulate that gene X.
    - Additionally (independently of the first point), if Pavg is within 700 bp (or the int set by the const `forward_proximity`) of the **beginning** of a gene Y, then it regulates that gene Y.
    - Also to be applied independently: if Pavg is within 300 bp (or the const `antisense_proximity`) of the **end** of a gene Z, then it regulates that gene Z.
    Note that "beginning" is different for a gene in the + vs - direction. For + direction, it means leftmost (smaller) coordinate. For - direction, it is the opposite. Same deal for "end" of a gene.
    """
    genes = gene_table["genes"]
    starts = genes["gene_start"].to_numpy()
    ends = genes["gene_end"].to_numpy()
    directions = genes["direction"].to_numpy()

    inside_gene = (starts <= peak_average) & (peak_average <= ends)
    near_beginning = (
        ((directions == "+") & (starts > peak_average) & (starts - peak_average <= forward_proximity))
        | ((directions == "-") & (ends < peak_average) & (peak_average - ends <= forward_proximity))
    )
    near_end = (
        ((directions == "+") & (ends < peak_average) & (peak_average - ends <= antisense_proximity))
        | ((directions == "-") & (starts > peak_average) & (starts - peak_average <= antisense_proximity))
    )

    candidates = genes.loc[
        inside_gene | near_beginning | near_end, "gene_name"
    ].astype(str).tolist()

    return "-" if not candidates else "/".join(candidates)


def annotate_peaks(
    peaks_df,
    mab_df,
    forward_proximity=700,
    antisense_proximity=300,
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
    
    tqdm.pandas(
        desc=(
            f"Finding regulated genes with forward_proximity={forward_proximity}, "
            f"antisense_proximity={antisense_proximity}"
        )
    )
    result["Gene(s) Regulated"] = result.progress_apply(
        lambda row: find_regulated_genes(
            row["Paverage"],
            gene_table,
            forward_proximity=forward_proximity,
            antisense_proximity=antisense_proximity,
        ),
        axis=1,
    )

    result['Coordinates'] = result.apply(lambda r: f"{r['P1']}_{r['P2']}", axis=1)
    result.drop(columns=['P1', 'P2', 'Paverage'], inplace=True)
    # move Coordinates to be first column
    result = result[['Coordinates'] + [col for col in result.columns if col != 'Coordinates']]
    
    return result
