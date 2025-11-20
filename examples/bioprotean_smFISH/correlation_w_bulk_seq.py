#!/usr/bin/env python3
# correlation.py
#
# Positional usage:
#   python correlation.py file1a file1b file2 gene_col1a gene_col1b gene_col2 fpkm_col
#
# Example:
#   python correlation.py spots_a.csv spots_b.csv bulk.txt target_molecule_name target_gene gene_id FPKM \
#     --file2-sep "\t" --min-fpkm 0.1 --only-cell-id-positive
#
# Notes:
# - file1a/file1b: .txt / .csv / .parquet (only the specified gene column is used; counts = occurrences)
# - file2: .txt (contains gene + FPKM)
# - Optional: filter file1a/file1b rows to keep only those with cell_id > 0

from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import typer

app = typer.Typer(add_completion=False, pretty_exceptions_enable=False)

def _load_file1(path, sep_override):
    p = Path(path)
    if not p.exists():
        typer.secho("ERROR: file1 not found: {}".format(p), fg=typer.colors.RED)
        raise typer.Exit(code=2)
    ext = p.suffix.lower()
    try:
        if ext == ".parquet":
            df = pd.read_parquet(p)
        elif ext == ".csv":
            df = pd.read_csv(p, sep=sep_override if sep_override else ",")
        elif ext == ".txt":
            if sep_override:
                df = pd.read_csv(p, sep=sep_override)
            else:
                df = pd.read_csv(p, sep=None, engine="python")  # sniff
        else:
            typer.secho("ERROR: file1 must be .txt, .csv, or .parquet (got {}).".format(ext), fg=typer.colors.RED)
            raise typer.Exit(code=2)
    except Exception as e:
        typer.secho("ERROR: failed to read file1 ({}): {}".format(p, e), fg=typer.colors.RED)
        raise typer.Exit(code=2)
    return df

def _load_file2_txt(path, sep_override):
    p = Path(path)
    if not p.exists():
        typer.secho("ERROR: file2 not found: {}".format(p), fg=typer.colors.RED)
        raise typer.Exit(code=2)
    if p.suffix.lower() != ".txt":
        typer.secho("ERROR: file2 must be .txt (got {}).".format(p.suffix), fg=typer.colors.RED)
        raise typer.Exit(code=2)
    try:
        if sep_override:
            return pd.read_csv(p, sep=sep_override)
        return pd.read_csv(p, sep=None, engine="python")  # sniff
    except Exception as e:
        typer.secho("ERROR: failed to read file2 ({}): {}".format(p, e), fg=typer.colors.RED)
        raise typer.Exit(code=2)

def _counts_vs_fpkm(df1, df2, gene_col1, gene_col2, fpkm_col, only_pos_cell, cell_col, label_for_errors):
    if gene_col1 not in df1.columns:
        typer.secho("ERROR: '{}' not in {}. Columns: {}".format(gene_col1, label_for_errors, list(df1.columns)), fg=typer.colors.RED)
        raise typer.Exit(code=2)
    for col in (gene_col2, fpkm_col):
        if col not in df2.columns:
            typer.secho("ERROR: '{}' not in file2. Columns: {}".format(col, list(df2.columns)), fg=typer.colors.RED)
            raise typer.Exit(code=2)

    df1_work = df1
    if only_pos_cell:
        if cell_col not in df1.columns:
            typer.secho("ERROR: Requested cell_id filtering but '{}' not found in {}. Columns: {}".format(cell_col, label_for_errors, list(df1.columns)), fg=typer.colors.RED)
            raise typer.Exit(code=2)
        df1_work = df1.copy()
        # robust numeric filter: cell_id > 0
        df1_work[cell_col] = pd.to_numeric(df1_work[cell_col], errors="coerce")
        df1_work = df1_work[df1_work[cell_col] > 0]
        if df1_work.empty:
            typer.secho("ERROR: After '{}' > 0 filtering, {} has no rows.".format(cell_col, label_for_errors), fg=typer.colors.RED)
            raise typer.Exit(code=2)

    counts_df = (
        df1_work[gene_col1].dropna().astype(str).value_counts()
        .rename_axis("gene_id").reset_index(name="count")
    )
    counts_df["count"] = counts_df["count"].astype("int64")

    df2 = df2.copy()
    df2[gene_col2] = df2[gene_col2].astype(str)
    fpkm_by_gene = (
        df2.groupby(gene_col2, dropna=False, as_index=False)[fpkm_col]
        .mean()
        .rename(columns={gene_col2: "gene_id", fpkm_col: "fpkm"})
    )

    merged = counts_df.merge(fpkm_by_gene, on="gene_id", how="inner")
    merged = merged.dropna(subset=["count", "fpkm"])
    merged = merged[(merged["count"] > 0) & (merged["fpkm"] > 0)]
    if merged.empty:
        typer.secho("ERROR: No overlapping genes with positive counts and positive FPKM for {}.".format(label_for_errors), fg=typer.colors.RED)
        raise typer.Exit(code=2)
    return merged

def _pearson_loglog_x_fpkm_y_counts(merged):
    # Pearson correlation in log10 space with x=FPKM, y=counts
    x = merged["fpkm"].to_numpy(dtype=float)
    y = merged["count"].to_numpy(dtype=float)
    logx = np.log10(x)
    logy = np.log10(y)
    return float(np.corrcoef(logx, logy)[0, 1])

@app.command(help="Two file1 inputs (possibly different gene columns) vs one file2: FPKM on X, counts on Y; filter by FPKM and optional cell_id>0; no trend lines; 2× styling; only r.")
def main(
    file1a = typer.Argument(..., help="First file1 (.txt/.csv/.parquet)."),
    file1b = typer.Argument(..., help="Second file1 (.txt/.csv/.parquet)."),
    file2   = typer.Argument(..., help="file2 (.txt) with gene IDs and FPKM."),
    gene_col1a = typer.Argument(..., help="Column in file1a containing gene IDs."),
    gene_col1b = typer.Argument(..., help="Column in file1b containing gene IDs."),
    gene_col2 = typer.Argument(..., help="Column in file2 containing gene IDs."),
    fpkm_col  = typer.Argument(..., help="Column in file2 containing FPKM."),
    plot_out = typer.Option("fpkm_vs_counts_dual.png", help="Output plot filename."),
    file1a_sep = typer.Option("", help="Optional delimiter for file1a (.txt/.csv). E.g., '\\t' or ','."),
    file1b_sep = typer.Option("", help="Optional delimiter for file1b (.txt/.csv). E.g., '\\t' or ','."),
    file2_sep  = typer.Option("", help="Optional delimiter for file2 (.txt). E.g., '\\t' or ','."),
    min_fpkm = typer.Option(0.0, help="Minimum FPKM threshold to include (e.g., 0.1 for 10^-1). Genes with FPKM <= threshold are excluded."),
    only_cell_id_positive: bool = typer.Option(False, "--only-cell-id-positive / --no-only-cell-id-positive", help="If set, keep only rows with cell_id > 0 in file1a and file1b before counting."),
    cellid_col1a = typer.Option("cell_id", help="Cell ID column name in file1a (used only if --only-cell-id-positive)."),
    cellid_col1b = typer.Option("cell_id", help="Cell ID column name in file1b (used only if --only-cell-id-positive)."),
):
    # Load
    df1a = _load_file1(file1a, file1a_sep)
    df1b = _load_file1(file1b, file1b_sep)
    df2  = _load_file2_txt(file2, file2_sep)
    min_fpkm = float(min_fpkm)

    # Compute merged tables (apply cell_id filter if requested)
    merged_a = _counts_vs_fpkm(df1a, df2, gene_col1a, gene_col2, fpkm_col,
                               only_cell_id_positive, cellid_col1a, "file1a")
    merged_b = _counts_vs_fpkm(df1b, df2, gene_col1b, gene_col2, fpkm_col,
                               only_cell_id_positive, cellid_col1b, "file1b")

    # Apply FPKM threshold filter
    if min_fpkm > 0.0:
        merged_a = merged_a[merged_a["fpkm"] > min_fpkm]
        merged_b = merged_b[merged_b["fpkm"] > min_fpkm]
    if merged_a.empty:
        typer.secho("ERROR: After filtering (min_fpkm={}), no genes remain for file1a.".format(min_fpkm), fg=typer.colors.RED)
        raise typer.Exit(code=2)
    if merged_b.empty:
        typer.secho("ERROR: After filtering (min_fpkm={}), no genes remain for file1b.".format(min_fpkm), fg=typer.colors.RED)
        raise typer.Exit(code=2)

    # Correlations (log10 space; x=FPKM, y=counts)
    rA = _pearson_loglog_x_fpkm_y_counts(merged_a)
    rB = _pearson_loglog_x_fpkm_y_counts(merged_b)

    # Fixed labels
    labelA = "merlin"
    labelB = "merfish3d"

    # ---- Styling: ~2× bigger everything ----
    scale = 2.0
    base_font = 12.0
    plt.rcParams.update({
        "font.size": base_font * scale,
        "axes.labelsize": base_font * scale,
        "xtick.labelsize": base_font * scale,
        "ytick.labelsize": base_font * scale,
        "legend.fontsize": base_font * scale,
    })
    marker_area = 12 * (scale ** 2)
    spine_lw = 1.0 * scale
    tick_len = 3.5 * scale
    tick_w = 1.0 * scale
    ann_font = base_font * scale

    # Plot
    fig, ax = plt.subplots(figsize=(7.5 * scale, 5.5 * scale), dpi=150)

    colorA = "C0"
    colorB = "C1"

    # Scatter: x = FPKM, y = counts
    ax.scatter(merged_a["fpkm"], merged_a["count"], s=marker_area, alpha=0.7, marker="o",
               label=labelA, edgecolors="none", c=colorA)
    ax.scatter(merged_b["fpkm"], merged_b["count"], s=marker_area, alpha=0.7, marker="s",
               label=labelB, edgecolors="none", c=colorB)

    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Bulk RNA-Seq (FPKM)")
    ax.set_ylabel("MERFISH (total counts in ROI)")

    # Minimal look: remove top/right box, thicken remaining axes
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_linewidth(spine_lw)
    ax.spines["bottom"].set_linewidth(spine_lw)
    ax.tick_params(axis="both", which="both", length=tick_len, width=tick_w)

    # Only correlation coefficients on the plot
    ax.text(0.03, 0.97, "{}: r = {:.3f}".format(labelA, rA),
            transform=ax.transAxes, va="top", ha="left", color=colorA, fontsize=ann_font)
    ax.text(0.03, 0.90, "{}: r = {:.3f}".format(labelB, rB),
            transform=ax.transAxes, va="top", ha="left", color=colorB, fontsize=ann_font)

    # Legend outside
    ax.legend(frameon=False, loc="center left", bbox_to_anchor=(1.02, 0.5))

    fig.tight_layout()
    fig.savefig(plot_out, bbox_inches="tight")
    plt.close(fig)
    typer.secho("Wrote plot: {}".format(plot_out), fg=typer.colors.GREEN)

if __name__ == "__main__":
    app()
