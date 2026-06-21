#!/usr/bin/env python3

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp


def select_matrix(adata, layer):
    if layer == "X":
        return adata.X
    if layer == "raw":
        if adata.raw is None:
            raise ValueError("Requested --layer=raw but adata.raw is absent.")
        return adata.raw.X
    if layer not in adata.layers:
        raise ValueError(f"Requested layer {layer!r} is absent. Available layers: {list(adata.layers.keys())}")
    return adata.layers[layer]


def gene_names_for_layer(adata, layer):
    var = adata.raw.var if layer == "raw" and adata.raw is not None else adata.var
    for col in ["feature_name", "gene_name", "gene_symbols", "symbol"]:
        if col in var.columns:
            values = var[col].astype(str).values
            if len(set(values)) == len(values):
                return values
    return var.index.astype(str).values


def main():
    parser = argparse.ArgumentParser(description="Export h5ad counts to FibroDynMix MTX/TSV inputs.")
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--layer", default="X", help="X, raw, or a named AnnData layer")
    parser.add_argument("--cell-type-col", default=None)
    parser.add_argument("--cell-type-regex", default=None)
    parser.add_argument("--max-cells", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=1)
    args = parser.parse_args()

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)

    adata = ad.read_h5ad(args.h5ad)
    if args.cell_type_col and args.cell_type_regex:
        mask = adata.obs[args.cell_type_col].astype(str).str.contains(args.cell_type_regex, case=False, regex=True, na=False)
        adata = adata[mask].copy()

    if adata.n_obs == 0:
        raise ValueError("No cells remain after h5ad filtering.")

    if adata.n_obs > args.max_cells:
        rng = np.random.default_rng(args.seed)
        idx = rng.choice(adata.n_obs, args.max_cells, replace=False)
        adata = adata[idx].copy()

    x = select_matrix(adata, args.layer)
    if not sp.issparse(x):
        x = sp.csr_matrix(x)

    # AnnData stores cells x genes; FibroDynMix expects genes x cells.
    gene_by_cell = x.T.tocoo()
    scipy.io.mmwrite(out / "counts.mtx", gene_by_cell)

    pd.Series(gene_names_for_layer(adata, args.layer)).to_csv(out / "genes.tsv", sep="\t", header=False, index=False)
    pd.Series(adata.obs_names.astype(str)).to_csv(out / "cells.tsv", sep="\t", header=False, index=False)

    metadata = adata.obs.copy()
    metadata.insert(0, "cell_id", adata.obs_names.astype(str))
    metadata.to_csv(out / "metadata.tsv", sep="\t", index=False)

    manifest = pd.DataFrame({
        "field": ["h5ad", "layer", "n_cells", "n_genes"],
        "value": [args.h5ad, args.layer, adata.n_obs, adata.n_vars],
    })
    manifest.to_csv(out / "h5ad_export_manifest.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
