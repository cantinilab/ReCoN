"""Build the compact, deterministic data bundle used by the spatial tutorial."""

from __future__ import annotations

import csv
import hashlib
import heapq
import json
import shutil
from pathlib import Path

import anndata as ad
import numpy as np


SEED = 20260825
N_REFERENCE_CELLS = 4000
N_GENES = 4500
MAX_GRN_EDGES = 250_000
CELL_TYPES = [
    "Cardiomyocyte",
    "Endothelial",
    "Fibroblast",
    "Lymphoid",
    "Myeloid",
    "Pericyte",
    "vSMCs",
]

HERE = Path(__file__).resolve().parent
SHOWCASE = HERE.parents[5] / "recon_reproducibility" / "new_showcases"
RESULTS = SHOWCASE / "results" / "cardiomyocyte_subtypes"
# This is the count-scale copy produced immediately before reference-model fitting.
# The public heart-map object's X is log-normalized (its raw slot is restored in
# the showcase notebook), so subsetting that X directly would be incorrect.
REF_FILE = RESULTS / "reference_signatures" / "sc.h5ad"
SPATIAL_FILE = RESULTS / "cell2location_map" / "sp_Visium_18.h5ad"
GRN_FILE = SHOWCASE / "data" / "heart_model" / "data" / "hummus.hummus.hummus.hummus.mdl.csv"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def allocate_reference(obs, rng: np.random.Generator) -> np.ndarray:
    """Balance cell types while retaining every type-by-sample stratum."""
    eligible = obs["cell_type_original"].isin(CELL_TYPES)
    selected: list[int] = []
    # An equal type budget prevents abundant cardiomyocytes from crowding out
    # lymphoid/vSMC cells. Within a type, sampling remains proportional to sample.
    base, remainder = divmod(N_REFERENCE_CELLS, len(CELL_TYPES))
    for type_i, cell_type in enumerate(CELL_TYPES):
        type_budget = base + (type_i < remainder)
        positions = np.flatnonzero((obs["cell_type_original"] == cell_type).to_numpy() & eligible.to_numpy())
        samples = obs.iloc[positions]["sample"].astype(str).to_numpy()
        strata = {sample: positions[samples == sample] for sample in sorted(set(samples))}
        # Guarantee one cell per observed sample, then distribute proportionally.
        allocation = {sample: 1 for sample in strata}
        left = type_budget - len(allocation)
        weights = np.array([len(strata[s]) - 1 for s in strata], dtype=float)
        if weights.sum() > 0:
            raw = left * weights / weights.sum()
            floors = np.floor(raw).astype(int)
            for sample, count in zip(strata, floors):
                allocation[sample] += int(count)
            for i in np.argsort(-(raw - floors))[: left - floors.sum()]:
                allocation[list(strata)[i]] += 1
        for sample, candidates in strata.items():
            selected.extend(rng.choice(candidates, allocation[sample], replace=False).tolist())
    return np.array(sorted(selected), dtype=int)


def communication_genes() -> set[str]:
    genes: set[str] = set()
    for niche in (1, 5):
        path = RESULTS / f"cell_communication_niche_{niche}_Visium_18.csv"
        with path.open(newline="") as handle:
            for row in csv.DictReader(handle):
                genes.update((row["source"], row["target"]))
    return genes


def main() -> None:
    rng = np.random.default_rng(SEED)
    ref_backed = ad.read_h5ad(REF_FILE, backed="r")
    spatial_backed = ad.read_h5ad(SPATIAL_FILE, backed="r")

    ref_symbols = ref_backed.var["feature_name"].astype(str)
    spatial_symbols = spatial_backed.var["SYMBOL"].astype(str)
    ref_by_symbol = dict(zip(ref_symbols, ref_backed.var_names))
    spatial_by_symbol = dict(zip(spatial_symbols, spatial_backed.var_names))
    shared = set(ref_by_symbol) & set(spatial_by_symbol)
    ccc = communication_genes() & shared

    # Spatial detection is a transparent, stable proxy for informative genes.
    ranked = spatial_backed.var.loc[
        spatial_backed.var["SYMBOL"].isin(shared)
    ].sort_values(["n_cells_by_counts", "total_counts"], ascending=False)["SYMBOL"].astype(str)
    symbols = list(dict.fromkeys(list(sorted(ccc)) + ranked.tolist()))[:N_GENES]
    # The sorted final order makes the same feature axis across both AnnData files.
    symbols = sorted(symbols)
    ref_ids = [ref_by_symbol[s] for s in symbols]
    spatial_ids = [spatial_by_symbol[s] for s in symbols]

    selected = allocate_reference(ref_backed.obs, rng)
    reference = ref_backed[selected, ref_ids].to_memory()
    reference.var_names = symbols
    reference.var_names.name = "gene_symbol"
    reference.var = reference.var[["feature_name", "feature_biotype"]].copy()
    reference.obs = reference.obs[["sample", "donor_id", "patient_group", "major_labl", "cell_type_original"]].copy()
    reference.obsm.clear()
    reference.uns.clear()
    reference.write_h5ad(HERE / "reference_heart_4k.h5ad", compression="gzip")

    spatial = spatial_backed[:, spatial_ids].to_memory()
    spatial.var_names = symbols
    spatial.var_names.name = "gene_symbol"
    spatial.var = spatial.var[["SYMBOL", "n_cells_by_counts", "total_counts"]].copy()
    abundance = CELL_TYPES
    spatial.obs = spatial.obs[["array_row", "array_col", "sample", "celltype_niche", "molecular_niche", *abundance]].copy()
    coords = spatial.obsm["spatial"].copy()
    spatial.obsm.clear()
    spatial.obsm["spatial"] = coords
    spatial.layers.clear()
    retained_uns = {
        key: spatial.uns[key]
        for key in ("spatial", "celltype_niche_colors", "neighbors")
        if key in spatial.uns
    }
    spatial.uns.clear()
    spatial.uns.update(retained_uns)
    spatial.write_h5ad(HERE / "visium18_all_spots_4k_genes.h5ad", compression="gzip")
    ref_backed.file.close()
    spatial_backed.file.close()

    for name in [
        "proportions.csv",
        "compositions.csv",
        "avg_expression_niche_1_Visium_18.csv",
        "avg_expression_niche_5_Visium_18.csv",
        "cell_communication_niche_1_Visium_18.csv",
        "cell_communication_niche_5_Visium_18.csv",
    ]:
        shutil.copy2(RESULTS / name, HERE / name)

    # Keep the strongest induced GRN edges. This preserves complete edges rather
    # than copying the 510 MB network, and its node set matches the tutorial genes.
    heap: list[tuple[float, str, str]] = []
    gene_set = set(symbols)
    with GRN_FILE.open(newline="") as handle:
        for row in csv.DictReader(handle):
            if row["source"] not in gene_set or row["target"] not in gene_set:
                continue
            item = (float(row["score"]), row["source"], row["target"])
            if len(heap) < MAX_GRN_EDGES:
                heapq.heappush(heap, item)
            elif item[0] > heap[0][0]:
                heapq.heapreplace(heap, item)
    with (HERE / "heart_grn_top250k.csv").open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["source", "target", "score"])
        writer.writerows(
            (source, target, score)
            for score, source, target in sorted(heap, reverse=True)
        )

    artifacts = {}
    for path in sorted(HERE.iterdir()):
        if path.is_file() and path.name not in {"build_mini_data.py", "manifest.json", "README.md"}:
            artifacts[path.name] = {"bytes": path.stat().st_size, "sha256": sha256(path)}
    manifest = {
        "seed": SEED,
        "reference_cells": reference.n_obs,
        "spatial_spots": spatial.n_obs,
        "genes": len(symbols),
        "cell_types": reference.obs["cell_type_original"].value_counts().sort_index().to_dict(),
        "reference_samples": reference.obs["sample"].nunique(),
        "spatial_niches": spatial.obs["celltype_niche"].value_counts().sort_index().to_dict(),
        "ccc_genes_present": len(ccc),
        "grn_edges": len(heap),
        "artifacts": artifacts,
    }
    (HERE / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
