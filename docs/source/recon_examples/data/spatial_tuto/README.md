# Mini data for the spatial ReCoN tutorial

This Zenodo-hosted bundle is a compact derivative of the spatial-heart
reproducibility showcase. It is downloaded by the tutorial while retaining
the complete Visium 18 tissue section and every modeled cell type. This
directory keeps its deterministic build recipe and provenance manifest; the
generated datasets themselves are not stored in Git.

## Contents

- `reference_heart_4k.h5ad`: raw counts for 4,000 single nuclei and 4,500
  genes. The seven cell types used by the cell2location model are balanced
  (571–572 cells each), and every cell type contains cells from all 29 samples.
- `visium18_all_spots_4k_genes.h5ad`: raw Visium counts for all 3,572 spots,
  with spatial coordinates, the spatial neighbor graph, image metadata,
  cell2location abundance estimates, and cell-type/molecular niche labels.
- `avg_expression_niche_{1,5}_Visium_18.csv`: precomputed niche averages for
  the tutorial's quick path and result checking.
- `cell_communication_niche_{1,5}_Visium_18.csv`: precomputed niche-specific
  communication networks.
- `proportions.csv` and `compositions.csv`: niche-level cell-type summaries.
- `heart_grn_top250k.csv`: strongest 250,000 heart-GRN edges whose two genes
  are present on the reduced feature axis.
- `manifest.json`: dimensions, coverage, byte sizes, and SHA-256 checksums.
- `build_mini_data.py`: deterministic recipe used to create all files.

The precomputed spot–cell-type expression profiles are intentionally not
copied: the full object is 1.17 GB and is precisely the intermediate that the
tutorial teaches users to generate. The much smaller niche-average tables
provide a fast path for users who cannot run cell2location.

## Selection method

The recipe uses NumPy seed `20260825`. It assigns an equal reference-cell
budget to each of the seven modeled cell types. Within each cell type, it first
retains at least one cell from every sample in which that type occurs and then
allocates remaining cells in proportion to the original type-by-sample sizes.
This avoids losing rare populations or donor coverage while keeping the model
fit suitable for a tutorial.

The complete slide is retained because subsampling spots would alter spatial
adjacency and could erase the small niches (niche 2 has only 24 spots). Both
AnnData objects use the same 4,500 gene symbols. The axis includes **all 758
ligand/receptor genes shared by the niche-1 and niche-5 communication tables
and both source objects**, followed by genes ranked by spatial detection
(`n_cells_by_counts`, then `total_counts`). Files are gzip-compressed.

## Provenance

Paths below are relative to `recon_reproducibility/new_showcases/`:

- count-scale scRNA-seq reference:
  `results/cardiomyocyte_subtypes/reference_signatures/sc.h5ad`
  (185,137 cells × 13,994 genes). This is used instead of the main heart-map
  object's log-normalized `X`.
- processed Visium 18:
  `results/cardiomyocyte_subtypes/cell2location_map/sp_Visium_18.h5ad`
  (3,572 spots × 11,613 genes).
- full spot–cell-type profiles, audited but not copied:
  `results/cardiomyocyte_subtypes/cell2location_map/spot_celltype_pseudobulks_Visium_18.h5ad`
  (25,004 spot–cell-type profiles × 11,613 genes).
- CCC, expression, proportions, and compositions:
  `results/cardiomyocyte_subtypes/`.
- full heart GRN:
  `data/heart_model/data/hummus.hummus.hummus.hummus.mdl.csv` (510 MB).

Run the recipe from the ReCoN repository root with an environment containing
compatible `anndata`, `pandas`, `numpy`, and `scipy`. It expects the
`recon_reproducibility` repository beside ReCoN.

## Validation

Both generated AnnData files were reopened after writing. Observation and
feature names are unique, matrices contain no missing or negative values, the
reference matrix is integer count-scale, and the Visium object contains all
3,572 unique spots and all nine original cell-type niches. See the manifest
for exact checksums and coverage counts.
