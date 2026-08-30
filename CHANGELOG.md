# Changelog

All notable changes to ReCoN are documented here.

This changelog follows the spirit of
[Keep a Changelog](https://keepachangelog.com/en/1.1.0/) and uses semantic
version headings. New versions should be added above older versions.

## Maintainer Notes

When adding a release entry, keep the structure stable enough that future
maintainers and agents can update it safely:

- Use `## vX.Y.Z - YYYY-MM-DD` for released versions.
- Group changes with conventional headings when they apply:
  - `Added`
  - `Changed`
  - `Fixed`
  - `Deprecated`
  - `Removed`
  - `Tests`
  - `Documentation`
- Write user-facing changes first, then internal/testing details.
- Prefer concise bullets that mention the affected public API, command, or
  workflow.
- Keep GitHub release drafts shorter than the changelog. The release page can
  summarize highlights and link back here for detail.

## v0.0.4 - 2026-08-30

ReCoN v0.0.4 adds guided spatial and cascade-plot workflows and improves the
customization, layout, and diagnostics of biological cascade figures.

### Added

- Added a Visium spatial-niche tutorial covering cell2location deconvolution,
  spot–cell-type expression reconstruction, niche-specific communication and
  ReCoN models, and comparison of focal-cell responses across niches.
- Added a biological cascade-plot tutorial covering intracellular, ligand-to-
  gene, complete intercellular, and two-run contrast views.
- Added `start_layer` and `stop_layer` controls for selecting a cascade window.
- Added `return_edges` so users can inspect the exact paths retained in a plot.
- Added configurable molecular node colors through `node_type_colors`.
- Added warnings when a requested cascade transition has no retained edges,
  including guidance on the relevant `top_*_n` parameter.
- Added a GRNBoost2 tutorial for constructing a GRN from scRNA-seq data alone.

### Changed

- Improved receiver-only cascade geometry and expanded the intracellular
  drawing area.
- Improved gene placement within the receiving-cell nucleus, using compact
  circular layouts for small programs and centered multi-row layouts for
  larger programs.
- Updated the default molecular palette to a ggplot2-inspired color scheme.
- Reduced duplicate labels and improved connected-path selection across legacy
  GRN orientations and layer names.
- Extended tutorial-data upload and download support for the spatial example.

### Fixed

- Fixed TF-to-gene edges that could be selected but omitted from rendering.
- Fixed duplicate receiving-gene nodes and labels in cascade plots.
- Fixed contrast construction when result tables contain duplicate indices.
- Fixed empty sender-layer edge sets causing Sankey construction to fail.

### Tests

- Added regression coverage for layer windows, empty-layer warnings, contrast
  plots, custom node colors, receiver layouts, legacy GRN formats, and empty
  sender networks.

### Documentation

- Added spatial-workflow illustrations, dependency and runtime guidance, and
  links to the corresponding reproducibility showcases.
- Reorganized and clarified the GRN-generation tutorial navigation.

## v0.0.3 - 2026-04-29

ReCoN v0.0.3 expands molecular cascade visualization, improves downstream
cascade logic, adds tutorial data helpers, and substantially strengthens the
test suite and coverage reporting.

### Highlights

- Added biological cascade plotting:
  - `cascade_plot`
  - `contrast_cascade_plot`
- Made Sankey/cascade path construction direction-aware.
- Added one-step exploration helpers for common ReCoN workflows.
- Added packaged receptor-gene prior resources.
- Added tutorial-data download helpers.
- Added Codecov-based coverage reporting and broad unit/golden regression tests.

### Cascade and Sankey Visualization

#### New cascade plots

ReCoN now includes a Matplotlib-based cascade plotting API under `recon.plot`:

- `cascade_plot`: single-run biological cascade diagram with score-based node
  and edge coloring.
- `contrast_cascade_plot`: two-run comparison diagram with diverging
  delta-based coloring.

These functions share a new cascade core for:

- network table formatting
- node set collection
- biological layout geometry
- sender/receiver cell rendering
- legends and contrast color utilities

#### Direction-aware Sankey construction

The Sankey network construction is now direction-aware.

Previously, `flow="downstream"` mostly reversed displayed arrows of an upstream
cascade. The underlying selected paths were still built as upstream paths from
receptors/ligands toward seed genes.

Now, downstream plots build paths using downstream biology from the seed side of
the network:

- `flow="upstream"` keeps the existing behavior:
  - intracellular: Receptor -> TF -> Gene
  - ligand: Ligand -> Receptor -> TF -> Gene
  - intercellular: Upstream Receptor -> Upstream TF -> Ligand -> Receptor -> TF -> Gene

- `flow="downstream"` selects downstream targets from seed genes, then displays:
  - intracellular: Gene -> TF -> Receptor
  - ligand: Gene -> TF -> Receptor -> Ligand
  - intercellular: Gene -> TF -> Receptor -> Ligand -> Upstream TF -> Upstream Receptor

Seed handling is also more flexible. Sankey construction now accepts plain gene
names such as `Tnf` and suffixed names such as `Tnf::Macrophage`, matching the
tutorial notebook usage.

A downstream extraction path was added so seed genes can be matched to TF-layer
syntax such as:

```text
SEED::CellType
SEED_TF::CellType
```

This lets downstream cascades follow edges from seed TF-like nodes to downstream
gene nodes while preserving existing upstream matching logic.

### Explore Workflow Improvements

Added convenience methods for one-step exploration:

- `Celltype.explore`
- `Multicell.explore`

These methods build the HuMMuS/Multixrank object, run random walk with restart,
store the multilayer object and results on the ReCoN object, and return the
result table.

Additional helper improvements:

- Added validation for `set_lambda(direction=...)`:
  - allowed values: `"upstream"`, `"downstream"`
- Added validation for `set_lambda(strategy=...)`:
  - allowed values: `"intracell"`, `"intercell"`
- Fixed `set_lambda(..., multicell=..., celltypes=...)` so it warns and uses
  `multicell` instead of raising a `TypeError`.
- Fixed receptor-prior loading in `multicell_targets` to use
  `load_receptor_genes` and `receptor_gene_resources`.
- Fixed `ccc_to_celltype_proba` indexing after cell types are converted to a
  dictionary.
- Fixed `combine_effects(..., cell_comm_matrix=...)` so weighted indirect
  effects are returned correctly.
- Removed a pandas `SettingWithCopyWarning` in `format_multicell_results`.

### Data and Tutorial Helpers

New packaged receptor-gene priors are included:

- `human_receptor_gene_from_NichenetPKN`
- `mouse_receptor_gene_from_NichenetPKN`

New tutorial data helpers were added:

- `fetch_tutorial_data`
- `fetch_all_tutorial_data`
- `download_tutorial`

The missing-`pooch` error message now points users toward the tutorial extra:

```bash
pip install "recon[tutorials]"
```

### Plotting Utilities

Added result plotting utilities in `plot_results.py`.

Updated the public plotting namespace in `recon.plot` to include:

- multicell illustration tools
- Sankey path plotting tools
- result plotting helpers
- cascade plotting helpers

### Documentation

Documentation and tutorials were expanded, including:

- API pages for explore, plot, inference, and data loading
- molecular cascade tutorial updates
- receptor-gene database documentation
- ReadTheDocs configuration
- docs tooling for Plotly output patching

### Testing and CI

This release adds extensive test coverage across ReCoN:

- cascade plot internals and public plotting functions
- Sankey upstream/downstream construction
- `Celltype` and `Multicell` construction and exploration
- `set_lambda`
- `multicell_targets`
- tutorial data loading helpers
- receptor-gene packaged resources
- plotting utilities
- GRN inference helper functions

Golden regression tests were added for tiny deterministic explore workflows:

- `Celltype.explore`
- `Multicell.explore`
- `multicell_targets`

The golden outputs are stored in:

```text
tests/golden/
```

They can be regenerated manually with:

```bash
python scripts/generate_golden_explore_outputs.py
```

Tutorial download tests are now offline-safe by mocking downloads instead of
requiring live Zenodo access.

CI now includes coverage upload through:

```text
.github/workflows/codecov.yml
```

### Packaging

- Version bumped to `0.0.3`.
- Receptor-gene parquet resources are included as package data.
- Workflow and wheel build configuration were updated.
