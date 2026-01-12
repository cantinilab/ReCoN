# ReCoN: AI Coding Agent Instructions

## Project Overview

**ReCoN** reconstructs multicellular models by integrating **gene regulatory networks (GRNs)** and **cell communication networks (CCNs)**. It uses heterogeneous multilayer network analysis to explore molecular coordination between cell types from single-cell RNA-seq data.

Core architecture: Multilayer network → Random Walk with Restart (RWR) → Treatment effect prediction & molecular cascade reconstruction.

## Architecture

### Core Components (src/recon/)

1. **explore/recon.py** - Central module defining two key classes:
   - `Celltype`: Represents a single cell type's multilayer network (receptor layer → GRN layer)
   - `Multicell`: Extends `Celltype` to model interactions across multiple cell types via cell communication layer
   
2. **infer_grn/layers.py** - GRN inference using:
   - HuMMuS methodology (multilayer TF→DNA→gene)
   - CellOracle for TF-DNA motif links
   - Arboreto/GRNBoost2 for TF-gene relationships
   
3. **plot/** - Visualization utilities:
   - `plot_multicell.py`: 3D multilayer network visualization
   - `sankey_paths.py`: Molecular cascade/pathway Sankey diagrams

### Data Flow Pattern

```
scRNA-seq data → GRN inference → Celltype objects → Multicell integration → 
RWR exploration (via hummuspy/multixrank) → Treatment effect scores → Visualization
```

## Critical Conventions

### Network Representation

- All networks are **pandas DataFrames** with columns: `source`, `target`, `weight` (optional)
- Node naming conventions:
  - **Multicell genes**: `{gene_name}::{celltype_name}` (double colon) to avoid conflicts between cell types
  - **Receptors**: `{receptor_name}_receptor` suffix (underscore)
  - **Cell communication ligands**: `{ligand_name}-{celltype_source}` (hyphen) format
- These suffixes are automatically added during object construction - do NOT add them manually to input DataFrames

### Graph Types Encoding

Graph types encoded as 2-digit strings (e.g., "01", "11"):
- First digit: `1` if directed, `0` if undirected  
- Second digit: `1` if weighted, `0` if unweighted
- Example: `"01"` = undirected, weighted (typical for GRN)

### Layer Structure

Each `Celltype` has two layers in multilayer network:
1. **receptor layer**: Receptor-receptor interactions (often sparse/fake)
2. **grn layer**: Gene regulatory network (TF→gene edges)
3. **bipartite**: Receptor→GRN connections

`Multicell` adds:
4. **cell_communication layer**: Ligand-receptor pairs across cell types

### HuMMuS Layer Roles

ReCoN uses the **HuMMuS methodology** (Trimbour et al., 2024) for GRN inference, which structures regulation across multiple biological layers:

1. **TF layer**: Transcription factors that regulate gene expression
2. **DNA layer**: Regulatory regions (ATAC-seq peaks) where TFs bind
3. **Gene layer**: Target genes regulated by TFs through DNA binding

The multilayer structure captures:
- **TF → DNA links**: From motif scanning (CellOracle) - where TFs can bind
- **DNA → Gene links**: From peak-gene associations (CIRCE) - which peaks regulate which genes
- **TF → Gene links**: From expression correlation (GRNBoost2/Arboreto) - overall regulatory relationships

Random walk across these layers integratesindirect (TF→DNA→Gene) regulatory evidence, improving GRN quality compared to single-layer approaches. It also allows walks within layers (e.g., TF-TF interactions) to capture co-regulation.

### Transition Matrices (lamb & eta)

- **lamb**: Layer-to-layer transition probability matrix (DataFrame)
  - Format: `lamb[i, j]` = probability of transition **FROM layer i TO layer j**
  - Default for Celltype: Allows transitions within GRN layer and from receptor→GRN
  - Example customization:
    ```python
    # Block transitions from GRN to receptor layer
    lamb.loc["Tcell_grn", "Tcell_receptor"] = 0
    # Increase self-loop probability in receptor layer
    lamb.loc["Tcell_receptor", "Tcell_receptor"] = 0.5
    ```
  - Rows must sum to =1 (restart is considered outside this matrix)
  These probabilties govern how random walks traverse between layers, and can be changed dpeending on biological assumptions and questions.
- **eta**: Restart probability per layer (Series)
  - Controls where RWR "restarts" (typically uniform across layers)
  - Must sum to 1.0 across all layers

## Key Dependencies

- **hummuspy**: Multilayer network RWR engine (custom fork: `cantinilab/HuMMuS`)
- **multixrank**: Core RWR algorithm implementation
- **celloracle**: TF-DNA motif scanning (use `cantinilab/celloracle@lite` branch)
- **arboreto**: GRNBoost2 for TF-gene inference
- **circe-py**: ATAC-seq peak processing
- **liana**: Cell-cell communication inference (external, used in examples)

## Development Workflows

### Installation

```bash
# Standard installation (recommended)
conda create -n recon python=3.10
conda activate recon
pip install recon[grn-lite]

# Without GRN dependencies (Python 3.8-3.12)
pip install recon
```

**Python version constraints**:
- **Minimum**: Python ≥3.8 (for core ReCoN functionality)
- **GRN inference**: Python ≥3.10 AND ≤3.10 (tight constraint)
  - `circe-py` (ATAC-seq peak-gene links) requires Python ≥3.10
  - `celloracle` + `gimmemotifs` (TF-DNA motif scanning) require Python ≤3.10
  - **Result**: Must use exactly Python 3.10 for full GRN inference with ATAC-seq
- **Python 3.11+**: Works for core ReCoN (no GRN inference) but celloracle not compatible

**Common installation issues**:

1. **llvmlite build failure** (macOS): numba→llvmlite requires LLVM
   ```bash
   # Solution: Install via conda (includes prebuilt binaries)
   conda install -c conda-forge llvmlite numba
   # Then retry: pip install recon
   ```

2. **CellOracle dependency**: GRN inference requires older dependencies → use custom CellOracle fork from cantinilab

### Building Documentation

Uses Sphinx with reStructuredText:
```bash
cd docs
make html  # or make.bat on Windows
```

Examples are Jupyter notebooks in `docs/source/recon_examples/`.

### Common Patterns

**Creating a Celltype object:**
```python
ct = recon.explore.Celltype(
    celltype_name="Tcell",
    grn_graph=grn_df,  # source, target, weight
    receptor_grn_bipartite=receptor_df,  # source (receptor), target (gene)
    receptor_graph=None  # often None, creates fake receptor
)
```

**Creating a Multicell object:**
```python
# From Celltype objects
mc = recon.explore.Multicell(
    celltypes=[ct_a, ct_b],  # or dict: {"NewName": ct_a, ...}
    cell_communication_graph=ccn_df  # columns: source, target, celltype_source, celltype_target, lr_means
)

# From dictionaries (celltypes created internally)
mc = recon.explore.Multicell(
    celltypes=[
        {"celltype_name": "Tcell", "grn_graph": grn_df, "receptor_grn_bipartite": rec_df},
        {"celltype_name": "Bcell", "grn_graph": grn_df2, "receptor_grn_bipartite": rec_df2}
    ],
    cell_communication_graph=ccn_df
)
```

**Running RWR:**
```python
ct.seeds = ["GENE1", "GENE2"]  # or dict with weights
multilayer = ct.Multixrank(restart_proba=0.7)
scores = multilayer.extract_multilayer_pagerank_scores()
```

**Treatment effect prediction:**
- Direct effect: α weight on receptors
- Indirect effect: (1-α) weight on other cell types' ligands
- Typical α = 0.8 (empirically validated)

## Testing & Quality

### Running Tests

```bash
pytest tests/              # Run all tests
pytest tests/ -v          # Verbose output
pytest tests/test_celltype.py  # Specific file
```

**Optional dependencies**: Some tests (e.g., `test_infer_grn.py` ATAC-seq tests) require `celloracle` + reference genomes and are skipped when unavailable:
```bash
# Install celloracle
pip install 'git+https://github.com/cantinilab/celloracle@lite'

# Install reference genome (mm10 for mouse)
pip install genomepy
genomepy install mm10 UCSC --annotation  # Stored in ~/.local/share/genomes/

# Now previously-skipped tests will run
pytest tests/test_infer_grn.py -v
```

### Test Structure
- **tests/test_celltype.py**: Celltype construction, graph types, seeds, renaming
- **tests/test_multicell.py**: Multicell integration, cell communication, bipartites
- **tests/conftest.py**: Shared fixtures (simple_grn, simple_receptor_grn, simple_cell_communication)
- **tests/README.md**: Detailed testing guide

### Testing Patterns
Tests focus on data structure correctness without mocking:
- **Network construction**: Verify DataFrame formats, multiplexes/bipartites structure
- **Node naming**: Validate `::celltype` suffix, `_receptor` suffix, `-celltype` for ligands
- **Graph types**: Test encoding matches DataFrame properties (e.g., "01" = undirected weighted)
- **Edge cases**: Missing receptor_graph, various seed formats (list/dict), renaming
- **Small synthetic data**: 3-5 node networks for fast deterministic tests

### What NOT to Test
- Full RWR execution (requires integration testing)
- Large realistic networks (slow)
- External data files (use in-memory fixtures)

### Key Test Scenarios
1. Celltype/Multicell creation with minimal valid inputs
2. Seed format validation (list vs dict)
3. lamb/eta matrix shape and normalization
4. Node renaming propagation through multiplexes/bipartites
5. Column renaming (receptor_grn_bipartite: source/target → col2/col1)

## File Organization

```
ReCoN/
├── src/recon/           # Main package
│   ├── explore/         # Core Celltype/Multicell classes
│   ├── infer_grn/       # GRN inference (HuMMuS)
│   └── plot/            # Visualization tools
├── docs/                # Sphinx documentation
│   └── source/
│       ├── recon_examples/    # Tutorial notebooks
│       └── recon_explained/   # Conceptual guides
└── figures/             # README/paper figures
```

## Important Notes

- **Alpha parameter**: For treatment effects, `α=0.8` (direct) vs `1-α=0.2` (indirect) is the validated default
- **Memory considerations**: Set `copy_graphs=False` in constructors for large networks to avoid duplication
- **Node suffixes**: The `::celltype` suffix is critical for `Multicell` - don't strip it manually
- **Column renaming**: receptor_grn_bipartite internally renames `source→col2`, `target→col1` for multixrank compatibility

## Common Pitfalls

1. Mixing graph types: Ensure `graph_type` strings match actual DataFrame properties
2. Missing receptor_graph: Always OK to pass `None` - fake receptor created automatically  
3. Seeds format: Can be list `["GENE1"]` or dict `{"GENE1": 0.5}` but not mixed
4. Sparse warnings: `hummuspy` may warn about graph sparsity - usually safe to ignore

## Packaging & Distribution

### PyPI Publishing Workflow

ReCoN uses GitHub Actions (`.github/workflows/wheels.yml`) for automated PyPI publishing:
- **Trigger**: Automatically publishes on GitHub releases
- **Build**: Creates universal wheel (`py3-none-any`) and source distribution
- **Authentication**: Uses PyPI API token (stored in `secrets.pypi_password`)

```yaml
# Key workflow steps:
python -m build              # Build distributions
twine check dist/*           # Validate metadata
gh-action-pypi-publish       # Upload to PyPI
```

### Dependency Management Constraint

**Critical**: PyPI does not allow direct Git dependencies. The project uses `celloracle-lite>=0.21.0` (PyPI package) instead of Git URLs.

Current `pyproject.toml` structure:
```toml
[project.optional-dependencies]
grn-lite = [
    "celloracle-lite>=0.21.0"
]
# grn = [
#   "celloracle @ git+https://github.com/cantinilab/celloracle@lite"
# ]

[project.urls]
Homepage = "https://recon.readthedocs.io"
Repository = "https://github.com/cantinilab/recon"
```

**Status**: PyPI now accepts `celloracle-lite` as a PyPI package dependency (changed from Git URL in v0.1.0). The commented-out `grn` extras with Git URL is deprecated.

## Citations

When modifying code, maintain compatibility with:
- Trimbour et al. 2025 (ReCoN paper): bioRxiv 2025.09.20.461080
- Trimbour et al. 2024 (HuMMuS): Bioinformatics 40(3), btae143
