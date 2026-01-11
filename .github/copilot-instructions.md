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
- Node naming convention for multicell: `{gene_name}::{celltype_name}` to avoid conflicts
- Receptor nodes: `{receptor_name}_receptor` suffix
- Cell communication nodes: `{ligand_name}-{celltype_source}` format

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

### Transition Matrices (lamb & eta)

- **lamb**: Layer-to-layer transition probability matrix (DataFrame)
  - Default: Allow transitions within GRN and receptor→GRN
- **eta**: Restart probability per layer (Series)
  - Controls where RWR "restarts" (typically uniform across layers)

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
# Standard installation
conda create -n recon python=3.10
conda activate recon
pip install recon[grn]

# Without GRN dependencies (newer Python versions)
pip install recon
```

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

### PyPI Publishing Constraint

**Critical**: PyPI does not allow direct Git dependencies. The `celloracle@git+https://github.com/cantinilab/celloracle@lite` dependency in `[grn]` extras causes upload failures.

Current error from wheels.yml GitHub Action:
```
ERROR HTTPError: 400 Bad Request from https://upload.pypi.org/legacy/
Can't have direct dependency: celloracle@ git+https://github.com/cantinilab/celloracle@lite
```

**Workarounds**:
1. Document celloracle installation separately in README (current approach)
2. Users install via: `pip install recon` then manually `pip install git+https://...`
3. Consider contributing celloracle-lite as separate PyPI package
4. Use Trusted Publisher workflow (no API tokens) - see warning in GitHub Actions logs

When modifying dependencies, **never add Git URLs to required or optional dependencies** if PyPI publishing is needed.

## Citations

When modifying code, maintain compatibility with:
- Trimbour et al. 2025 (ReCoN paper): bioRxiv 2025.09.20.461080
- Trimbour et al. 2024 (HuMMuS): Bioinformatics 40(3), btae143
