# ReCoN Test Suite

## Running Tests

```bash
# Run all tests
pytest tests/

# Run with verbose output
pytest tests/ -v

# Run with coverage report
pytest tests/ --cov=src/recon --cov-report=term-missing

# Run a specific test file
pytest tests/test_celltype.py

# Run a specific test class
pytest tests/test_celltype.py::TestCelltypeConstruction

# Run a specific test
pytest tests/test_celltype.py::TestCelltypeConstruction::test_create_celltype_basic
```

### Optional Tests

Some tests require optional dependencies (like celloracle) and will be automatically skipped when those dependencies are not installed. These tests use `@pytest.mark.skipif` to gracefully skip when the dependency is unavailable.

#### GRN Inference Tests (test_infer_grn.py)

**Requirements:**
- Basic tests: No extra dependencies (always run)
- ATAC-seq tests: Require celloracle + reference genome

**To run all GRN tests including ATAC-seq:**

```bash
# 1. Install celloracle
pip install 'git+https://github.com/cantinilab/celloracle@lite'

# 2. Install reference genome (mm10 for mouse)
pip install genomepy
genomepy install mm10 UCSC --annotation

# 3. Run tests (previously skipped tests will now run)
pytest tests/test_infer_grn.py -v
```

**Where are genomes stored?**

Genomes are downloaded to `~/.local/share/genomes/` (user-specific, not portable across machines). Each developer/CI environment needs to download genomes separately. Available genomes:

```bash
# List available genomes
genomepy list

# Install specific genomes
genomepy install hg38 UCSC --annotation  # Human
genomepy install mm10 UCSC --annotation  # Mouse
genomepy install mm39 UCSC --annotation  # Mouse (latest)

# Check installed genomes
genomepy genomes
```

**Why tests are skipped:**

- `celloracle not installed` → Install with `pip install 'git+https://github.com/cantinilab/celloracle@lite'`
- `mm10 genome not installed` → Install with `genomepy install mm10 UCSC --annotation`
- Both conditions must be met for ATAC-seq tests to run

## Test Structure

### `test_celltype.py`
Tests for the `Celltype` class covering:
- Object construction with various graph configurations
- Graph type encoding (directed/undirected, weighted/unweighted)
- Default lamb/eta matrix generation
- Seed handling (list vs dict formats)
- Celltype renaming
- DataFrame format preservation and column renaming

### `test_multicell.py`
Tests for the `Multicell` class covering:
- Integration of multiple celltypes
- Cell communication layer creation
- Node suffix addition (`::celltype`)
- Bipartite connections (ligands and receptors)
- Transition matrix (lamb) for multicellular systems
- Celltype renaming in multicellular context

### `test_infer_grn.py`
Tests for GRN inference functions covering:
- TF network creation (`compute_tf_network`)
- RNA network inference with GRNBoost2-style methods (`compute_rna_network`)
- **Optional**: ATAC-seq TF-to-peak links (requires celloracle, auto-skipped if not installed)
- Method selection (GBM vs RF)
- Input validation (DataFrame vs AnnData)

### `test_lambda.py`
Tests for transition matrix generation:
- Downstream vs upstream direction handling
- Intracell vs intercell transition strategies
- Row normalization verification

### `test_plot.py`
Tests for visualization functions:
- `illustrate_multicell` input validation
- Plot generation with matplotlib integration

### `test_utils.py`
Tests for utility functions:
- `split_layer_name` node name parsing
- Separator handling

### `conftest.py`
Shared fixtures:
- `simple_grn`: 3-node gene regulatory network
- `simple_receptor_grn`: 2-receptor to gene connections
- `simple_cell_communication`: 2-ligand/receptor pairs across 2 cell types

## Test Philosophy

Tests focus on:
1. **Core data structures**: Verify multiplexes, bipartites, lamb, eta are correctly constructed
2. **Format compliance**: Ensure DataFrames have expected columns and node naming conventions
3. **Edge cases**: Test with/without receptor graphs, various seed formats, celltype renaming
4. **No mocking**: Tests use actual Celltype/Multicell objects with minimal synthetic data

Tests intentionally avoid:
- Full RWR execution (requires hummuspy integration testing)
- Large realistic networks (keep tests fast)
- External data files (use in-memory fixtures)

## Adding New Tests

When adding features, add tests that verify:
1. Data structure correctness (keys in dicts, columns in DataFrames)
2. Node naming conventions (`::celltype` suffixes, `_receptor` suffixes, `-celltype` for ligands)
3. Graph type encoding matches actual graph properties
4. Error handling for invalid inputs

## Notes

- CellOracle is optional - tests run without it (import wrapped in try/except)
- Tests generate UserWarnings about missing receptor_graph - this is expected behavior
- scipy.sparse deprecation warnings from multixrank are not our concern
