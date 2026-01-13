Installation troubleshooting & FAQ
====================================

Installation
--------------

Basic installation
~~~~~~~~~~~~~~~~~~

For basic ReCoN functionality (without GRN inference from ATAC-seq):

..  code-block:: bash

    pip install recon

Or install from source for development:

..  code-block:: bash

    git clone https://github.com/cantinilab/recon.git
    cd recon
    pip install -e .

Installation with GRN inference (optional)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

You'll need gimmemotifs<=0.17.2, celloracle (lite branch), and llvmlite to install ReCoN with full GRN inference capabilities.

..  code-block:: bash

    # Create environment with required dependencies
    conda create -n recon  python=3.10
    conda activate recon
    
    # Install ReCoN with GRN extras
    pip install recon[grn-lite]


Installation with macOS
~~~~~~~~~~~~~~~~~~~~~~~

Some packages may be tricky to install with pip on macOS due to system library dependencies. We recommend using conda to manage these dependencies: gimmemotifs and llvmlite

..  code-block:: bash

    # Create environment with required dependencies
    conda create -n recon -c bioconda -c conda-forge python=3.10 gimmemotifs llvmlite cmake

    conda activate recon
    
    # Install ReCoN with GRN extras
    pip install recon[grn-lite]



You'll need cmake, gimmemotifs, and llvmlite to install ReCoN with full GRN inference capabilities including ATAC-seq motif scanning.

..  code-block:: bash

    # Create environment with required dependencies
    conda create -n recon -c bioconda -c conda-forge python=3.10 gimmemotifs llvmlite cmake
    conda activate recon
    
    # Install ReCoN with GRN extras
    pip install recon[grn-lite]

**Why these dependencies?**

- ``gimmemotifs``: TF motif scanning for ATAC-seq peaks (requires pre-compiled binaries from conda)
- ``llvmlite``: Required by numba for JIT compilation (system LLVM libraries needed on macOS)
- ``cmake``: Build tool for compiling C/C++ extensions

Installing reference genomes for ATAC-seq
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

If you plan to use ATAC-seq data for TF-to-peak motif scanning, you need to install reference genomes:

..  code-block:: bash

    # Install genomepy (included with celloracle)
    pip install genomepy
    
    # Install mouse genome (mm10)
    genomepy install mm10 -p UCSC -a
    
    # Install human genome (hg38)
    genomepy install hg38 -p UCSC -a
    
    # Check installed genomes
    genomepy genomes
    
    # List available genomes
    genomepy search mouse

**Where are genomes stored?**

Genomes are downloaded to ``~/.local/share/genomes/`` by default. 
There are usually **large files**: Genomes are typically 1-3 GB each.

You can customize the location:

..  code-block:: bash

    # Install to custom directory
    genomepy install mm10 -p UCSC -a -g /path/to/genomes
    
    # Or set environment variable
    export GENOMES_DIR=/path/to/genomes
    genomepy install mm10 -p UCSC -a

**Available genome providers:**

- ``UCSC``: University of California Santa Cruz (recommended)
- ``Ensembl``: European Bioinformatics Institute
- ``NCBI``: National Center for Biotechnology Information

**Common genomes:**

- ``mm10``: Mouse (GRCm38/mm10)
- ``mm39``: Mouse (GRCm39, latest)
- ``hg38``: Human (GRCh38)
- ``hg19``: Human (GRCh37, older)

Common installation issues
~~~~~~~~~~~~~~~~~~~~~~~~~~~

**Problem: llvmlite build fails on macOS**

..  code-block:: bash

    # Install via conda instead of pip
    conda install -c conda-forge llvmlite numba

**Problem: gimmemotifs installation fails**

Gimmemotifs has C dependencies that may not compile on all systems. Install from conda-forge and bioconda channels:

..  code-block:: bash

    conda install -c conda-forge -c bioconda gimmemotifs

**Problem: genomepy install command fails with "Got unexpected extra argument"**

The correct syntax uses ``-p`` flag for provider:

..  code-block:: bash

    # WRONG: genomepy install mm10 UCSC --annotation
    # RIGHT:
    genomepy install mm10 -p UCSC -a

**Problem: "Genomes_dir does not exist" error**

Create the directory first or specify a custom location:

..  code-block:: bash

    # Option 1: Create default directory
    mkdir -p ~/.local/share/genomes
    
    # Option 2: Use custom directory
    genomepy install mm10 -p UCSC -a -g /path/to/genomes

**Problem: ATAC tests are skipped with "mm10 genome not installed"**

This is expected when the genome isn't downloaded. To run ATAC-seq tests:

..  code-block:: bash

    # Install celloracle
    pip install 'git+https://github.com/cantinilab/celloracle@lite'
    
    # Install genome
    genomepy install mm10 -p UCSC -a
    
    # Verify installation
    ls ~/.local/share/genomes/mm10/
    
    # Run tests
    pytest tests/test_infer_grn.py -v

**Problem: I cannot compute GRNs**

CellOracle is required for GRN inference with ATAC-seq data. Options:

1. Install our 'lite' branch direclty with recon: ``pip install recon[grn-lite]``
2. Install it separately: ``pip install 'git+https://github.com/cantinilab/celloracle@lite'``
3. Compute your GRN externally and provide it to ReCoN.

**Problem: Tests are skipped for celloracle functions**

This is expected behavior when celloracle is not installed. The tests use ``@pytest.mark.skipif`` to gracefully skip ATAC-seq tests when celloracle is unavailable. Install with ``[grn]`` extras to run all tests.

Python version compatibility
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- **Minimum**: Python 3.8
- **Recommended**: Python 3.10+
- **Note**: Some dependencies (like circe-py) use Python 3.10+ type syntax (``Type | None``). If you encounter ``TypeError: unsupported operand type(s) for |``, upgrade to Python 3.10+.

GRN inference
--------------

What GRN inference methods are available?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ReCoN supports multiple approaches:

1. **TF-to-gene (RNA-seq only)**: Uses GRNBoost2-style gradient boosting (GBM) or Random Forest (RF)
2. **TF-to-gene with ATAC-seq**: Adds TF-to-peak motif scanning via CellOracle
3. **Receptor-to-gene**: Custom connections for cell surface receptors

When should I use ATAC-seq integration?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Use ATAC-seq when:

- You want more accurate TF-gene regulatory links
- You have matched scRNA-seq + scATAC-seq data
- You're interested in chromatin accessibility effects

Skip ATAC-seq when:

- You only have scRNA-seq data
- Computational resources are limited (motif scanning is slow)
- You're doing quick exploratory analysis

How accurate is GRN inference?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

GRN inference is probabilistic and noisy. ReCoN combines:

- **Expression correlation** (GRNBoost2 importance scores)
- **Motif evidence** (TF binding site predictions from ATAC peaks)
- **Network propagation** (RWR to capture indirect effects)

**Best practices:**

- Use biological validation (ChIP-seq, literature)
- Focus on highly-ranked edges (top 10-20%)
- Combine with perturbation data when available

Can I use my own GRN?
~~~~~~~~~~~~~~~~~~~~~~

Yes! Provide a custom DataFrame with columns ``['source', 'target', 'weight']``:

..  code-block:: python

    import pandas as pd
    from recon.explore import Celltype
    
    # Custom GRN from literature or ChIP-seq
    custom_grn = pd.DataFrame({
        'source': ['TF1', 'TF1', 'TF2'],
        'target': ['GENE1', 'GENE2', 'GENE3'],
        'weight': [0.8, 0.6, 0.9]
    })
    
    celltype = Celltype(
        grn=custom_grn,
        receptor_grn=receptor_grn,
        name="MyCell"
    )

ReCoN results
--------------

What do the scores mean?
~~~~~~~~~~~~~~~~~~~~~~~~~

ReCoN outputs **Random Walk with Restart (RWR) scores** representing:

- **Treatment propagation**: How molecular perturbations flow through the network
- **Values 0-1**: Higher = more affected by treatment
- **Relative ranking**: Compare scores across genes/cells, not absolute magnitudes

The ``alpha`` parameter (default 0.8) controls:

- **High alpha (0.8-0.9)**: Treatment stays local to seeds
- **Low alpha (0.3-0.5)**: Treatment diffuses widely across network

What seeds should I use?
~~~~~~~~~~~~~~~~~~~~~~~~~

**Seeds** are your treatment entry points. Common choices:

1. **Differentially expressed genes** from treated vs control
2. **Drug targets** (e.g., receptor targeted by therapy)
3. **Pathway genes** (e.g., all genes in immune response pathway)

..  code-block:: python

    # Dictionary format: {gene: score}
    seeds = {'RECEPTOR1': 1.0, 'TF1': 0.8}
    
    # Or list format (all seeds weighted equally)
    seeds = ['RECEPTOR1', 'TF1', 'GENE1']

How do I interpret multicellular results?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In ``Multicell`` objects:

- **Node names**: Suffixed with ``::celltype`` (e.g., ``CD8::T_cell``)
- **Ligand-receptor connections**: ``LIGAND-celltype`` → ``RECEPTOR_receptor::celltype``
- **Cell communication layer**: Bipartite graphs between cell types
- **Lamb matrix**: Controls transition probabilities between layers

Higher scores in receiving cells indicate:

- Strong cell-cell communication effects
- Potential for coordinated responses
- Targets for combination therapies

Why are some scores zero?
~~~~~~~~~~~~~~~~~~~~~~~~~~

Possible reasons:

1. **Disconnected components**: Gene not reachable from seeds in network
2. **Low edge weights**: Weak connections filtered out
3. **High restart probability**: Treatment didn't diffuse far enough (try lower restart probability)
4. **Missing edges**: Incomplete GRN (add more regulatory links)

ReCoN interpretation
----------------------

How to validate ReCoN predictions?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Literature search**: Check if predicted genes are known treatment targets
2. **Pathway analysis**: Are high-scoring genes in expected pathways?
3. **Perturbation data**: Compare with experimental knockdown/overexpression
4. **Cross-validation**: Split cells into train/test, validate predictions
5. **Temporal data**: Do predictions match time-series gene expression?

What biological insights can ReCoN provide?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

ReCoN helps answer:

- **Which genes are affected** by a treatment beyond direct targets?
- **How do cell types coordinate** their responses?
- **What off-target effects** might occur?
- **Why do some cells respond** differently than others?
- **And maybe your own biological question!** :)

How to compare conditions?
~~~~~~~~~~~~~~~~~~~~~~~~~~~

Compare RWR scores between conditions:

..  code-block:: python

    # Run ReCoN on both conditions
    results_control = celltype.Multixrank(seeds=control_seeds, alpha=0.8)
    results_treated = celltype.Multixrank(seeds=treated_seeds, alpha=0.8)
    
    # Compare scores
    import pandas as pd
    comparison = pd.DataFrame({
        'control': results_control['GRN'],
        'treated': results_treated['GRN']
    })
    comparison['delta'] = comparison['treated'] - comparison['control']
    
    # Genes with largest changes
    top_changes = comparison.nlargest(20, 'delta')

ReCoN Visualization
----------------------

What visualization tools are available?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Sankey diagrams**: Trace treatment flow from seeds through network
2. **Network plots**: Visualize multicellular architecture
3. **Heatmaps**: Compare scores across cell types or conditions
4. **Custom plotting**: Export scores to pandas DataFrames for ggplot/matplotlib

Example Sankey diagram:

..  code-block:: python

    from recon.plot import plot_sankey
    
    plot_sankey(
        multicell=multicell,
        results=results,
        source_celltype='Tumor',
        target_celltype='T_cell',
        top_n=10
    )

How to export results?
~~~~~~~~~~~~~~~~~~~~~~

Results are pandas DataFrames - use standard methods:

..  code-block:: python

    # Save to CSV
    results['GRN'].to_csv('recon_results.csv')
    
    # Save to Excel with multiple sheets
    with pd.ExcelWriter('recon_results.xlsx') as writer:
        results['GRN'].to_excel(writer, sheet_name='GRN')
        results['Receptor'].to_excel(writer, sheet_name='Receptors')

Reproducibility
---------------

How to ensure reproducible results?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. **Set random seeds**: ReCoN uses deterministic algorithms, but upstream GRN inference may not
2. **Version dependencies**: Document versions of recon, scanpy, etc.
3. **Save parameters**: Record alpha, seeds, network sizes
4. **Archive networks**: Save GRN/receptor_grn DataFrames

..  code-block:: python

    # Save complete configuration
    import json
    
    config = {
        'recon_version': '0.1.0',
        'alpha': 0.8,
        'seeds': seeds,
        'graph_types': {'GRN': '01', 'Receptor': '01'},
        'n_genes': len(celltype.multiplexes['GRN'])
    }
    
    with open('recon_config.json', 'w') as f:
        json.dump(config, f, indent=2)

Performance & Scalability
-------------------------

How long does ReCoN take?
~~~~~~~~~~~~~~~~~~~~~~~~~~

**GRN inference** (slowest):

- GRNBoost2: 10-60 minutes for 5000 genes
- ATAC motif scanning: 1-10 hours depending on peak count

**RWR computation** (fast):

- Single celltype: Seconds to minutes
- Multicellular (3-5 celltypes): 1-5 minutes

**Tips for speed:**

- Pre-compute GRN once, reuse for multiple seed sets
- Limit GRN to top expressed genes (2000-5000)
- Use sparse matrices (automatically handled)

Memory requirements?
~~~~~~~~~~~~~~~~~~~~

- **Minimal**: ~2 GB for small networks (<1000 genes)
- **Typical**: 4-8 GB for realistic scRNA-seq data
- **Large**: 16+ GB for 10+ cell types with full GRNs

Reduce memory by:

- Filtering low-expressed genes before GRN inference
- Using fewer celltypes in multicellular models
- Clearing intermediate results: ``del results``

Can I parallelize ReCoN?
~~~~~~~~~~~~~~~~~~~~~~~~~

- **GRN inference**: Parallel by default (set ``n_cpu`` parameter)
- **RWR computation**: Single-threaded (already very fast)
- **Multiple conditions**: Run in parallel with multiprocessing

..  code-block:: python

    from multiprocessing import Pool
    
    def run_recon(seed_set):
        return celltype.Multixrank(seeds=seed_set, alpha=0.8)
    
    with Pool(4) as p:
        results = p.map(run_recon, [seeds1, seeds2, seeds3, seeds4])

Getting Help
------------

Where to find more information?
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

- **Documentation**: https://recon.readthedocs.io
- **GitHub Issues**: https://github.com/cantinilab/recon/issues
- **Examples**: See notebooks in ``docs/source/recon_examples/``
- **Paper**: [Add citation when published]

How to report bugs?
~~~~~~~~~~~~~~~~~~~

Open a GitHub issue with:

1. **Python/package versions**: ``pip list | grep recon``
2. **Minimal example**: Code that reproduces the error
3. **Error message**: Full traceback
4. **Expected behavior**: What should happen instead

..  code-block:: bash

    # Get version info for bug report
    python -c "import recon; print(recon.__version__)"
    python --version
    pip list | grep -E "recon|scanpy|numpy|pandas"

Contributing
~~~~~~~~~~~~

ReCoN is open source (GPL-3.0 license). Contributions welcome:

- **Code**: Submit pull requests on GitHub
- **Documentation**: Fix typos, add examples
- **Testing**: Report issues, suggest improvements
- **Citations**: Cite ReCoN in your publications

License: GPL-3.0 allows free use, modification, and distribution, but derivative works must also be open source.