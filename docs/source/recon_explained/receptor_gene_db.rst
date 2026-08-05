Receptor Gene Databases
=======================

The receptor-gene network is what connects a cell's receptors to its intracellular gene regulatory
network (GRN): once a ligand binds a receptor, this network tells ReCoN which intracellular genes are
downstream of that signaling event, so the effect can keep propagating through the GRN.

.. important::
   This receptor-gene network is **not computed from your own dataset**. It is a fixed, pre-compiled
   **prior knowledge resource**, shipped with ReCoN and derived from NicheNet's prior knowledge network
   (PKN). It only varies by **species** (human or mouse) — the same table is reused across every
   analysis, tutorial, and dataset for that species, exactly like a TF motif database or a pathway
   annotation resource. You do not need to rebuild it for each new study; you only choose which
   species' version to load.

You can also provide your **own** receptor-gene network instead, for example if:

- you work with a species other than human or mouse,
- you want receptor-gene links tailored to your own dataset rather than the generic prior,
- or you want a different receptor-gene network per cell type (technically supported, though we have
  not tested this setting).

Available Databases
-------------------

- **Human Receptor Genes**: Derived from NicheNet's PKN.
- **Mouse Receptor Genes**: Derived from NicheNet's PKN.

Usage
-----

To load receptor gene databases in ReCoN, use the `load_receptor_genes` function:

.. code-block:: python

   from recon.data import load_receptor_genes

   # Load human receptor genes
   human_receptors = load_receptor_genes("human_receptor_gene_from_NichenetPKN")

   # Load mouse receptor genes
   mouse_receptors = load_receptor_genes("mouse_receptor_gene_from_NichenetPKN")

These functions return a pandas DataFrame containing receptor gene information. Load it once per
species and reuse it across all cell types and analyses — there is no dataset-specific step here.

Structure: a receptor-to-gene bipartite network
-------------------------------------------------

Each database is a **bipartite edge table** with three columns:

- ``source``: a receptor gene
- ``target``: an intracellular gene regulated downstream of that receptor's signaling
- ``weight``: the strength of that receptor-to-gene link, derived from NicheNet's regression-based
  ligand-receptor-to-target-gene model (edges below a fixed weight threshold are already filtered out)

"Bipartite" here means edges only ever go from a receptor to a downstream gene — never between two
receptors or between two genes. This is what lets ReCoN plug a common, receptor-centric layer in front
of *any* cell type's GRN: the same receptor-gene prior is reused for every cell type of a given species,
while the GRN layer itself (target genes and their regulators) stays specific to each cell type and is
inferred from your own data.

Integration in Tutorials
-------------------------

Receptor gene databases are used in the following tutorials:

1. **Predicting Treatment Effects**: Identifying receptor-ligand interactions.
2. **Building GRNs with HuMMuS**: Using receptor genes to infer gene regulatory networks.

For more details, refer to the respective tutorials in the documentation.
