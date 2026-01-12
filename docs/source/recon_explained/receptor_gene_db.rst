Receptor Gene Databases
=======================

Receptor gene databases are essential for understanding cell communication networks. ReCoN provides pre-compiled receptor gene lists for human and mouse, derived from NicheNet's prior knowledge network (PKN). These databases are used to identify receptor-ligand interactions in single-cell RNA-seq data.

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

These functions return a pandas DataFrame containing receptor gene information.

Integration in Tutorials
-------------------------

Receptor gene databases are used in the following tutorials:

1. **Predicting Treatment Effects**: Identifying receptor-ligand interactions.
2. **Building GRNs with HuMMuS**: Using receptor genes to infer gene regulatory networks.

For more details, refer to the respective tutorials in the documentation.