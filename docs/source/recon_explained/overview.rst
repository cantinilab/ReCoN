ReCoN overview and algorithm
================================

What ReCoN builds
------------------

ReCoN assembles a **heterogeneous multilayer network** for the multicellular system under study.
Each cell type contributes its own intracellular layers — a gene regulatory network (GRN) layer and a
receptor layer — and cell types are connected to one another through a **ligand-receptor bipartite layer**,
so that a ligand produced in one cell type can act on the receptors of another (or the same) cell type.

The intracellular layers can be inferred from the user's own single-cell RNA-seq (and optionally
scATAC-seq) data. The receptor-to-gene links, however, come from a fixed, pre-compiled prior knowledge
resource shared across all analyses of a given species — see
:doc:`Receptor Gene Databases <receptor_gene_db>` for details on why this part of the network is *not*
recomputed per dataset.

Direct and indirect effects
----------------------------

Given a set of **seed** nodes (e.g. the targets of a treatment, or a set of differentially expressed
genes), ReCoN estimates how their perturbation propagates through the multicellular network using two
separate **Random Walk with Restart (RWR)** processes:

- **Direct effect**: an RWR that starts from the seeds and walks through the receptors of the cell type
  of interest. It captures the effect of treatment-receptor binding local to that cell type.
- **Indirect effect**: an RWR that starts from the genes of the *other* cell types, propagates to the
  ligands they produce, and continues into the receptors of the cell type of interest (and other cell
  types). It captures how surrounding cells respond and, in turn, feed back on the focal cell type
  through cell-cell communication.

The two RWR outputs are combined into a single score per gene and cell type:

.. math::

   \text{score} = (1 - \alpha) \cdot \text{direct effect} + \alpha \cdot \text{indirect effect}

:math:`\alpha \in [0, 1]` is the weight given to the **indirect** effect, and :math:`1 - \alpha` is the
weight given to the **direct** effect (default :math:`\alpha = 0.8`). See
:doc:`get_ready` for guidance on choosing this parameter, and the main
:doc:`../index` page for why an indirect-effect-dominant setting performed best in our benchmarks.

What the output scores mean
----------------------------

ReCoN outputs one RWR-derived score per gene (and, in multicellular runs, per target cell type). These
scores should be read as:

- **Relative, not absolute**: scores reflect how much a gene is reached by the propagated perturbation
  compared to other genes in the *same* network — they are not calibrated probabilities and should not be
  compared across networks of different sizes or densities.
- **Higher means more affected**: a higher score indicates the gene is more strongly reached by the
  simulated perturbation, either directly (through receptor binding) or indirectly (through cell
  communication feedback), depending on how ``alpha`` was set.
- **Zero is common and expected**: genes disconnected from the seeds in the network, or reached only
  through very low-weight edges, will score at or near zero. This does not necessarily mean the gene is
  biologically unaffected — only that the network used does not support a propagation path to it.
- **Direct vs. indirect breakdown**: when combining effects, the direct and indirect components remain
  available separately (before being merged by ``alpha``), so you can inspect which part of the network
  is driving a gene's combined score.

For a more detailed, use-case-oriented walkthrough of score interpretation (comparing conditions,
validating predictions, reading multicellular output), see the
:doc:`Troubleshooting and FAQs page <get_ready>`.
