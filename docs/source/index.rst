.. role:: raw-html(raw)
    :format: html
    
.. ReCoN documentation master file, created by
   sphinx-quickstart on Mon Sep 22 18:07:57 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

.. image:: ../../figures/ReCoN_logo_1000.png
   :alt: ReCoN-logo Remi-Trimbour 2025
   :align: center


.. toctree::
   :maxdepth: 1
   :hidden:
   :glob:
   :caption: ReCoN explained

   recon_explained/*

.. toctree::
   :maxdepth: 1
   :hidden:
   :glob:
   :caption: Examples

   recon_examples/*

.. toctree::
   :maxdepth: 1
   :hidden:
   :glob:
   :caption: API Reference

   API/*

Welcome to ReCoN's documentation!
==================================
**ReCoN** is a new tool for reconstructing multicellular models.

It combines both **gene regulatory networks** and **cell communication networks** to explore the molecular coordinations between multiple cell types — all at once.

ReCoN uses **heterogeneous multilayer networks** and integrates several layers of information into a **complex network**, ready to be explored and analyzed.  
Both the GRNs and intercellular networks are inferred from **single-cell RNA-seq data** (and optionally **scATAC-seq**).

.. image:: ../../figures/recon_abstract_1.png
   :alt: ReCoN-abstract Remi-Trimbour 2025
   :align: center

.. note::
   💡 The philosophy behind ReCoN:
   :raw-html:`<br />` 🧬 *Cells do not act in isolation, but in a coordinated, dynamic system.*

ReCoN can be used to address several biological questions, including:

.. image:: ../../figures/recon_outputs.png
   :alt: ReCoN-outputs Remi-Trimbour 2025
   :align: center
   :scale: 50%

:raw-html:`<br />`

.. admonition:: ReCoN use cases
   :class: tip

   - :doc:`Predicting treatment effects in multicellular systems <recon_examples/1.recon_molecular_treatment>`
   - :doc:`Understanding multicellular program coordination <recon_examples/2.recon_multicellular_coordination>`
   - :doc:`Exploring intracellular and intercellular regulatory mechanisms <recon_examples/3.recon_molecular_cascades>`
   - :doc:`Build GRNs through HuMMuS methodology <recon_examples/4.recon_hummus>`


📦 Installation
---------------

ReCoN is available as a Python package and can be installed through pip. 
⚠️ *Actually, not yet but very very soon ! for now, check the github :)*

..  code-block:: python

   conda create -n recon python=3.10
   conda activate recon
   pip install recon[grn-lite]

If you are generating your grn externally, you can install ReCoN without the GRN dependencies.:raw-html:`<br />`
*You should then be able to use more recent version of Python.*

..  code-block:: python

   pip install recon

⚠️ **To generate GRNs**, ReCoN requires **CellOracle** to be installed.
Since CellOracle requires quite old dependency versions, we propose to install our
`own lite branch <https://github.com/cantinilab/CellOracle>`_ that contains 
only the necessary functions through the code "recon[grn-lite]".

.. admonition:: Installation troubleshooting & Frequently asked questions
   :class: warning

   See the :doc:`Troubleshooting and FAQs page <recon_explained/get_ready>`


💊 How does a treatment affect the molecular state of multiple cell types?
---------------------------------------------------------------------------

ReCoN can be used to predict how a treatment (e.g., a drug) affects
the molecular state of each cell type in a multicellular context (e.g., organ, tumor microenvironment).
:raw-html:`<br />`
It represents both **direct** effect, through treatment - receptor bindings, and **indirect** effects of the treatment,
through **cell communication** interplays.

.. image:: ../../figures/indirect_effect_schema.png
   :alt: ReCoN-indirect-effect Remi-Trimbour 2025
   :align: center
   :scale: 50%


.. note::
   A treatment effect can be decomposed into two components:

   - **Direct effect** — The effect of the treatment caused by *direct* binding of the receptors of a cell type.

   - **Indirect effect** — The effect on a cell type **mediated by other cell types** that respond and, in turn, secrete ligands which modulate the molecular state of the cell type of interest.


ReCoN models these two components through different random walk with restart (RWR) processes
on the multicellular network (cf. :doc:`ReCoN overview and algorithm <recon_explained/overview>`). 
The parameter 
:math:`\alpha \in [0, 1]` allows to tune the relative importance of the direct and indirect effects.
:math:`\alpha` is the weight of the direct effect, while :math:`1 - \alpha` is the weight of the indirect effect.

.. image:: ../../figures/indirect_direct_effect_formula.png
   :alt: ReCoN-direct-indirect-effect-formula Remi-Trimbour 2025
   :align: center

**Why indirect effects matter**

.. admonition:: Importance of indirect effects
   :class: important

   Surrounding cells can secrete ligands in response to the treatment, which then feed back
   and alter signaling and regulation in the focal cell type.

   In our evaluation of ReCoN, we found that **giving more importance to the indirect effect**
   (:math:`\alpha = 0.8` — or an indirect effect 4 times stronger than the direct effect) led to the best performance
   in both individual and multiple perturbation setups. :raw-html:`<br />`
   *(Trimbour et al., 2025 — Immune Dictionary and Heart Failure showcases)*


See how to use ReCoN to model the effect of a drug here : :doc:`Predicting treatment effects in multicellular systems <recon_examples/1.recon_molecular_treatment>`.

🧫 Understanding multicellular program coordination
----------------------------------------------------
How are the surrounding cells regulating and impacted by the state of a given cell-type ?
ReCoN can be used to explore the causes and consequences of a given cell/multicellular state.

We can then identify the key molecules and cell types that are involved in the coordination of this state.

.. image:: ../../figures/recon_multicellular_programs.png
   :alt: ReCoN-multicellular-programs Remi-Trimbour 2025
   :align: center
   :scale: 50%

See how to use ReCoN to explore multicellular coordination here : :doc:`Exploring multicellular coordination <recon_examples/2.recon_multicellular_coordination>`.


⚙️ Visualizing multicellular molecular cascades modulating cell states
-----------------------------------------------------------------------

ReCoN can reconstruct the intercellular cascades driving a specific transcriptomic state. 
These cascades includes intracellular elements (receptors, transcription factors), but also ligands
and their own regulators. It offers you a comprehensive view of these interactions,
and the possibility to identify new potential targets at different regulatory levels.

See how to use ReCoN to explore intracellular regulatory elements here : :doc:`Exploring molecular cascades and identify regulators <recon_examples/3.recon_molecular_cascades>`.

🧬 Building GRNs through HuMMuS methodology
--------------------------------------------

We previously developped HuMMuS (Trimbour et al., 2024), an other method based on heterogeneous multilayer networks to build gene regulatory networks from single-cell RNA-seq and ATAC-seq data.
HuMMuS can be used as a standalone method, but it was initially developed as an hybrid package with R and Python.

ReCoN can now also be used to run a Python implementation of HuMMuS to infer GRNs from single-cell data.
This implementation leverages the functions of CellOracle to build the prior knowledge links between TF, DNA regions and target genes.
HuMMuS is then applied on a multilayer composed of a TF layer, a DNA region layer and a target gene layer to infer the final GRN.

See how to use ReCoN to build GRNs with HuMMuS here : :doc:`Building gene regulatory networks with HuMMuS <recon_examples/4.recon_hummus>`.


📖 Cite ReCoN
----------------

If you use ReCoN in your work, please cite:

.. admonition:: Cite ReCoN
   :class: seealso

   Trimbour R., Ramirez Flores R. O., Saez Rodriguez J., Cantini L. (2025). ReCoN: Reconstructing multicellular models by integrating gene regulatory and cell communication networks. *bioRxiv*. https://doi.org/___________

If you also use ReCoN to generate GRNs, please cite:

.. admonition:: Cite HuMMuS
   :class: seealso

   Trimbour R., Ramirez Flores R. O., Saez Rodriguez J., Cantini L. (2025). ReCoN: Reconstructing multicellular models by integrating gene regulatory and cell communication networks. *bioRxiv*. https://doi.org/___________
  
   Trimbour R., Deutschmann I. M., Cantini L. (2024). HuMMuS: Inferring gene regulatory networks through heterogeneous multilayer networks. *Bioinformatics*, 40(3), btae143. https://doi.org/10.1093/bioinformatics/btae143

.. note::
   To download the tutorial data, use the following commands:

   **All tutorial data:**
   
   .. code-block:: python

      from recon.data import fetch_all_tutorial_data
      fetch_all_tutorial_data(data_dir='./data')

   **Specific file (e.g., RNA data):**
   
   .. code-block:: python

      from recon.data import fetch_tutorial_data
      fetch_tutorial_data('perturbation_tuto/rna.h5ad', data_dir='./data')

