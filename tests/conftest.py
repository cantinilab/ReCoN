"""Shared pytest fixtures for ReCoN tests."""
import pytest
import pandas as pd
import numpy as np


@pytest.fixture
def simple_grn():
    """Minimal 3-node GRN for testing."""
    return pd.DataFrame({
        'source': ['TF1', 'TF1', 'TF2'],
        'target': ['GENE1', 'GENE2', 'GENE2'],
        'weight': [0.8, 0.6, 0.7],
        'network_key': ['gene', 'gene', 'gene']
    })


@pytest.fixture
def simple_receptor_grn():
    """Minimal receptor-to-GRN bipartite links.
    
    Note: In ReCoN, this connects receptors to genes.
    Column names are 'source' (receptor) and 'target' (gene).
    """
    return pd.DataFrame({
        'source': ['RECEPTOR1', 'RECEPTOR2'],
        'target': ['GENE1', 'GENE2'],
        'score': [1.0, 1.0]
    })


@pytest.fixture
def simple_cell_communication():
    """Minimal cell-cell communication network.
    
    Format matches LIANA output:
    - source/target: ligand/receptor names
    - celltype_source/target: cell type producing/receiving
    - lr_means: interaction strength
    """
    return pd.DataFrame({
        'source': ['LIGAND1', 'LIGAND2'],
        'target': ['RECEPTOR1', 'RECEPTOR2'],
        'celltype_source': ['CellA', 'CellB'],
        'celltype_target': ['CellB', 'CellA'],
        'lr_means': [0.5, 0.6],
        'network_key': ['cell_communication', 'cell_communication']
    })
