"""Tests for recon.explore.Multicell class."""
import pytest
import pandas as pd
import numpy as np
from recon.explore.recon import Celltype, Multicell


class TestMulticellConstruction:
    """Test Multicell object creation and integration."""
    
    def test_create_multicell_from_celltypes(
        self, simple_grn, simple_receptor_grn, simple_cell_communication
    ):
        """Test Multicell creation from Celltype objects."""
        ct_a = Celltype(
            celltype_name="CellA",
            grn_graph=simple_grn.copy(),
            receptor_grn_bipartite=simple_receptor_grn.copy()
        )
        
        ct_b = Celltype(
            celltype_name="CellB",
            grn_graph=simple_grn.copy(),
            receptor_grn_bipartite=simple_receptor_grn.copy()
        )
        
        mc = Multicell(
            celltypes=[ct_a, ct_b],
            cell_communication_graph=simple_cell_communication
        )
        
        assert "cell_communication" in mc.multiplexes
        assert "CellA_grn" in mc.multiplexes
        assert "CellB_grn" in mc.multiplexes
        assert len(mc.celltypes_names) == 2
    
    def test_create_multicell_from_dicts(
        self, simple_grn, simple_receptor_grn, simple_cell_communication
    ):
        """Test Multicell creation from celltype dictionaries."""
        celltypes = [
            {
                "celltype_name": "CellA",
                "grn_graph": simple_grn.copy(),
                "receptor_grn_bipartite": simple_receptor_grn.copy()
            },
            {
                "celltype_name": "CellB",
                "grn_graph": simple_grn.copy(),
                "receptor_grn_bipartite": simple_receptor_grn.copy()
            }
        ]
        
        mc = Multicell(
            celltypes=celltypes,
            cell_communication_graph=simple_cell_communication
        )
        
        assert mc.celltypes_names == ["CellA", "CellB"]
    
    def test_node_suffix_addition(
        self, simple_grn, simple_receptor_grn, simple_cell_communication
    ):
        """Test that nodes get ::celltype suffix in Multicell."""
        ct_a = Celltype(
            celltype_name="CellA",
            grn_graph=simple_grn.copy(),
            receptor_grn_bipartite=simple_receptor_grn.copy()
        )
        
        ct_b = Celltype(
            celltype_name="CellB",
            grn_graph=simple_grn.copy(),
            receptor_grn_bipartite=simple_receptor_grn.copy()
        )
        
        mc = Multicell(
            celltypes=[ct_a, ct_b],
            cell_communication_graph=simple_cell_communication
        )
        
        # Check that GRN nodes have celltype suffix
        grn_layer = mc.multiplexes["CellA_grn"]["layers"][0]
        # At least one node should have ::CellA suffix
        assert any('::CellA' in str(node) for col in ['source', 'target'] 
                   for node in grn_layer[col].values if col in grn_layer.columns)
    
    def test_cell_communication_layer_format(
        self, simple_grn, simple_receptor_grn, simple_cell_communication
    ):
        """Test cell communication layer has correct format."""
        ct_a = Celltype(
            celltype_name="CellA",
            grn_graph=simple_grn.copy(),
            receptor_grn_bipartite=simple_receptor_grn.copy()
        )
        
        ct_b = Celltype(
            celltype_name="CellB",
            grn_graph=simple_grn.copy(),
            receptor_grn_bipartite=simple_receptor_grn.copy()
        )
        
        mc = Multicell(
            celltypes=[ct_a, ct_b],
            cell_communication_graph=simple_cell_communication
        )
        
        cc_layer = mc.multiplexes["cell_communication"]["layers"][0]
        
        # Should have ligand-celltype format: LIGAND1-CellA
        assert any('-' in str(node) for col in ['source', 'target'] 
                   for node in cc_layer[col].values if col in cc_layer.columns)
    
    def test_bipartite_ligand_connections(
        self, simple_grn, simple_receptor_grn, simple_cell_communication
    ):
        """Test bipartite connections from GRN to ligands."""
        ct_a = Celltype(
            celltype_name="CellA",
            grn_graph=simple_grn.copy(),
            receptor_grn_bipartite=simple_receptor_grn.copy()
        )
        
        mc = Multicell(
            celltypes=[ct_a],
            cell_communication_graph=simple_cell_communication
        )
        
        # Should have bipartite for ligands
        assert "CellA_to_ligands" in mc.bipartites
        ligand_bipartite = mc.bipartites["CellA_to_ligands"]
        assert ligand_bipartite["source"] == "cell_communication"
        assert ligand_bipartite["target"] == "CellA_grn"
    
    def test_bipartite_receptor_connections(
        self, simple_grn, simple_receptor_grn, simple_cell_communication
    ):
        """Test bipartite connections from receptors to cell communication."""
        ct_a = Celltype(
            celltype_name="CellA",
            grn_graph=simple_grn.copy(),
            receptor_grn_bipartite=simple_receptor_grn.copy()
        )
        
        mc = Multicell(
            celltypes=[ct_a],
            cell_communication_graph=simple_cell_communication
        )
        
        # Should have bipartite for receptors
        assert "CellA_to_receptors" in mc.bipartites
        receptor_bipartite = mc.bipartites["CellA_to_receptors"]
        assert receptor_bipartite["source"] == "cell_communication"
        assert receptor_bipartite["target"] == "CellA_receptor"


class TestMulticellLambMatrix:
    """Test transition matrix (lamb) for Multicell."""
    
    def test_default_lamb_shape(
        self, simple_grn, simple_receptor_grn, simple_cell_communication
    ):
        """Test that default lamb matrix has correct dimensions."""
        ct_a = Celltype(
            celltype_name="CellA",
            grn_graph=simple_grn.copy(),
            receptor_grn_bipartite=simple_receptor_grn.copy()
        )
        
        ct_b = Celltype(
            celltype_name="CellB",
            grn_graph=simple_grn.copy(),
            receptor_grn_bipartite=simple_receptor_grn.copy()
        )
        
        mc = Multicell(
            celltypes=[ct_a, ct_b],
            cell_communication_graph=simple_cell_communication
        )
        
        # Should have layers for: cell_comm, CellA_receptor, CellA_grn, CellB_receptor, CellB_grn
        expected_layers = 5
        assert mc.lamb.shape == (expected_layers, expected_layers)
    
    def test_lamb_allows_grn_transitions(
        self, simple_grn, simple_receptor_grn, simple_cell_communication
    ):
        """Test that lamb allows transitions within GRN layers."""
        ct_a = Celltype(
            celltype_name="CellA",
            grn_graph=simple_grn.copy(),
            receptor_grn_bipartite=simple_receptor_grn.copy()
        )
        
        mc = Multicell(
            celltypes=[ct_a],
            cell_communication_graph=simple_cell_communication
        )
        
        # GRN should allow self-transitions (normalized, so > 0)
        assert mc.lamb.loc["CellA_grn", "CellA_grn"] > 0


class TestMulticellRenaming:
    """Test celltype renaming in Multicell context."""
    
    def test_rename_via_dict(
        self, simple_grn, simple_receptor_grn, simple_cell_communication
    ):
        """Test renaming celltypes using dictionary keys."""
        ct_a = Celltype(
            celltype_name="OldA",
            grn_graph=simple_grn.copy(),
            receptor_grn_bipartite=simple_receptor_grn.copy()
        )
        
        ct_b = Celltype(
            celltype_name="OldB",
            grn_graph=simple_grn.copy(),
            receptor_grn_bipartite=simple_receptor_grn.copy()
        )
        
        mc = Multicell(
            celltypes={"NewA": ct_a, "NewB": ct_b},
            cell_communication_graph=simple_cell_communication
        )
        
        assert "NewA" in mc.celltypes_names
        assert "NewB" in mc.celltypes_names
        assert "OldA" not in mc.celltypes_names
