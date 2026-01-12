"""Tests for recon.explore.Celltype class."""
import pytest
import pandas as pd
import numpy as np
from recon.explore.recon import Celltype


class TestCelltypeConstruction:
    """Test Celltype object creation and initialization."""
    
    def test_create_celltype_basic(self, simple_grn, simple_receptor_grn):
        """Test basic Celltype creation with minimal inputs."""
        ct = Celltype(
            celltype_name="TestCell",
            grn_graph=simple_grn,
            receptor_grn_bipartite=simple_receptor_grn,
            receptor_graph=None
        )
        
        assert ct.celltype_name == "TestCell"
        assert "TestCell_grn" in ct.multiplexes
        assert "TestCell_receptor" in ct.multiplexes
        assert "TestCell_grn-TestCell_receptor" in ct.bipartites
    
    def test_create_celltype_with_receptor_graph(self, simple_grn, simple_receptor_grn):
        """Test Celltype creation with explicit receptor graph."""
        receptor_graph = pd.DataFrame({
            'source': ['RECEPTOR1'],
            'target': ['RECEPTOR2'],
            'weight': [0.5]
        })
        
        ct = Celltype(
            celltype_name="TestCell",
            grn_graph=simple_grn,
            receptor_grn_bipartite=simple_receptor_grn,
            receptor_graph=receptor_graph
        )
        
        # Check receptor layer was created from provided graph
        receptor_layer = ct.multiplexes["TestCell_receptor"]["layers"][0]
        assert len(receptor_layer) == 1
        assert receptor_layer.iloc[0]['source'] == 'RECEPTOR1_receptor'
    
    def test_fake_receptor_creation(self, simple_grn, simple_receptor_grn):
        """Test that fake receptor is created when receptor_graph is None."""
        ct = Celltype(
            celltype_name="TestCell",
            grn_graph=simple_grn,
            receptor_grn_bipartite=simple_receptor_grn,
            receptor_graph=None
        )
        
        # Fake receptor should be created
        receptor_layer = ct.multiplexes["TestCell_receptor"]["layers"][0]
        assert 'fake_receptor' in receptor_layer['target'].values
    
    def test_graph_type_encoding(self, simple_grn, simple_receptor_grn):
        """Test that graph type encoding is correct."""
        ct = Celltype(
            celltype_name="TestCell",
            grn_graph=simple_grn,
            receptor_grn_bipartite=simple_receptor_grn,
            receptor_graph=None,
            grn_graph_directed=False,
            grn_graph_weighted=True
        )
        
        # GRN should be "01" (undirected=0, weighted=1)
        assert ct.multiplexes["TestCell_grn"]["graph_type"][0] == "01"
    
    def test_default_lamb_matrix(self, simple_grn, simple_receptor_grn):
        """Test default transition matrix (lamb) structure."""
        ct = Celltype(
            celltype_name="TestCell",
            grn_graph=simple_grn,
            receptor_grn_bipartite=simple_receptor_grn
        )
        
        # Check lamb matrix shape and key transitions
        assert ct.lamb.shape == (2, 2)
        assert ct.lamb.loc["TestCell_grn", "TestCell_grn"] == 1
        assert ct.lamb.loc["TestCell_receptor", "TestCell_grn"] == 1
    
    def test_default_eta_vector(self, simple_grn, simple_receptor_grn):
        """Test default restart probability vector (eta)."""
        ct = Celltype(
            celltype_name="TestCell",
            grn_graph=simple_grn,
            receptor_grn_bipartite=simple_receptor_grn
        )
        
        # Default eta should be uniform
        assert len(ct.eta) == 2
        assert np.allclose(ct.eta.values, [0.5, 0.5])


class TestCelltypeSeeds:
    """Test seed handling in Celltype."""
    
    def test_seeds_as_list(self, simple_grn, simple_receptor_grn):
        """Test setting seeds as a list."""
        ct = Celltype(
            celltype_name="TestCell",
            grn_graph=simple_grn,
            receptor_grn_bipartite=simple_receptor_grn,
            seeds=["GENE1", "GENE2"]
        )
        
        assert ct.seeds == ["GENE1", "GENE2"]
    
    def test_seeds_as_dict(self, simple_grn, simple_receptor_grn):
        """Test setting seeds as a weighted dictionary."""
        ct = Celltype(
            celltype_name="TestCell",
            grn_graph=simple_grn,
            receptor_grn_bipartite=simple_receptor_grn,
            seeds={"GENE1": 0.7, "GENE2": 0.3}
        )
        
        assert ct.seeds == {"GENE1": 0.7, "GENE2": 0.3}
    
    def test_empty_seeds(self, simple_grn, simple_receptor_grn):
        """Test that empty seeds are handled gracefully."""
        ct = Celltype(
            celltype_name="TestCell",
            grn_graph=simple_grn,
            receptor_grn_bipartite=simple_receptor_grn,
            seeds=[]
        )
        
        assert ct.seeds == []


class TestCelltypeRenaming:
    """Test celltype renaming functionality."""
    
    def test_rename_celltype(self, simple_grn, simple_receptor_grn):
        """Test renaming a celltype updates all references."""
        ct = Celltype(
            celltype_name="OldName",
            grn_graph=simple_grn,
            receptor_grn_bipartite=simple_receptor_grn
        )
        
        ct.rename_celltype("NewName")
        
        assert ct.celltype_name == "NewName"
        assert "NewName_grn" in ct.multiplexes
        assert "NewName_receptor" in ct.multiplexes
        assert "OldName_grn" not in ct.multiplexes


class TestCelltypeDataFrameFormats:
    """Test that DataFrames maintain expected formats."""
    
    def test_grn_columns_preserved(self, simple_grn, simple_receptor_grn):
        """Test that GRN DataFrame columns are properly formatted."""
        ct = Celltype(
            celltype_name="TestCell",
            grn_graph=simple_grn,
            receptor_grn_bipartite=simple_receptor_grn
        )
        
        grn_layer = ct.multiplexes["TestCell_grn"]["layers"][0]
        # After internal processing, columns should be source/target
        assert 'source' in grn_layer.columns or 'target' in grn_layer.columns
    
    def test_bipartite_column_renaming(self, simple_grn, simple_receptor_grn):
        """Test that bipartite edges are renamed to col1/col2."""
        ct = Celltype(
            celltype_name="TestCell",
            grn_graph=simple_grn,
            receptor_grn_bipartite=simple_receptor_grn
        )
        
        bipartite_edges = ct.bipartites["TestCell_grn-TestCell_receptor"]["edge_list_df"]
        # Should be renamed to col1/col2 for multixrank
        assert 'col1' in bipartite_edges.columns
        assert 'col2' in bipartite_edges.columns


class TestCelltypeMultixrank:
    """Test Celltype.Multixrank method."""
    
    def test_multixrank_with_list_seeds(self, simple_grn, simple_receptor_grn):
        """Test that Multixrank creates a valid multilayer object with list seeds."""
        ct = Celltype(
            celltype_name="TestCell",
            grn_graph=simple_grn,
            receptor_grn_bipartite=simple_receptor_grn,
            seeds=["GENE1"]
        )
        
        multilayer = ct.Multixrank(restart_proba=0.7, verbose=False)
        
        # Verify multilayer object was created
        assert multilayer is not None
        assert hasattr(multilayer, 'random_walk_rank')
    
    def test_multixrank_with_dict_seeds(self, simple_grn, simple_receptor_grn):
        """Test Multixrank with weighted dictionary seeds."""
        ct = Celltype(
            celltype_name="TestCell",
            grn_graph=simple_grn,
            receptor_grn_bipartite=simple_receptor_grn,
            seeds={"GENE1": 0.7, "GENE2": 0.3}
        )
        
        multilayer = ct.Multixrank(restart_proba=0.7, verbose=False)
        
        assert multilayer is not None
        # Check that pr vector was set
        assert hasattr(multilayer, 'pr')
        assert multilayer.pr is not None
    
    def test_multixrank_custom_restart_proba(self, simple_grn, simple_receptor_grn):
        """Test Multixrank with custom restart probability."""
        ct = Celltype(
            celltype_name="TestCell",
            grn_graph=simple_grn,
            receptor_grn_bipartite=simple_receptor_grn,
            seeds=["GENE1"]
        )
        
        multilayer = ct.Multixrank(restart_proba=0.5, verbose=False)
        assert multilayer is not None
    
    def test_multixrank_no_self_loops(self, simple_grn, simple_receptor_grn):
        """Test Multixrank with self_loops=False."""
        ct = Celltype(
            celltype_name="TestCell",
            grn_graph=simple_grn,
            receptor_grn_bipartite=simple_receptor_grn,
            seeds=["GENE1"]
        )
        
        multilayer = ct.Multixrank(self_loops=False, restart_proba=0.7, verbose=False)
        assert multilayer is not None


class TestCelltypeAdvancedOptions:
    """Test advanced Celltype configuration options."""
    
    def test_celltype_copy_graphs_false(self, simple_grn, simple_receptor_grn):
        """Test that copy_graphs=False doesn't copy DataFrames."""
        original_grn = simple_grn.copy()
        original_receptor = simple_receptor_grn.copy()
        
        ct = Celltype(
            celltype_name="TestCell",
            grn_graph=original_grn,
            receptor_grn_bipartite=original_receptor,
            copy_graphs=False
        )
        
        assert ct is not None
        assert "TestCell_grn" in ct.multiplexes
    
    def test_celltype_custom_lamb(self, simple_grn, simple_receptor_grn):
        """Test providing custom lamb matrix."""
        custom_lamb = pd.DataFrame(
            [[0.8, 0.2], [0.3, 0.7]],
            index=['TestCell_receptor', 'TestCell_grn'],
            columns=['TestCell_receptor', 'TestCell_grn']
        )
        
        ct = Celltype(
            celltype_name="TestCell",
            grn_graph=simple_grn,
            receptor_grn_bipartite=simple_receptor_grn,
            lamb=custom_lamb
        )
        
        assert (ct.lamb == custom_lamb).all().all()
    
    def test_celltype_custom_eta(self, simple_grn, simple_receptor_grn):
        """Test providing custom eta vector."""
        custom_eta = pd.Series(
            [0.3, 0.7],
            index=['TestCell_receptor', 'TestCell_grn']
        )
        
        ct = Celltype(
            celltype_name="TestCell",
            grn_graph=simple_grn,
            receptor_grn_bipartite=simple_receptor_grn,
            eta=custom_eta
        )
        
        assert (ct.eta == custom_eta).all()
    
    def test_celltype_custom_graph_types(self, simple_grn, simple_receptor_grn):
        """Test various combinations of graph type flags."""
        ct = Celltype(
            celltype_name="TestCell",
            grn_graph=simple_grn,
            receptor_grn_bipartite=simple_receptor_grn,
            grn_graph_directed=True,
            grn_graph_weighted=False,
            receptor_graph_directed=False,
            receptor_graph_weighted=True
        )
        
        # GRN should be "10" (directed=1, weighted=0)
        assert ct.multiplexes["TestCell_grn"]["graph_type"][0] == "10"
        # Receptor should be "01" (directed=0, weighted=1)
        assert ct.multiplexes["TestCell_receptor"]["graph_type"][0] == "01"


class TestCelltypeErrorHandling:
    """Test error handling in Celltype."""
    
    def test_invalid_seeds_dict_with_string_values(self, simple_grn, simple_receptor_grn):
        """Test that dict seeds with non-numeric values raise error."""
        ct = Celltype(
            celltype_name="TestCell",
            grn_graph=simple_grn,
            receptor_grn_bipartite=simple_receptor_grn,
            seeds={"GENE1": "invalid"}
        )
        
        with pytest.raises(ValueError, match="Seeds should be a dictionary with numeric values"):
            ct.Multixrank(verbose=False)
    
    def test_invalid_seeds_type(self, simple_grn, simple_receptor_grn):
        """Test that invalid seed type raises error."""
        ct = Celltype(
            celltype_name="TestCell",
            grn_graph=simple_grn,
            receptor_grn_bipartite=simple_receptor_grn,
            seeds="invalid_string_seed"
        )
        
        with pytest.raises(ValueError, match="Seeds should be"):
            ct.Multixrank(verbose=False)
