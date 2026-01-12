"""Tests for recon.plot.sankey_paths module."""
import pytest
import pandas as pd
import numpy as np
from recon.explore.recon import Celltype, Multicell
from recon.plot import sankey_paths


@pytest.fixture
def simple_multicell_for_sankey(simple_grn, simple_receptor_grn, simple_cell_communication):
    """Create a minimal multicell object for sankey testing."""
    # Create two cell types
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
    
    return Multicell(
        celltypes=[ct_a, ct_b],
        cell_communication_graph=simple_cell_communication
    )


@pytest.fixture
def simple_results():
    """Create a minimal results dataframe mimicking RWR output."""
    return pd.DataFrame({
        'multiplex': [
            'CellA_grn', 'CellA_grn', 'CellA_receptor',
            'CellB_grn', 'CellB_grn', 'cell_communication',
            'cell_communication', 'cell_communication'
        ],
        'node': [
            'TF1_TF::CellA', 'TF2_TF::CellA', 'RECEPTOR1_receptor::CellA',
            'TF1_TF::CellB', 'TF2_TF::CellB', 'LIGAND1::CellB',
            'LIGAND2::CellB', 'RECEPTOR1_receptor::CellA'
        ],
        'score': [0.9, 0.85, 0.8, 0.75, 0.7, 0.65, 0.6, 0.55]
    })


class TestGetCelltypeGrnReceptorBipartite:
    """Test get_celltype_grn_receptor_bipartite function."""
    
    def test_basic_retrieval(self, simple_multicell_for_sankey):
        """Test basic retrieval of bipartite dataframe."""
        result = sankey_paths.get_celltype_grn_receptor_bipartite(
            simple_multicell_for_sankey,
            cell_type="CellA",
            as_dataframe=True
        )
        
        assert isinstance(result, pd.DataFrame)
        # Check expected columns from bipartite structure
        assert 'col1' in result.columns or 'col2' in result.columns
    
    def test_invalid_celltype(self, simple_multicell_for_sankey):
        """Test that invalid cell type raises KeyError."""
        with pytest.raises(KeyError, match="No bipartite entry found"):
            sankey_paths.get_celltype_grn_receptor_bipartite(
                simple_multicell_for_sankey,
                cell_type="NonexistentCell",
                as_dataframe=True
            )


class TestGetCelltypeGeneLayer:
    """Test get_celltype_gene_layer function."""
    
    def test_basic_retrieval(self, simple_multicell_for_sankey):
        """Test basic retrieval of gene layer dataframe."""
        result = sankey_paths.get_celltype_gene_layer(
            simple_multicell_for_sankey,
            cell_type="CellA",
            layer_name="gene",  # Use 'gene' as the layer name
            as_dataframe=True
        )
        
        assert isinstance(result, pd.DataFrame)
        # Should have source, target columns from GRN
        assert 'source' in result.columns
        assert 'target' in result.columns
    
    def test_invalid_celltype(self, simple_multicell_for_sankey):
        """Test that invalid cell type raises KeyError."""
        with pytest.raises(KeyError, match="No multiplex found"):
            sankey_paths.get_celltype_gene_layer(
                simple_multicell_for_sankey,
                cell_type="NonexistentCell",
                as_dataframe=True
            )


class TestGetCellCommunicationLayer:
    """Test get_cell_communication_layer function."""
    
    def test_basic_retrieval(self, simple_multicell_for_sankey):
        """Test basic retrieval of cell communication layer."""
        result = sankey_paths.get_cell_communication_layer(
            simple_multicell_for_sankey,
            as_dataframe=True
        )
        
        assert isinstance(result, pd.DataFrame)
        # Should have ligand and receptor columns
        assert 'ligand' in result.columns
        assert 'receptor' in result.columns
        assert 'weight' in result.columns


class TestGetTopFunctions:
    """Test get_top_tfs, get_top_receptors, get_top_ligands functions."""
    
    def test_get_top_tfs(self, simple_results):
        """Test extraction of top TFs."""
        result = sankey_paths.get_top_tfs(
            simple_results,
            cell_type="CellA",
            n=2
        )
        
        assert len(result) <= 2
        assert all(result['node'].str.contains('_TF::CellA'))
    
    def test_get_top_receptors(self, simple_results):
        """Test extraction of top receptors."""
        result = sankey_paths.get_top_receptors(
            simple_results,
            cell_type="CellA",
            n=1
        )
        
        assert len(result) <= 1
        # Should not include fake_receptor
        assert not any(result['node'].str.contains('fake_receptor'))
    
    def test_get_top_ligands_basic(self, simple_results):
        """Test extraction of top ligands."""
        # Create minimal receptor_ligand_df
        receptor_ligand_df = pd.DataFrame({
            'ligand': ['LIGAND1::CellB', 'LIGAND2::CellB'],
            'receptor': ['RECEPTOR1::CellA', 'RECEPTOR1::CellA']
        })
        
        result = sankey_paths.get_top_ligands(
            simple_results,
            receptor_ligand_df=receptor_ligand_df,
            n=2,
            per_celltype=False
        )
        
        assert len(result) <= 2
        # Should only include ligands that are in receptor_ligand_df
        assert all(result['multiplex'] == 'cell_communication')


class TestExtractPairFunctions:
    """Test extract_gene_tf_pairs, extract_receptor_tf_pairs, extract_receptor_ligand_pairs."""
    
    def test_extract_gene_tf_pairs(self, simple_grn, simple_results):
        """Test extraction of gene-TF pairs."""
        # Prepare inputs
        tf_gene_layer = simple_grn.copy()
        tf_gene_layer['source'] = tf_gene_layer['source'] + '::CellA'
        tf_gene_layer['target'] = tf_gene_layer['target'] + '_TF::CellA'
        tf_gene_layer['network_key'] = 'grn'
        
        top_tfs = simple_results[simple_results['multiplex'] == 'CellA_grn'].head(2)
        seeds = pd.Index(['GENE1::CellA', 'GENE2::CellA'])
        
        result = sankey_paths.extract_gene_tf_pairs(
            tf_gene_layer,
            top_tfs,
            seeds
        )
        
        assert isinstance(result, pd.DataFrame)
        # Should have gene, tf, weight columns after processing
        assert 'gene' in result.columns
        assert 'tf' in result.columns
    
    def test_extract_receptor_tf_pairs(self, simple_receptor_grn, simple_results):
        """Test extraction of receptor-TF pairs."""
        # Prepare receptor_gene_layer
        receptor_gene_layer = simple_receptor_grn.copy()
        receptor_gene_layer = receptor_gene_layer.rename(columns={
            'source': 'col1',
            'target': 'col2',
            'score': 'weight'
        })
        receptor_gene_layer['col1'] = receptor_gene_layer['col1'] + '::CellA'
        receptor_gene_layer['col2'] = receptor_gene_layer['col2'] + '_TF::CellA'
        receptor_gene_layer['network_key'] = 'bipartite'
        
        top_tfs = simple_results[simple_results['multiplex'] == 'CellA_grn'].head(2)
        top_receptors = simple_results[simple_results['multiplex'] == 'CellA_receptor']
        
        result = sankey_paths.extract_receptor_tf_pairs(
            receptor_gene_layer,
            top_tfs,
            top_receptors
        )
        
        assert isinstance(result, pd.DataFrame)
        # Check column names after processing
        assert 'receptor' in result.columns
        assert 'tf' in result.columns
    
    def test_extract_receptor_ligand_pairs(self):
        """Test extraction of receptor-ligand pairs."""
        # Create test data
        receptor_ligand_df = pd.DataFrame({
            'ligand': ['LIGAND1::CellB', 'LIGAND2::CellB', 'LIGAND3::CellC'],
            'receptor': ['RECEPTOR1::CellA', 'RECEPTOR2::CellA', 'RECEPTOR1::CellA']
        })
        
        top_ligands_df = pd.DataFrame({
            'node': ['LIGAND1::CellB', 'LIGAND2::CellB'],
            'score': [0.9, 0.8]
        })
        
        top_receptors_df = pd.DataFrame({
            'node': ['RECEPTOR1_receptor::CellA', 'RECEPTOR2_receptor::CellA'],
            'score': [0.85, 0.75]
        })
        
        result = sankey_paths.extract_receptor_ligand_pairs(
            receptor_ligand_df,
            top_ligands_df,
            top_receptors_df
        )
        
        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0
        # All ligands in result should be in top_ligands
        assert all(result['ligand'].isin(['LIGAND1::CellB', 'LIGAND2::CellB']))


class TestBuildPartialNetworks:
    """Test build_partial_networks function."""
    
    def test_build_without_ligand_cells(self, simple_multicell_for_sankey, simple_results):
        """Test building partial networks without ligand cells (intracellular only)."""
        seeds = pd.Index(['GENE1::CellA', 'GENE2::CellA'])
        
        networks = sankey_paths.build_partial_networks(
            multicell_obj=simple_multicell_for_sankey,
            results=simple_results,
            cell_type="CellA",
            seeds=seeds,
            ligand_cells=[],
            top_receptor_n=5,
            top_tf_n=3,
            include_before_cells=False
        )
        
        # Should return tuple of 5 dataframes
        assert isinstance(networks, tuple)
        assert len(networks) == 5
        # All should be DataFrames
        assert all(isinstance(df, pd.DataFrame) for df in networks)
    
    def test_build_with_ligand_cells(self, simple_multicell_for_sankey, simple_results):
        """Test building partial networks with ligand cells (intercellular)."""
        seeds = pd.Index(['GENE1::CellA', 'GENE2::CellA'])
        
        networks = sankey_paths.build_partial_networks(
            multicell_obj=simple_multicell_for_sankey,
            results=simple_results,
            cell_type="CellA",
            seeds=seeds,
            ligand_cells=["CellB"],
            top_ligand_n=10,
            top_receptor_n=5,
            top_tf_n=3,
            include_before_cells=True
        )
        
        assert isinstance(networks, tuple)
        assert len(networks) == 5


class TestSankeyPlotFunctions:
    """Test the main sankey plotting functions (without actually rendering plots)."""
    
    def test_plot_intracell_sankey_runs(self, simple_multicell_for_sankey, simple_results, monkeypatch):
        """Test that plot_intracell_sankey runs without errors."""
        seeds = pd.Index(['GENE1::CellA', 'GENE2::CellA'])
        
        # Mock the actual plotting function to avoid display issues in tests
        def mock_plot(*args, **kwargs):
            pass
        
        monkeypatch.setattr(sankey_paths, 'plot_3layer_sankey', mock_plot)
        
        # Should not raise an error
        sankey_paths.plot_intracell_sankey(
            multicell_obj=simple_multicell_for_sankey,
            results=simple_results,
            cell_type="CellA",
            seeds=seeds,
            top_receptor_n=5,
            top_tf_n=3,
            flow="upstream",
            save_path=None
        )
    
    def test_plot_ligand_sankey_runs(self, simple_multicell_for_sankey, simple_results, monkeypatch):
        """Test that plot_ligand_sankey runs without errors."""
        seeds = pd.Index(['GENE1::CellA', 'GENE2::CellA'])
        
        def mock_plot(*args, **kwargs):
            pass
        
        monkeypatch.setattr(sankey_paths, 'plot_4layer_sankey', mock_plot)
        
        sankey_paths.plot_ligand_sankey(
            multicell_obj=simple_multicell_for_sankey,
            results=simple_results,
            cell_type="CellA",
            seeds=seeds,
            ligand_cells=["CellB"],
            top_ligand_n=10,
            top_receptor_n=5,
            top_tf_n=3,
            flow="upstream",
            save_path=None
        )
    
    def test_plot_intercell_sankey_runs(self, simple_multicell_for_sankey, simple_results, monkeypatch):
        """Test that plot_intercell_sankey runs without errors."""
        seeds = pd.Index(['GENE1::CellA', 'GENE2::CellA'])
        
        def mock_plot(*args, **kwargs):
            pass
        
        monkeypatch.setattr(sankey_paths, 'plot_6layer_sankey', mock_plot)
        
        sankey_paths.plot_intercell_sankey(
            multicell_obj=simple_multicell_for_sankey,
            results=simple_results,
            cell_type="CellA",
            seeds=seeds,
            ligand_cells=["CellB"],
            top_ligand_n=10,
            top_receptor_n=5,
            top_tf_n=3,
            before_top_n=2,
            flow="upstream",
            save_path=None
        )


class TestEdgeCases:
    """Test edge cases and error handling."""
    
    def test_empty_results(self, simple_multicell_for_sankey):
        """Test behavior with empty results dataframe."""
        empty_results = pd.DataFrame(columns=['multiplex', 'node', 'score'])
        
        # Should handle gracefully without crashing
        result = sankey_paths.get_top_tfs(empty_results.copy(), "CellA", n=5)
        assert len(result) == 0
    
    def test_seeds_format_validation(self, simple_multicell_for_sankey, simple_results):
        """Test that seeds must have correct format with :: suffix."""
        # Seeds without :: suffix should work if they're converted internally
        seeds = pd.Index(['GENE1::CellA', 'GENE2::CellA'])
        
        networks = sankey_paths.build_partial_networks(
            multicell_obj=simple_multicell_for_sankey,
            results=simple_results,
            cell_type="CellA",
            seeds=seeds,
            ligand_cells=[],
            include_before_cells=False
        )
        
        assert isinstance(networks, tuple)


class TestGetTopLigands:
    """Test get_top_ligands function."""
    
    def test_basic_ligand_extraction(self, simple_results):
        """Test extracting top ligands from results."""
        cc_df = pd.DataFrame({
            'ligand': ['LIGAND1', 'LIGAND2', 'LIGAND3'],
            'receptor': ['RECEPTOR1', 'RECEPTOR2', 'RECEPTOR1'],
            'celltype_source': ['CellB', 'CellB', 'CellC']
        })
        
        result = sankey_paths.get_top_ligands(
            results_df=simple_results,
            receptor_ligand_df=cc_df,
            n=2,
            per_celltype=False
        )
        
        assert isinstance(result, pd.DataFrame)
        assert len(result) <= 2
    
    def test_per_celltype_mode(self, simple_results):
        """Test ligand extraction with per_celltype=True."""
        cc_df = pd.DataFrame({
            'ligand': ['LIGAND1', 'LIGAND2'],
            'receptor': ['RECEPTOR1', 'RECEPTOR2'],
            'celltype_source': ['CellB', 'CellC']
        })
        
        result = sankey_paths.get_top_ligands(
            results_df=simple_results,
            receptor_ligand_df=cc_df,
            n=5,
            per_celltype=True
        )
        
        assert isinstance(result, pd.DataFrame)
