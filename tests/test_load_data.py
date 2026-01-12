"""
Tests for recon.data.load_data module.

Tests the receptor gene loading functionality.
"""

import pytest
import pandas as pd
from recon.data.load_data import load_receptor_genes, receptor_gene_resources, fetch_tutorial_data
import os


class TestReceptorGeneResources:
    """Test the receptor_gene_resources constant."""
    
    def test_resources_list_exists(self):
        """Test that receptor_gene_resources is defined."""
        assert isinstance(receptor_gene_resources, list)
        assert len(receptor_gene_resources) > 0
    
    def test_resources_list_content(self):
        """Test that expected resources are in the list."""
        assert "human_receptor_gene_from_NichenetPKN" in receptor_gene_resources
        assert "mouse_receptor_gene_from_NichenetPKN" in receptor_gene_resources


class TestLoadReceptorGenes:
    """Test load_receptor_genes function."""
    
    def test_load_human_receptor_genes(self):
        """Test loading human receptor genes."""
        df = load_receptor_genes("human_receptor_gene_from_NichenetPKN")
        
        # Should return a DataFrame
        assert isinstance(df, pd.DataFrame)
        
        # Should have required columns
        assert "source" in df.columns
        assert "target" in df.columns
        
        # Should not be empty
        assert len(df) > 0
        
        # Source should be receptors (genes), target should be downstream genes
        assert all(isinstance(val, str) for val in df["source"].head())
        assert all(isinstance(val, str) for val in df["target"].head())
    
    def test_load_mouse_receptor_genes(self):
        """Test loading mouse receptor genes."""
        df = load_receptor_genes("mouse_receptor_gene_from_NichenetPKN")
        
        # Should return a DataFrame
        assert isinstance(df, pd.DataFrame)
        
        # Should have required columns
        assert "source" in df.columns
        assert "target" in df.columns
        
        # Should not be empty
        assert len(df) > 0
    
    def test_human_vs_mouse_different(self):
        """Test that human and mouse data are different."""
        human_df = load_receptor_genes("human_receptor_gene_from_NichenetPKN")
        mouse_df = load_receptor_genes("mouse_receptor_gene_from_NichenetPKN")
        
        # Should have different data (different gene names)
        # Not all genes should be identical
        human_genes = set(human_df["source"].unique())
        mouse_genes = set(mouse_df["source"].unique())
        
        # There should be some differences
        assert human_genes != mouse_genes
    
    def test_invalid_resource_name(self):
        """Test that invalid resource names raise ValueError."""
        with pytest.raises(ValueError, match="must be in"):
            load_receptor_genes("invalid_receptor_list")
    
    def test_typo_in_resource_name(self):
        """Test that typos in resource names raise ValueError."""
        with pytest.raises(ValueError):
            load_receptor_genes("human_receptor_gene_from_Nichenet")  # Missing PKN
    
    def test_none_as_input(self):
        """Test that None as input raises appropriate error."""
        with pytest.raises((ValueError, TypeError)):
            load_receptor_genes(None)
    
    def test_empty_string_as_input(self):
        """Test that empty string raises ValueError."""
        with pytest.raises(ValueError):
            load_receptor_genes("")
    
    def test_dataframe_structure(self):
        """Test the structure of returned DataFrames."""
        df = load_receptor_genes("human_receptor_gene_from_NichenetPKN")
        
        # Should have at least source and target columns
        assert len(df.columns) >= 2
        
        # No null values in critical columns
        assert df["source"].notna().all()
        assert df["target"].notna().all()
    
    def test_load_multiple_times_same_result(self):
        """Test that loading the same resource multiple times gives consistent results."""
        df1 = load_receptor_genes("human_receptor_gene_from_NichenetPKN")
        df2 = load_receptor_genes("human_receptor_gene_from_NichenetPKN")
        
        # Should be identical - compare shapes and values
        assert df1.shape == df2.shape
        assert list(df1.columns) == list(df2.columns)
        assert df1.equals(df2)
    
    def test_receptor_genes_are_unique_per_target(self):
        """Test that receptor-gene pairs are properly structured."""
        df = load_receptor_genes("mouse_receptor_gene_from_NichenetPKN")
        
        # Each receptor can have multiple target genes
        # Reset index to avoid RangeIndex assertion issues in pandas groupby
        receptor_counts = df.reset_index(drop=True).groupby("source").size()
        
        # Should have multiple targets per receptor on average
        assert receptor_counts.mean() >= 1
        
        # But total pairs should be reasonable
        assert len(df) > len(df["source"].unique())


class TestDataIntegration:
    """Test integration with ReCoN's data structures."""
    
    def test_compatible_with_celltype_construction(self):
        """Test that loaded data is compatible with Celltype construction."""
        df = load_receptor_genes("mouse_receptor_gene_from_NichenetPKN")
        
        # Should be usable as receptor_grn_bipartite
        # Check it has the minimal required columns
        required_cols = ["source", "target"]
        for col in required_cols:
            assert col in df.columns
        
        # Data types should be strings (gene names)
        assert df["source"].dtype == object
        assert df["target"].dtype == object
    
    def test_receptor_names_format(self):
        """Test that receptor names don't have special suffixes yet."""
        df = load_receptor_genes("human_receptor_gene_from_NichenetPKN")
        
        # Source should be plain gene names (no _receptor suffix yet)
        # That suffix is added by Celltype constructor
        sample_receptors = df["source"].head(10)
        
        # Most should not end with _receptor (that's added later)
        non_suffixed = sum(1 for r in sample_receptors if not r.endswith("_receptor"))
        assert non_suffixed >= 5  # Most should be plain names


class TestFetchTutorialData:
    """Test fetch_tutorial_data function with the smallest file."""

    def test_fetch_smallest_file(self, tmp_path):
        """Test downloading the smallest tutorial file."""
        # Define the smallest file to test
        smallest_file = "perturbation_tuto/rna.h5ad"
        
        # Temporary directory for testing
        test_data_dir = tmp_path / "data"
        test_data_dir.mkdir()

        # Fetch the file
        file_path = fetch_tutorial_data(smallest_file, data_dir=str(test_data_dir))

        # Check if the file exists
        assert os.path.exists(file_path)

        # Check if the file is not empty
        assert os.path.getsize(file_path) > 0
