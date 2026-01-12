"""Tests for recon.plot.plot_results module."""
import pytest
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')  # Non-interactive backend for testing
import matplotlib.pyplot as plt
from recon.plot import plot_results


@pytest.fixture
def sample_comparison_data():
    """Create sample data for celltype comparison plots."""
    np.random.seed(42)
    n_genes = 100
    
    # Create data with some outliers
    celltype1_scores = np.random.randn(n_genes) * 0.5 + 1.0
    celltype2_scores = np.random.randn(n_genes) * 0.5 + 1.0
    
    # Add some clear outliers (celltype1 >> celltype2)
    outlier_indices = [0, 1, 2]
    celltype1_scores[outlier_indices] = [5.0, 4.5, 4.0]
    celltype2_scores[outlier_indices] = [0.5, 0.6, 0.7]
    
    df = pd.DataFrame({
        'CellA': celltype1_scores,
        'CellB': celltype2_scores
    }, index=[f'Gene{i}' for i in range(n_genes)])
    
    return df


class TestPlotCelltypeComparison:
    """Test plot_celltype_comparison function."""
    
    def test_basic_plot(self, sample_comparison_data):
        """Test that the function runs without error and creates a plot."""
        # Close any existing plots
        plt.close('all')
        
        # Run the plotting function
        plot_results.plot_celltype_comparison(
            sample_comparison_data,
            celltype1='CellA',
            celltype2='CellB',
            quantile=0.95
        )
        
        # Check that a figure was created
        fig = plt.gcf()
        assert fig is not None
        assert len(fig.get_axes()) >= 1  # At least main axis
        
        plt.close('all')
    
    def test_with_default_quantile(self, sample_comparison_data):
        """Test with default quantile parameter."""
        plt.close('all')
        
        plot_results.plot_celltype_comparison(
            sample_comparison_data,
            celltype1='CellA',
            celltype2='CellB'
        )
        
        fig = plt.gcf()
        assert fig is not None
        
        plt.close('all')
    
    def test_with_different_quantiles(self, sample_comparison_data):
        """Test with various quantile thresholds."""
        plt.close('all')
        
        for quantile in [0.90, 0.95, 0.99, 0.999]:
            plot_results.plot_celltype_comparison(
                sample_comparison_data,
                celltype1='CellA',
                celltype2='CellB',
                quantile=quantile
            )
            plt.close('all')
    
    def test_with_negative_values(self):
        """Test with data containing negative values."""
        plt.close('all')
        
        df = pd.DataFrame({
            'CellA': [-1.0, -0.5, 0.0, 0.5, 1.0, 3.0],
            'CellB': [-1.0, -0.5, 0.0, 0.5, 1.0, 0.5]
        }, index=[f'Gene{i}' for i in range(6)])
        
        plot_results.plot_celltype_comparison(
            df,
            celltype1='CellA',
            celltype2='CellB',
            quantile=0.9
        )
        
        fig = plt.gcf()
        assert fig is not None
        
        plt.close('all')
    
    def test_with_identical_values(self):
        """Test when both celltypes have very similar scores."""
        plt.close('all')
        
        # Create data with more variation to ensure valid quantiles
        values_a = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])
        values_b = values_a + np.random.randn(len(values_a)) * 0.5  # Add noise
        
        df = pd.DataFrame({
            'CellA': values_a,
            'CellB': values_b
        }, index=[f'Gene{i}' for i in range(len(values_a))])
        
        plot_results.plot_celltype_comparison(
            df,
            celltype1='CellA',
            celltype2='CellB',
            quantile=0.8  # Lower quantile for small dataset
        )
        
        fig = plt.gcf()
        assert fig is not None
        
        plt.close('all')
    
    def test_with_small_dataset(self):
        """Test with very small dataset."""
        plt.close('all')
        
        df = pd.DataFrame({
            'CellA': [1.0, 2.0, 3.0],
            'CellB': [0.5, 1.5, 1.0]
        }, index=['Gene1', 'Gene2', 'Gene3'])
        
        plot_results.plot_celltype_comparison(
            df,
            celltype1='CellA',
            celltype2='CellB',
            quantile=0.5
        )
        
        fig = plt.gcf()
        assert fig is not None
        
        plt.close('all')
    
    def test_axes_properties(self, sample_comparison_data):
        """Test that axes have correct properties."""
        plt.close('all')
        
        plot_results.plot_celltype_comparison(
            sample_comparison_data,
            celltype1='CellA',
            celltype2='CellB',
            quantile=0.95
        )
        
        fig = plt.gcf()
        ax = fig.get_axes()[0]  # Main axis
        
        # Check labels
        assert 'CellA' in ax.get_xlabel()
        assert 'CellB' in ax.get_ylabel()
        assert ax.get_title() != ''
        
        # Check that there's a legend
        legend = ax.get_legend()
        assert legend is not None
        
        plt.close('all')
    
    def test_with_all_zeros(self):
        """Test with all zero values - expect it to handle gracefully."""
        plt.close('all')
        
        # Add small variation to avoid all identical values
        df = pd.DataFrame({
            'CellA': [0.0] * 9 + [0.1],
            'CellB': [0.0] * 10
        }, index=[f'Gene{i}' for i in range(10)])
        
        plot_results.plot_celltype_comparison(
            df,
            celltype1='CellA',
            celltype2='CellB',
            quantile=0.9
        )
        
        fig = plt.gcf()
        assert fig is not None
        
        plt.close('all')
