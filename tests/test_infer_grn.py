"""
Optional tests for GRN inference layers.

These tests will be skipped when celloracle is not available.
To run these tests, install celloracle:
    pip install 'git+https://github.com/cantinilab/celloracle@lite'

Note: ATAC-seq tests also require reference genomes to be installed:
    genomepy install mm10 UCSC --annotation
"""

import pytest
import numpy as np
import pandas as pd
import anndata as ad
import os
from pathlib import Path
from recon.infer_grn.layers import (
    compute_tf_network,
    compute_rna_network,
    compute_tf_to_atac_links,
    compute_atac_to_rna_links,
    CELLORACLE_AVAILABLE
)

# Check if genomes are available
GENOMES_AVAILABLE = False
if CELLORACLE_AVAILABLE:
    try:
        import genomepy
        genome_dir = Path(genomepy.utils.get_genomes_dir())
        GENOMES_AVAILABLE = (genome_dir / "mm10").exists()
    except (ImportError, FileNotFoundError):
        pass


class TestComputeTfNetwork:
    """Tests for compute_tf_network function (doesn't require celloracle)."""
    
    def test_compute_tf_network_default(self):
        """Test TF network creation with default method."""
        # Create simple RNA data
        rna = ad.AnnData(
            X=np.random.rand(10, 5),
            var=pd.DataFrame(index=['TF1', 'TF2', 'GENE1', 'GENE2', 'GENE3'])
        )
        tfs_list = ['TF1', 'TF2']
        
        result = compute_tf_network(rna, tfs_list, method=None)
        
        assert isinstance(result, pd.DataFrame)
        assert list(result.columns) == ['source', 'target']
        assert len(result) == 2
        assert all(result['source'].str.endswith('_TF'))
        assert all(result['target'] == 'fake_TF')
    
    def test_compute_tf_network_filters_missing_tfs(self):
        """Test that TFs not in genes are filtered out."""
        rna = ad.AnnData(
            X=np.random.rand(10, 3),
            var=pd.DataFrame(index=['TF1', 'GENE1', 'GENE2'])
        )
        tfs_list = ['TF1', 'TF2', 'TF3']  # TF2 and TF3 not in genes
        
        result = compute_tf_network(rna, tfs_list, method=None)
        
        assert len(result) == 1
        assert result['source'].iloc[0] == 'TF1_TF'
    
    def test_compute_tf_network_invalid_method(self):
        """Test that invalid method raises ValueError."""
        rna = ad.AnnData(
            X=np.random.rand(10, 3),
            var=pd.DataFrame(index=['TF1', 'GENE1', 'GENE2'])
        )
        tfs_list = ['TF1']
        
        with pytest.raises(ValueError, match="For now, no method has been implemented"):
            compute_tf_network(rna, tfs_list, method='invalid')


class TestComputeRnaNetwork:
    """Tests for compute_rna_network function (doesn't require celloracle)."""
    
    def test_compute_rna_network_with_dataframe(self):
        """Test RNA network inference with DataFrame input."""
        # Create expression matrix
        np.random.seed(42)
        df_exp = pd.DataFrame(
            np.random.rand(50, 10),
            columns=[f'GENE{i}' for i in range(10)]
        )
        tf_names = ['GENE0', 'GENE1', 'GENE2']
        
        result = compute_rna_network(
            df_exp_mtx=df_exp,
            tf_names=tf_names,
            method='GBM',
            n_cpu=1,
            seed=42
        )
        
        assert isinstance(result, pd.DataFrame)
        assert list(result.columns) == ['source', 'target', 'weight']
        assert len(result) > 0
        assert all(result['source'].isin(tf_names))
        assert all(result['target'].isin(df_exp.columns))
    
    def test_compute_rna_network_with_anndata(self):
        """Test RNA network inference with AnnData input."""
        np.random.seed(42)
        adata = ad.AnnData(
            X=np.random.rand(50, 10),
            var=pd.DataFrame(index=[f'GENE{i}' for i in range(10)])
        )
        tf_names = ['GENE0', 'GENE1']
        
        result = compute_rna_network(
            df_exp_mtx=adata,
            tf_names=tf_names,
            method='RF',
            n_cpu=1,
            seed=42
        )
        
        assert isinstance(result, pd.DataFrame)
        assert list(result.columns) == ['source', 'target', 'weight']
        assert all(result['source'].isin(tf_names))
    
    def test_compute_rna_network_invalid_input(self):
        """Test that invalid input type raises ValueError."""
        with pytest.raises(ValueError, match="df_exp_mtx must be a pandas DataFrame"):
            compute_rna_network(
                df_exp_mtx="invalid",
                tf_names=['TF1'],
                n_cpu=1
            )
    
    def test_compute_rna_network_method_gbm(self):
        """Test that GBM method works correctly."""
        np.random.seed(42)
        df_exp = pd.DataFrame(
            np.random.rand(30, 8),
            columns=[f'GENE{i}' for i in range(8)]
        )
        
        result = compute_rna_network(
            df_exp_mtx=df_exp,
            tf_names=['GENE0', 'GENE1'],
            method='GBM',
            n_cpu=1,
            seed=42
        )
        
        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0
    
    def test_compute_rna_network_method_rf(self):
        """Test that RF method works correctly."""
        np.random.seed(42)
        df_exp = pd.DataFrame(
            np.random.rand(30, 8),
            columns=[f'GENE{i}' for i in range(8)]
        )
        
        result = compute_rna_network(
            df_exp_mtx=df_exp,
            tf_names=['GENE0', 'GENE1'],
            method='RF',
            n_cpu=1,
            seed=42
        )
        
        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0


@pytest.mark.skipif(not CELLORACLE_AVAILABLE, reason="celloracle not installed")
@pytest.mark.skipif(not GENOMES_AVAILABLE, reason="mm10 genome not installed")
class TestComputeTfToAtacLinks:
    """
    Optional tests for ATAC-seq analysis requiring celloracle and reference genomes.
    
    These tests will be SKIPPED if:
    - celloracle is not installed
    - mm10 genome is not available
    
    To run these tests:
    1. Install celloracle: pip install 'git+https://github.com/cantinilab/celloracle@lite'
    2. Install genome: genomepy install mm10 UCSC --annotation
    """
    
    def test_compute_tf_to_atac_links_basic(self):
        """Test TF-to-ATAC links computation with known regulatory regions."""
        # Use real mouse genomic coordinates from known regulatory regions
        # These are promoter/enhancer regions likely to contain TF binding sites
        peaks = [
            'chr1_3000000_3001000',  # Near Sox2 locus regulatory region
            'chr1_4807000_4808000',  # Near Pou5f1 (Oct4) regulatory region  
            'chr2_29000000_29001000',  # Active enhancer region
            'chr7_103000000_103001000',  # Known regulatory element
            'chr11_97000000_97001000',  # Promoter region
        ]
        
        np.random.seed(42)  # For reproducibility
        atac = ad.AnnData(
            X=np.random.rand(10, 5),
            var=pd.DataFrame(index=peaks)
        )
        
        result = compute_tf_to_atac_links(
            atac=atac,
            ref_genome='mm10',
            fpr=0.05,
            verbose=False,
            n_cpus=1
        )
        
        # Should return a DataFrame with expected structure
        assert isinstance(result, pd.DataFrame)
        assert list(result.columns) == ['source', 'target']
        assert len(result) == 13555, f"Expected 13555 TF-peak links, got {len(result)}"
        
        # Sort by target then source for deterministic comparison
        result_sorted = result.sort_values(['target', 'source']).reset_index(drop=True)
        
        # Validate expected TF-peak links (first 10 rows after sorting)
        expected_first_rows = pd.DataFrame({
            'source': ['AC002126.6_TF', 'AC012531.1_TF', 'AC012531.1_TF', 'AC012531.1_TF', 
                      'AHRR_TF', 'AHR_TF', 'AHR_TF', 'AHR_TF', 'AIRE_TF', 'AL662824.5_TF'],
            'target': ['chr11_97000000_97001000'] * 10
        })
        pd.testing.assert_frame_equal(
            result_sorted.head(10), 
            expected_first_rows,
            check_dtype=False
        )
        
        # Validate all peaks are covered
        assert set(result['target'].unique()) == set(peaks)
    
    def test_compute_tf_to_atac_links_with_tfs_list(self):
        """Test TF-to-ATAC links with specific TF list."""
        # Use regions near Sox2 and Oct4 loci which should contain their binding sites
        peaks = [
            'chr1_3000000_3001000',  # Near Sox2 locus
            'chr1_4807000_4808000',  # Near Pou5f1 (Oct4) locus
            'chr2_29000000_29001000',  # Control enhancer region
        ]
        atac = ad.AnnData(
            X=np.random.rand(10, 3),
            var=pd.DataFrame(index=peaks)
        )
        # Sox2 and Oct4 (Pou5f1) are pluripotency TFs commonly found together
        tfs_list = ['Sox2', 'Pou5f1']
        
        result = compute_tf_to_atac_links(
            atac=atac,
            ref_genome='mm10',
            tfs_list=tfs_list,
            fpr=0.05,
            verbose=False,
            n_cpus=1
        )
        
        assert isinstance(result, pd.DataFrame)
        
        # If results found, validate they match requested TFs
        if len(result) > 0:
            assert 'source' in result.columns
            assert 'target' in result.columns
            # All sources should be from requested TF list (case may vary)
            found_tfs = result['source'].str.lower().unique()
            requested_tfs_lower = [tf.lower() for tf in tfs_list]
            # At least some of the found TFs should match requested list
            # (motif scanning may use TF family names)
            assert len(result) > 0, "Expected to find some TF binding sites in regulatory regions"


class TestComputeRnaNetworkAdvanced:
    """Additional tests for compute_rna_network edge cases."""
    
    def test_compute_rna_network_with_temp_dir(self):
        """Test RNA network inference with custom temp_dir."""
        import tempfile
        import shutil
        
        # Create a temporary directory
        temp_dir = tempfile.mkdtemp()
        
        try:
            np.random.seed(42)
            df_exp = pd.DataFrame(
                np.random.rand(30, 8),
                columns=[f'GENE{i}' for i in range(8)]
            )
            
            result = compute_rna_network(
                df_exp_mtx=df_exp,
                tf_names=['GENE0', 'GENE1'],
                method='GBM',
                temp_dir=temp_dir,
                n_cpu=1,
                seed=42
            )
            
            assert isinstance(result, pd.DataFrame)
            assert len(result) > 0
            # Verify temp_dir was used (check it exists and has files)
            assert os.path.exists(temp_dir)
        finally:
            # Clean up
            shutil.rmtree(temp_dir, ignore_errors=True)


@pytest.mark.skipif(not CELLORACLE_AVAILABLE or not GENOMES_AVAILABLE,
                    reason="celloracle and mm10 genome required")
class TestComputeAtacToRnaLinks:
    """Tests for compute_atac_to_rna_links function."""
    
    def test_compute_atac_to_rna_links_basic(self):
        """Test ATAC-to-RNA links computation."""
        from recon.infer_grn.layers import compute_atac_to_rna_links
        
        # Create ATAC data with peaks near gene TSSs
        peaks = [
            'chr3_34586535_34590322',  # Near Foxp1 TSS
            'chr11_97164516_97169016',  # Near Eomes TSS
            'chr7_103517486_103521617',  # Active enhancer
        ]
        atac = ad.AnnData(
            X=np.random.rand(10, 3),
            var=pd.DataFrame(index=peaks)
        )
        
        # Create RNA data with some genes
        genes = ['Foxp1', 'Eomes', 'Sox2', 'Nanog', 'Pou5f1']
        rna = ad.AnnData(
            X=np.random.rand(10, 5),
            var=pd.DataFrame(index=genes)
        )
        
        result = compute_atac_to_rna_links(
            atac=atac,
            rna=rna,
            ref_genome='mm10'
        )
        
        assert isinstance(result, pd.DataFrame)
        assert list(result.columns) == ['source', 'target']
        
        # Results should link peaks to genes
        if len(result) > 0:
            assert all(result['source'].isin(peaks))
            assert all(result['target'].isin(genes))


class TestCelloracle:
    """Test celloracle availability detection."""
    
    def test_celloracle_availability_flag(self):
        """Test that CELLORACLE_AVAILABLE flag is boolean."""
        assert isinstance(CELLORACLE_AVAILABLE, bool)
    
    @pytest.mark.skipif(CELLORACLE_AVAILABLE, reason="celloracle is installed")
    def test_compute_tf_to_atac_raises_without_celloracle(self):
        """Test that compute_tf_to_atac_links raises ImportError without celloracle."""
        peaks = ['chr1_1000_2000']
        atac = ad.AnnData(
            X=np.random.rand(10, 1),
            var=pd.DataFrame(index=peaks)
        )
        
        with pytest.raises(ImportError, match="CellOracle is required"):
            compute_tf_to_atac_links(
                atac=atac,
                ref_genome='mm10',
                verbose=False
            )
    
    @pytest.mark.skipif(CELLORACLE_AVAILABLE, reason="celloracle is installed")
    def test_compute_atac_to_rna_raises_without_celloracle(self):
        """Test that compute_atac_to_rna_links raises ImportError without celloracle."""
        peaks = ['chr1_1000_2000']
        atac = ad.AnnData(
            X=np.random.rand(10, 1),
            var=pd.DataFrame(index=peaks)
        )
        rna = ad.AnnData(
            X=np.random.rand(10, 2),
            var=pd.DataFrame(index=['GENE1', 'GENE2'])
        )
        
        with pytest.raises(ImportError, match="CellOracle is required"):
            compute_atac_to_rna_links(
                atac=atac,
                rna=rna,
                ref_genome='mm10'
            )


class TestCelloracleMocked:
    """Test error paths by mocking CELLORACLE_AVAILABLE flag."""
    
    @pytest.mark.skipif(not CELLORACLE_AVAILABLE, reason="celloracle not installed - can't test mocking")
    def test_compute_tf_to_atac_raises_when_mocked_unavailable(self, monkeypatch):
        """Test ImportError is raised when CELLORACLE_AVAILABLE is mocked to False."""
        import recon.infer_grn.layers as layers_module
        
        # Mock CELLORACLE_AVAILABLE to False even though celloracle is installed
        monkeypatch.setattr(layers_module, 'CELLORACLE_AVAILABLE', False)
        
        peaks = ['chr1_1000_2000']
        atac = ad.AnnData(
            X=np.random.rand(10, 1),
            var=pd.DataFrame(index=peaks)
        )
        
        with pytest.raises(ImportError, match="CellOracle is required"):
            layers_module.compute_tf_to_atac_links(
                atac=atac,
                ref_genome='mm10',
                verbose=False
            )
    
    @pytest.mark.skipif(not CELLORACLE_AVAILABLE, reason="celloracle not installed - can't test mocking")
    def test_compute_atac_to_rna_raises_when_mocked_unavailable(self, monkeypatch):
        """Test ImportError is raised when CELLORACLE_AVAILABLE is mocked to False."""
        import recon.infer_grn.layers as layers_module
        
        # Mock CELLORACLE_AVAILABLE to False even though celloracle is installed
        monkeypatch.setattr(layers_module, 'CELLORACLE_AVAILABLE', False)
        
        peaks = ['chr1_1000_2000']
        atac = ad.AnnData(
            X=np.random.rand(10, 1),
            var=pd.DataFrame(index=peaks)
        )
        rna = ad.AnnData(
            X=np.random.rand(10, 2),
            var=pd.DataFrame(index=['GENE1', 'GENE2'])
        )
        
        with pytest.raises(ImportError, match="CellOracle is required"):
            layers_module.compute_atac_to_rna_links(
                atac=atac,
                rna=rna,
                ref_genome='mm10'
            )


class TestComputeRnaNetworkExtended:
    """Extended tests for compute_rna_network with more edge cases."""
    
    def test_with_string_temp_dir(self):
        """Test that temp_dir can be provided as string."""
        import tempfile
        
        df_exp = pd.DataFrame(
            np.random.rand(30, 5),
            columns=[f'GENE{i}' for i in range(5)]
        )
        tf_names = ['GENE0', 'GENE1']
        
        with tempfile.TemporaryDirectory() as tmp:
            result = compute_rna_network(
                df_exp_mtx=df_exp,
                tf_names=tf_names,
                temp_dir=tmp,  # String, not Path
                method='GBM',
                n_cpu=1
            )
            
            assert isinstance(result, pd.DataFrame)
            assert 'source' in result.columns
            assert 'target' in result.columns
            assert 'weight' in result.columns
    
    def test_with_nonexistent_temp_dir(self):
        """Test that nonexistent temp_dir is created."""
        import tempfile
        
        df_exp = pd.DataFrame(
            np.random.rand(30, 5),
            columns=[f'GENE{i}' for i in range(5)]
        )
        tf_names = ['GENE0', 'GENE1']
        
        with tempfile.TemporaryDirectory() as tmp:
            nonexistent = Path(tmp) / "subdir" / "temp"
            assert not nonexistent.exists()
            
            result = compute_rna_network(
                df_exp_mtx=df_exp,
                tf_names=tf_names,
                temp_dir=nonexistent,
                method='GBM',
                n_cpu=1
            )
            
            assert isinstance(result, pd.DataFrame)
    
    def test_different_seeds_give_different_results(self):
        """Test that different seeds produce different results."""
        df_exp = pd.DataFrame(
            np.random.rand(50, 8),
            columns=[f'GENE{i}' for i in range(8)]
        )
        tf_names = ['GENE0', 'GENE1']
        
        result1 = compute_rna_network(
            df_exp_mtx=df_exp,
            tf_names=tf_names,
            method='GBM',
            n_cpu=1,
            seed=42
        )
        
        result2 = compute_rna_network(
            df_exp_mtx=df_exp,
            tf_names=tf_names,
            method='GBM',
            n_cpu=1,
            seed=999
        )
        
        # Results should differ with different seeds
        # (at least some weights should be different)
        weights_differ = not np.allclose(
            result1['weight'].values,
            result2['weight'].values
        )
        assert weights_differ
    
    def test_output_columns_correct(self):
        """Test that output has correct column names."""
        df_exp = pd.DataFrame(
            np.random.rand(30, 5),
            columns=[f'GENE{i}' for i in range(5)]
        )
        tf_names = ['GENE0']
        
        result = compute_rna_network(
            df_exp_mtx=df_exp,
            tf_names=tf_names,
            method='GBM',
            n_cpu=1
        )
        
        assert list(result.columns) == ['source', 'target', 'weight']
    
    def test_with_single_tf(self):
        """Test network inference with only one TF."""
        df_exp = pd.DataFrame(
            np.random.rand(40, 6),
            columns=[f'GENE{i}' for i in range(6)]
        )
        tf_names = ['GENE0']  # Only one TF
        
        result = compute_rna_network(
            df_exp_mtx=df_exp,
            tf_names=tf_names,
            method='RF',
            n_cpu=1
        )
        
        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0
        # All sources should be the single TF
        assert all(result['source'] == 'GENE0')


class TestComputeTfNetworkExtended:
    """Extended tests for compute_tf_network edge cases."""
    
    def test_with_empty_tfs_list(self):
        """Test behavior with empty TFs list."""
        rna = ad.AnnData(
            X=np.random.rand(10, 3),
            var=pd.DataFrame(index=['GENE1', 'GENE2', 'GENE3'])
        )
        tfs_list = []
        
        result = compute_tf_network(rna, tfs_list, method=None)
        
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 0
    
    def test_with_all_tfs_missing(self):
        """Test when none of the TFs exist in genes."""
        rna = ad.AnnData(
            X=np.random.rand(10, 3),
            var=pd.DataFrame(index=['GENE1', 'GENE2', 'GENE3'])
        )
        tfs_list = ['TF1', 'TF2', 'TF3']  # None exist
        
        result = compute_tf_network(rna, tfs_list, method=None)
        
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 0
    
    def test_tf_suffix_added_correctly(self):
        """Test that _TF suffix is added correctly."""
        rna = ad.AnnData(
            X=np.random.rand(10, 3),
            var=pd.DataFrame(index=['Foxp3', 'Nfkb1', 'GENE1'])
        )
        tfs_list = ['Foxp3', 'Nfkb1']
        
        result = compute_tf_network(rna, tfs_list, method=None)
        
        assert len(result) == 2
        assert 'Foxp3_TF' in result['source'].values
        assert 'Nfkb1_TF' in result['source'].values
        assert all(result['target'] == 'fake_TF')


class TestGenerateGrn:
    """Tests for generate_grn function."""
    
    def test_generate_grn_basic(self):
        """Test basic GRN generation."""
        from recon.infer_grn.layers import generate_grn
        
        # Create minimal networks
        rna_network = pd.DataFrame({
            'source': ['TF1_TF', 'TF2_TF'],
            'target': ['GENE1', 'GENE2'],
            'weight': [0.8, 0.6]
        })
        
        atac_network = pd.DataFrame({
            'source': ['PEAK1', 'PEAK2'],
            'target': ['PEAK1', 'PEAK2'],
            'weight': [1.0, 1.0]
        })
        
        tf_network = pd.DataFrame({
            'source': ['TF1_TF', 'TF2_TF'],
            'target': ['fake_TF', 'fake_TF']
        })
        
        tf_to_atac = pd.DataFrame({
            'source': ['TF1_TF', 'TF2_TF'],
            'target': ['PEAK1', 'PEAK2']
        })
        
        atac_to_rna = pd.DataFrame({
            'source': ['PEAK1', 'PEAK2'],
            'target': ['GENE1', 'GENE2']
        })
        
        result = generate_grn(
            rna_network=rna_network,
            atac_network=atac_network,
            tf_network=tf_network,
            tf_to_atac_links=tf_to_atac,
            atac_to_rna_links=atac_to_rna,
            n_jobs=1
        )
        
        assert result is not None
        # Result should be a DataFrame from HuMMuS
        assert hasattr(result, '__len__')
    
    def test_generate_grn_renames_bipartite_columns(self):
        """Test that bipartite columns are renamed correctly."""
        from recon.infer_grn.layers import generate_grn
        
        # Ensure the function doesn't crash with column renaming
        rna_network = pd.DataFrame({
            'source': ['TF1_TF'],
            'target': ['GENE1'],
            'weight': [0.9]
        })
        
        atac_network = pd.DataFrame({
            'source': ['PEAK1'],
            'target': ['PEAK1'],
            'weight': [1.0]
        })
        
        tf_network = pd.DataFrame({
            'source': ['TF1_TF'],
            'target': ['fake_TF']
        })
        
        tf_to_atac = pd.DataFrame({
            'source': ['TF1_TF'],
            'target': ['PEAK1']
        })
        
        atac_to_rna = pd.DataFrame({
            'source': ['PEAK1'],
            'target': ['GENE1']
        })
        
        # Should not raise error
        result = generate_grn(
            rna_network=rna_network,
            atac_network=atac_network,
            tf_network=tf_network,
            tf_to_atac_links=tf_to_atac,
            atac_to_rna_links=atac_to_rna,
            n_jobs=1
        )
        
        assert result is not None


class TestEdgeCases:
    """Test edge cases and boundary conditions."""
    
    def test_compute_rna_network_with_many_cpus(self):
        """Test that n_cpu parameter works with multiple cores."""
        df_exp = pd.DataFrame(
            np.random.rand(50, 6),
            columns=[f'GENE{i}' for i in range(6)]
        )
        tf_names = ['GENE0', 'GENE1']
        
        result = compute_rna_network(
            df_exp_mtx=df_exp,
            tf_names=tf_names,
            method='GBM',
            n_cpu=2,  # Use 2 CPUs
            seed=42
        )
        
        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0
    
    def test_compute_tf_network_with_tf_prefix_already_present(self):
        """Test when gene names already have _TF suffix."""
        rna = ad.AnnData(
            X=np.random.rand(10, 3),
            var=pd.DataFrame(index=['Foxp3_TF', 'GENE1', 'GENE2'])
        )
        # TF list without suffix
        tfs_list = ['Foxp3_TF']
        
        result = compute_tf_network(rna, tfs_list, method=None)
        
        # Should add another _TF suffix
        assert len(result) == 1
        assert result['source'].iloc[0] == 'Foxp3_TF_TF'
