# for python 3.10
# pip install circe-py pyarrow "git+https://github.com/cantinilab/HuMMuS.git@dask_update#subdirectory=hummuspy" arboreto
import os
import pathlib
import tempfile
import joblib
from typing import List, Literal, Union

from hummuspy.loader import load_tfs  # import tf names for grnboost2
import hummuspy as hummus
from arboreto.algo import _prepare_input
from arboreto.core import (EARLY_STOP_WINDOW_LENGTH, RF_KWARGS, SGBM_KWARGS,
                           infer_partial_network, to_tf_matrix)
import circe as ci # atac layer

# CellOracle is optional (only needed for GRN inference with ATAC)
try:
    from celloracle import motif_analysis as ma
    CELLORACLE_AVAILABLE = True
except ImportError:
    ma = None
    CELLORACLE_AVAILABLE = False

from tqdm import tqdm
import numpy as np
import pandas as pd
import anndata as ad


def compute_tf_network(
    rna,
    tfs_list,
    method=None
):

    genes = rna.var_names
    if method is None:
        tfs = [tf+"_TF" for tf in tfs_list if tf in genes]
        tf_network = pd.DataFrame([tfs, tfs], index=["source", "target"]).T
        tf_network["target"] = "fake_TF"
        return tf_network
    else:
        raise ValueError("""For now, no method has been implemented.
                   Do not precise a method, or use another function.""")


def compute_rna_network(
        df_exp_mtx: Union[pd.DataFrame, ad.AnnData],
        tf_names: List[str],
        temp_dir: Union[pathlib.Path, None] = None,
        method: Literal['GBM', 'RF'] = 'GBM',
        n_cpu: int = 1,
        seed: int = 666) -> pd.DataFrame:
    """
    # Inspired from SCENICPLUS:
    https://github.com/aertslab/scenicplus/blob/main/src/scenicplus/TF_to_gene.py

    Calculate TF-to-gene relationships using either
    Gradient Boosting Machine (GBM) or Random Forest (RF) regression.

    It is a wrapper around the `infer_partial_network` function
    from the arboreto package, similarly to GRNBoost2.
    It uses joblib to parallelize the inference of the relationships
    for each target gene.
    It returns a DataFrame with the TF-to-gene relationships and
    their importance scores.

    Parameters
    ----------
    df_exp_mtx : pd.DataFrame, ad.AnnData
        Gene expression matrix with genes as columns and cells as rows.
        If an AnnData object is provided, the expression matrix is extracted
        using the `to_df()` method.
    tf_names : List[str]
        List of transcription factor names to consider as potential regulators.
    temp_dir : pathlib.Path
        Path to a temporary directory to store intermediate files during
        parallel processing.
        If None, a temporary directory will be created and deleted after use.
    method : Literal['GBM', 'RF'], optional
        Method to use for regression. Either 'GBM' for Gradient
        Boosting Machine or 'RF' for Random Forest. Default is 'GBM'.
    n_cpu : int, optional
        Number of CPU cores to use for parallel processing. Default is 1.
    seed : int, optional
        Random seed for reproducibility. Default is 666.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns ['tf', 'target', 'importance'] 
        representing the TF-to-gene relationships and their importance scores.

    """

    if isinstance(df_exp_mtx, ad.AnnData):
        df_exp_mtx = df_exp_mtx.to_df()
    elif not isinstance(df_exp_mtx, pd.DataFrame):
        raise ValueError(
            "df_exp_mtx must be a pandas DataFrame or an AnnData object.")

    if temp_dir is None:
        with tempfile.TemporaryDirectory() as tmp:
            temp_dir = pathlib.Path(tmp)
            # use temp_dir here
    else:
        temp_dir = pathlib.Path(temp_dir)
        temp_dir.mkdir(parents=True, exist_ok=True)

    if (method == 'GBM'):
        method_params = [
            'GBM',      # regressor_type
            SGBM_KWARGS  # regressor_kwargs
        ]
    elif (method == 'RF'):
        method_params = [
            'RF',       # regressor_type
            RF_KWARGS   # regressor_kwargs
        ]

    exp_mtx, gene_names, tf_names = _prepare_input(
        expression_data=df_exp_mtx, gene_names=None, tf_names=tf_names)
    tf_matrix, tf_matrix_gene_names = to_tf_matrix(
        exp_mtx,  gene_names, tf_names)

    print('Calculating TF-to-gene importance')
    if temp_dir is not None:
        if type(temp_dir) is str:
            temp_dir = pathlib.Path(temp_dir)
        if not temp_dir.exists():
            Warning(f"{temp_dir} does not exist, creating it.")
            os.makedirs(temp_dir)

    TF_to_genes = joblib.Parallel(
        n_jobs=n_cpu,
        temp_folder=temp_dir)(
            joblib.delayed(infer_partial_network)(
                target_gene_name=gene,
                target_gene_expression=exp_mtx[:, gene_names.index(gene)],
                regressor_type=method_params[0],
                regressor_kwargs=method_params[1],
                tf_matrix=tf_matrix,
                tf_matrix_gene_names=tf_matrix_gene_names,
                include_meta=False,
                early_stop_window_length=EARLY_STOP_WINDOW_LENGTH,
                seed=seed)
            for gene in tqdm(
                gene_names,
                total=len(gene_names),
                desc=f'Running using {n_cpu} cores'))

    adj = pd.concat(TF_to_genes).sort_values(by='importance', ascending=False)

    return adj[["TF", "target", "importance"]].rename(columns={
        "TF": "source",
        "target": "target",
        "importance": "weight"})


def compute_tf_to_atac_links(
    atac,
    ref_genome,
    tfs_list: Union[List[str], None] = None,
    genomes_dir=None,
    motifs=None,
    fpr=0.02,
    verbose=True,
    indirect=True,
    n_cpus=-1
):

    """
    Compute TF-to-ATAC peak links using motif scanning.
    It uses the CellOracle motif_analysis module to scan for motifs
    in the provided ATAC peaks.
    It returns a DataFrame with the TF-to-peak links.

    Parameters
    ----------
    atac : anndata.AnnData
        AnnData object containing the ATAC-seq data.
        The peak names should be in the format 'chr_start_end'.
    ref_genome : str
        Reference genome to use for motif scanning.
        E.g., 'hg38', 'mm10', etc.
    genomes_dir : str, optional
        Directory containing the reference genomes.
        If None, the default CellOracle genomes directory will be used.
    motifs : list, optional
        List of motifs to use for scanning.
        If None, the default CellOracle motifs will be used.
    fpr : float, optional
        False positive rate for motif scanning.
        Default is 0.02.
    verbose : bool, optional
        Whether to print progress messages.
        Default is True.
    indirect : bool, optional
        Whether to include TF-to-peak links from indirect evidences.
        Default is True.
    n_cpus : int, optional
        Number of CPUs to use for parallel processing.
        Default is -1 (use all available CPUs).

    Returns
    -------
    pd.DataFrame
        DataFrame with columns ['source', 'target'] representing the
        TF-to-peak links.
    """
    
    if not CELLORACLE_AVAILABLE:
        raise ImportError(
            "CellOracle is required for ATAC-seq analysis. "
            "Install with: pip install 'git+https://github.com/cantinilab/celloracle@lite'"
        )

    # Create a DataFrame around the peak ids to run CellOracle motif search
    peaks_df = pd.concat([
        pd.Series(atac.var_names),
        pd.Series(['None']*len(atac.var_names))
        ], axis=1)
    peaks_df.columns = ['peak_id', 'gene_short_name']

    # remove wrongly formatted peaks
    peaks_df = ma.check_peak_format(peaks_df, ref_genome, genomes_dir=None)

    # Instantiate TFinfo object
    tfi = ma.TFinfo(
        peak_data_frame=peaks_df,
        ref_genome=ref_genome,
        genomes_dir=genomes_dir)

    # Scan motifs. !CAUTION! This step may take several hours with many peaks!
    tfi.scan(fpr=fpr,
             motifs=motifs,  # If None, default motifs will be loaded.
             verbose=verbose,
             n_cpus=n_cpus)

    # Reset filtering
    tfi.reset_filtering()

    # Do filtering
    tfi.filter_motifs_by_score(threshold=10)

    # Format post-filtering results.
    tfi.make_TFinfo_dataframe_and_dictionary(verbose=True)

    # Split seqname into a peak column
    tfi.scanned_df = tfi.scanned_df.rename(
        columns={"seqname": "peak"})

    # Split the factor lists into Python lists
    if isinstance(tfi.scanned_df["factors_direct"].iloc[0], str):
        tfi.scanned_df["factors_direct"] = tfi.scanned_df[
            "factors_direct"].dropna().str.split(r",\s*")
    if isinstance(tfi.scanned_df["factors_indirect"].iloc[0], str):
        tfi.scanned_df["factors_indirect"] = tfi.scanned_df[
            "factors_indirect"].dropna().str.split(r",\s*")

    # Explode both direct and indirect into long format
    direct = tfi.scanned_df.explode("factors_direct")[[
        "peak", "factors_direct"]].rename(columns={"factors_direct": "factor"})

    if indirect:
        indirect = tfi.scanned_df.explode("factors_indirect")[[
            "peak", "factors_indirect"]].rename(columns={
                "factors_indirect": "factor"})
        # Combine
        result = pd.concat([direct, indirect]).dropna().reset_index(drop=True)
        if tfs_list is not None:
            result = result[result["factor"].isin(tfs_list)]

        result["factor"] = result["factor"] + '_TF'
        return result[["factor", "peak"]].rename(
            columns={
                "factor": "source",
                "peak": "target",
                "score": "weight"})

    else:
        direct["factor"] = direct["factor"] + '_TF'
        if tfs_list is not None:
            direct = direct[direct["factor"].str.replace('_TF', '').isin(tfs_list)]
        return direct[["factor", "peak"]].rename(
            columns={
                "factor": "source",
                "peak": "target",
                "score": "weight"})


def compute_atac_to_rna_links(atac, rna, ref_genome):

    """
    Compute ATAC peak-to-gene links using TSS annotation.
    It uses the CellOracle motif_analysis module to get TSS information
    for the provided ATAC peaks.
    It returns a DataFrame with the peak-to-gene links.

    Parameters
    ----------
    atac : anndata.AnnData
        AnnData object containing the ATAC-seq data.
        The peak names should be in the format 'chr_start_end'.
    rna : anndata.AnnData
        AnnData object containing the RNA-seq data.
        The gene names should match those in the ATAC-seq data.
    ref_genome : str
        Reference genome to use for TSS annotation.
        E.g., 'hg38', 'mm10', etc.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns ['source', 'target'] representing the
        peak-to-gene links.
    """
    
    if not CELLORACLE_AVAILABLE:
        raise ImportError(
            "CellOracle is required for ATAC-seq analysis. "
            "Install with: pip install 'git+https://github.com/cantinilab/celloracle@lite'"
        )

    # Get TSS information for the peaks
    peaks_annotated = ma.get_tss_info(
        peak_str_list=atac.var_names.values,
        ref_genome=ref_genome)

    peaks_annotated["peak"] = peaks_annotated["chr"] + "_" +\
        peaks_annotated["start"].astype(str) + "_" +\
        peaks_annotated["end"].astype(str)

    peaks_annotated = peaks_annotated.rename(
        columns={"gene_short_name": "gene"})
    peaks_annotated = peaks_annotated[peaks_annotated["gene"].isin(
        rna.var_names.values)]

    return peaks_annotated[["peak", "gene"]].rename(columns={
        "peak": "source",
        "gene": "target"})


def generate_grn(
    rna_network,
    atac_network,
    tf_network,
    tf_to_atac_links,
    atac_to_rna_links,
    n_jobs=1
):
    
    """
    Generate a Gene Regulatory Network (GRN) by integrating
    TF-to-gene, peak-to-gene, and TF-to-peak relationships.
    It uses the HuMMuS package to create a multiplex network
    and perform random walks to rank the nodes.
    It returns a DataFrame with the ranked nodes.

    Parameters
    ----------
    rna_network : pd.DataFrame
        DataFrame with columns ['source', 'target', 'weight']
        representing the TF-to-gene relationships.
    atac_network : pd.DataFrame
        DataFrame with columns ['source', 'target', 'weight']
        representing the peak-to-gene relationships.
    tf_network : pd.DataFrame
        DataFrame with columns ['source', 'target', 'weight']
        representing the TF-to-TF relationships.
    tf_to_atac_links : pd.DataFrame
        DataFrame with columns ['source', 'target']
        representing the TF-to-peak relationships.
    atac_to_rna_links : pd.DataFrame
        DataFrame with columns ['source', 'target']
        representing the peak-to-gene relationships.
    n_jobs : int, optional
        Number of jobs to use for parallel processing.
        Default is 1.

    Returns
    -------
    pd.DataFrame
        DataFrame with ranked nodes from the GRN.
    """
    # check that networks have source, target and weight

    # modify source and target for bipartites
    tf_to_atac_links = tf_to_atac_links.rename(columns={
        "source": "col1", "target": "col2"})
    atac_to_rna_links = atac_to_rna_links.rename(columns={
        "source": "col1", "target": "col2"})

    config_grn = hummus.core_grn.init_GRN(
        TF_layer=tf_network,
        ATAC_layer=atac_network,
        RNA_layer=rna_network,
        TF_ATAC_bipartite=tf_to_atac_links,
        ATAC_RNA_bipartite=atac_to_rna_links,
    )

    seeds = np.unique(np.concatenate(
        [tf_network["source"].unique(), tf_network["target"].unique()]))
    seeds = [seed for seed in seeds if seed != "fake_TF"]

    ranking = hummus.explore_network.compute_multiple_RandomWalk(
        multiplex=config_grn[0],
        bipartite=config_grn[1],
        eta=config_grn[2],
        lamb=config_grn[3].values.tolist(),
        seeds=seeds,
        save=False,
        n_jobs=n_jobs
    )

    return ranking.compute()  # only if dask branch
