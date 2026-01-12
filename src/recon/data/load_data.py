from importlib.resources import files
import pandas as pd
import os

receptor_gene_resources = [
    "human_receptor_gene_from_NichenetPKN",
    "mouse_receptor_gene_from_NichenetPKN",
]

# Updated Zenodo record URL
TUTORIAL_DATA_URL = "https://zenodo.org/record/18223725/files/"

# Updated tutorial data registry with new Zenodo record
TUTORIAL_DATA_REGISTRY = {
    # Perturbation tutorial data
    "perturbation_tuto/rna.h5ad": "sha256:12be5576beccc26b286dfca8e1ea489a5a9c8f96b003a178022a929bd209af2e",
    "perturbation_tuto/rna_treated.h5ad": "sha256:0fbd658fca102ca24a0b1965e442abe704acd5ce0e26cdbbacd0453367cec42b",
    "perturbation_tuto/grn.csv": "sha256:0d4d7857d5ddbf023326b9f0041c5db9c4c2c0a3720dee2c5dad24aee3e00bf9",
    # GRN inference tutorial data
    "build_grn_tuto/pbmc10x.h5mu": "sha256:b12ab3b142315c297d198274792b8c55d74986c14b76148aba6409a76ae1c23c",
}


def load_receptor_genes(receptor_gene_list) -> "pd.DataFrame":

    if receptor_gene_list not in receptor_gene_resources:
       raise ValueError(f"The name of the receptor gene list must be in {receptor_gene_resources}")

    path = files("recon.data.receptor_genes").joinpath(receptor_gene_list+".parquet")
    return pd.read_parquet(path)


def fetch_tutorial_data(filename: str, data_dir: str = "./data", force: bool = False) -> str:
    """
    Download tutorial data files from Zenodo.
    
    Parameters
    ----------
    filename : str
        Name of the file to download. Available files:
        
        **Perturbation tutorial** (tutorials 1-3):
        - "perturbation_tuto/rna.h5ad": scRNA-seq data (24 MB)
        - "perturbation_tuto/rna_treated.h5ad": treated scRNA-seq data (1.7 GB)
        - "perturbation_tuto/grn.csv": pre-computed GRN (168 MB)
        
        **GRN inference tutorial** (tutorial 4):
        - "build_grn_tuto/pbmc10x.h5mu": multimodal PBMC data (748 MB)
        
    data_dir : str, default="./data"
        Base directory to save the downloaded file.
    force : bool, default=False
        If True, re-download even if file exists.
    
    Returns
    -------
    str
        Path to the downloaded file.
    
    Examples
    --------
    >>> from recon.data import fetch_tutorial_data
    >>> # Perturbation tutorial
    >>> rna_path = fetch_tutorial_data("perturbation_tuto/rna.h5ad")
    >>> import scanpy as sc
    >>> rna = sc.read_h5ad(rna_path)
    >>> 
    >>> # GRN inference tutorial
    >>> mdata_path = fetch_tutorial_data("build_grn_tuto/pbmc10x.h5mu")
    >>> import muon as mu
    >>> mdata = mu.read(mdata_path)
    """
    try:
        import pooch
    except ImportError:
        raise ImportError(
            "pooch is required to download tutorial data. "
            "Install it with: pip install pooch"
        )
    
    if filename not in TUTORIAL_DATA_REGISTRY:
        raise ValueError(
            f"Unknown file: {filename}. "
            f"Available files: {list(TUTORIAL_DATA_REGISTRY.keys())}"
        )
    
    # Create full path including subdirectory
    filepath = os.path.join(data_dir, filename)
    filedir = os.path.dirname(filepath)
    os.makedirs(filedir, exist_ok=True)
    
    # Check if file already exists
    if os.path.exists(filepath) and not force:
        print(f"File already exists: {filepath}")
        return filepath
    
    # Download the file (use basename for Zenodo URL)
    basename = os.path.basename(filename)
    print(f"Downloading {filename} from Zenodo...")
    url = TUTORIAL_DATA_URL + basename
    known_hash = TUTORIAL_DATA_REGISTRY[filename]
    
    downloaded_path = pooch.retrieve(
        url=url,
        known_hash=known_hash,
        fname=filename,
        path=data_dir,
        progressbar=True
    )
    
    print(f"Downloaded to: {downloaded_path}")
    return downloaded_path


def fetch_all_tutorial_data(data_dir: str = "./data/perturbation_tuto", force: bool = False) -> dict:
    """
    Download all tutorial data files.
    
    Parameters
    ----------
    data_dir : str, default="./data/perturbation_tuto"
        Directory to save the downloaded files.
    force : bool, default=False
        If True, re-download even if files exist.
    
    Returns
    -------
    dict
        Dictionary mapping filenames to their local paths.
    """
    paths = {}
    for filename in TUTORIAL_DATA_REGISTRY:
        paths[filename] = fetch_tutorial_data(filename, data_dir=data_dir, force=force)
    return paths