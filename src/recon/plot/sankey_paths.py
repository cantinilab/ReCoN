import pandas as pd
import networkx as nx
from typing import Union, List, Tuple
import plotly.graph_objects as go
import hashlib


def get_celltype_grn_receptor_bipartite(
    multicell_obj,
    cell_type: str,
    as_dataframe: bool = True
):
    """
    Retrieve the “{cell_type}_grn-{cell_type}_receptor” bipartite edges from multicell_obj.bipartites.

    Parameters
    ----------
    multicell_obj : any
        An object with attribute .bipartites, which is a dict whose keys include
        "{cell_type}_grn-{cell_type}_receptor" for various cell types.
    cell_type : str
        The exact cell-type name, e.g. "Endothelial cell" or "Fibroblast".
    as_dataframe : bool, default=True
        If True, return the raw pandas DataFrame ...
        If False, return a NetworkX Graph ...

    Returns
    -------
    pandas.DataFrame if as_dataframe=True
    networkx.Graph     if as_dataframe=False

    Raises
    ------
    KeyError
        If no key of the form "{cell_type}_grn-{cell_type}_receptor" is found.
    """
    # Build the expected bipartites key
    bipartite_key = f"{cell_type}_grn-{cell_type}_receptor"

    # Try to fetch that entry from multicell_obj.bipartites
    bip_info = multicell_obj.bipartites.get(bipartite_key)
    if bip_info is None:
        raise KeyError(f"No bipartite entry found for key '{bipartite_key}' "
                       f"in multicell_obj.bipartites. Available cell types are:\n"
                       f"{multicell_obj.celltypes_names}")

    # Extract the DataFrame
    edge_df = bip_info["edge_list_df"]

    if as_dataframe:
        return edge_df.copy()

    # Otherwise, build a NetworkX Graph from that DataFrame
    G = nx.Graph()
    # We assume edge_df has columns 'col1' and 'col2'; copy 'network_key' (i.e. name of this layer) as edge-attrib
    for _, row in edge_df.iterrows():
        u = row["col1"]
        v = row["col2"]
        data = {}
        if "network_key" in row.index:
            data["layer"] = row["network_key"]
        G.add_edge(u, v, **data)
    return G


def get_celltype_gene_layer(
    multicell_obj,
    cell_type: str,
    layer_name: str = "grn",
    as_dataframe: bool = True
) -> Union[pd.DataFrame, nx.DiGraph, nx.Graph]:
    """
    Retrieve a single sublayer (by its multiplex key) from multicell_obj.multiplexes.
    The multiplex key is constructed as "{cell_type}_{multiplex_suffix}".  We assume
    that this key corresponds to exactly one entry in .multiplexes, and that entry’s
    "layers" list contains one DataFrame (or, if multiple, we always take layers[0]).

    Parameters
    ----------
    multicell_obj : any
        An object with attribute .multiplexes, a dict whose keys include
        "{cell_type}_{multiplex_suffix}".
    cell_type : str
        The cell‐type portion of the key, e.g. "Endothelial cell" or "Fibroblast".
    layer_name : str, default="grn"
        The suffix naming the sublayer, e.g. "grn", selected inside the {celltype} gene multiplex.
    as_dataframe : bool, default=True
        If True, return the raw pandas DataFrame (layers[0]) for that key.
        If False, convert that DataFrame into a NetworkX Graph (directed if
        graph_type[0] is "01" or "10", otherwise undirected).

    Returns
    -------
    pandas.DataFrame
        If as_dataframe=True, returns a copy of the DataFrame under
        `multiplexes["{cell_type}_{multiplex_suffix}"]["layers"][0]`.
    networkx.DiGraph or networkx.Graph
        If as_dataframe=False, returns a Graph built from that DataFrame’s
        columns ["source", "target", "weight"], directed if graph_type[0]
        is in ("01","10"), else undirected.

    Raises
    ------
    KeyError
        If the key "{cell_type}_grn" is not present in .multiplexes.
    """
    # 1. Construct the multiplex key
    multiplex_key = f"{cell_type}_grn"

    if multiplex_key not in multicell_obj.multiplexes:
        available = list(multicell_obj.multiplexes.keys())
        raise KeyError(
            f"No multiplex found for key '{multiplex_key}'. Available keys:\n"
            f"{available}"
        )

    mux_info = multicell_obj.multiplexes[multiplex_key]

    # We assume "layers" is a list of DataFrames; take the first one
    layers_list = mux_info.get("layers", [])
    if not layers_list:
        raise KeyError(f"Multiplex '{multiplex_key}' has no layers listed.")

    # Step 1: build a dict mapping each name → its DataFrame
    name_to_layer = dict(zip(mux_info["names"], mux_info["layers"]))
    # Step 2: extract just the DataFrames whose names are in names_of_interest
    if layer_name in name_to_layer:
        layer = name_to_layer[layer_name] 
    else:
        raise KeyError(f"Multiplex '{multiplex_key}' has no layer {layer_name} listed.")

    if as_dataframe:
        return layer

    # Otherwise, build a NetworkX graph from df
    # Decide directed vs. undirected based on graph_type[0]
    graph_types = mux_info.get("graph_type", [])
    gt = graph_types[0] if graph_types else "00"
    is_directed = gt in ("01", "10")

    G = nx.DiGraph() if is_directed else nx.Graph()

    # DataFrame columns: ["source", "target", "weight", ...]
    for _, row in layer.iterrows():
        u = row["source"]
        v = row["target"]
        w = row.get("weight", 1.0)
        G.add_edge(u, v, weight=w)
        # If the multiplex is undirected but we created a DiGraph (because gt in ("01","10")?),
        # this code already only adds (u→v) for directed; for undirected we used nx.Graph() above.

    return G


def get_cell_communication_layer(
    multicell_obj,
    receptor_cells: Union[list, None] = None,
    ligand_cells: Union[list, None] = None,
    as_dataframe: Union[bool, None] = True
) -> Union[pd.DataFrame, nx.DiGraph, nx.Graph]:
    """
    Retrieve the “cell_communication” multiplex from multicell_obj.multiplexes.

    Parameters
    ----------
    multicell_obj : any
        An object with attribute .multiplexes, which is a dict whose keys include
        "cell_communication".
    as_dataframe : bool, default=True
        If True, return the raw pandas DataFrame (layers[0]).
        If False, return a networkx.Graph built from that DataFrame.

    Returns
    -------
    pandas.DataFrame if as_dataframe=True
    networkx.DiGraph or networkx.Graph if as_dataframe=False

    Raises
    ------
    KeyError
        If no key "cell_communication" is found in multicell_obj.multiplexes,
        or if its "layers" list is empty.
    receptor_cells : list of str, optional
        If provided, only keep pairs where the receptor’s “cell‐type suffix” is in this list.
        For example ["Fibroblast", "Endothelial cell"]. If None, do not filter on receptor‐celltype.
    ligand_cells : list of str, optional
        If provided, only keep pairs where the ligand’s “cell‐type suffix” is in this list.
        If None, do not filter on ligand‐celltype.

    Returns
    -------
    pandas.DataFrame
        A filtered DataFrame of receptor–ligand pairs, with columns at least:
          ["ligand", "receptor", "receptor_clean", "weight", …],
        where:
          • "ligand"         is "GENE::CellType"
          • "receptor"       is "GENE::CellType"
          • "receptor_clean" is "GENE_receptor::CellType"

    """
    key = "cell_communication"
    if key not in multicell_obj.multiplexes:
        available = list(multicell_obj.multiplexes.keys())
        raise KeyError(
            f"No multiplex found for key '{key}'. Available keys:\n{available}"
        )

    mux_info = multicell_obj.multiplexes[key]
    layers_list = mux_info.get("layers", [])
    if not layers_list:
        raise KeyError(f"Multiplex '{key}' has no layers listed.")

    df = layers_list[0]  # assume the first (and only) DataFrame is the cell_communication layer
    df = df.loc[:, ["source", "target", "weight", "celltype_source", "celltype_target", "network_key"]]

    # 1) Normalize both source and target to "GENE::CellType" form
    def _standardize_name(x: str) -> str:
        # If name contains '-' and not '::', replace last '-' with '::'
        if "-" in x and "::" not in x:
            gene_part, cell_part = x.rsplit("-", 1)
            return f"{gene_part}::{cell_part}"
        else:
            return x

    df.loc[:, "source_std"] = df.loc[:, "source"].astype(str).apply(_standardize_name)
    df.loc[:, "target_std"] = df.loc[:, "target"].astype(str).apply(_standardize_name)

    # 2) By convention:  source_std → ligand,  target_std → receptor (before adding suffix)
    df.loc[:, "receptor"] = df.loc[:, "source_std"]
    df.loc[:, "ligand"] = df.loc[:, "target_std"]

    # 3) Build receptor_clean by inserting "_receptor" before "::"
    def _add_receptor_suffix(x: str) -> str:
        if "::" in x:
            gene, cell = x.split("::", 1)
            return f"{gene}_receptor::{cell}"
        else:
            # If no "::", just append "_receptor"
            return f"{x}_receptor"
    df["receptor_clean"] = df["receptor"].apply(_add_receptor_suffix)

    # 4) Ensure there is a "weight" column
    if "weight" not in df.columns:
        df["weight"] = 1.0

    # 5) Filter by receptor_cells if provided
    if receptor_cells is not None:
        pattern = "|".join([f"::{rc}$" for rc in receptor_cells])
        df = df.loc[
            df["receptor"].str.contains(pattern, regex=True),
            :
        ]

    # 6) Filter by ligand_cells if provided
    if ligand_cells is not None:
        pattern = "|".join([f"::{lc}$" for lc in ligand_cells])
        df = df.loc[
            df["ligand"].str.contains(pattern, regex=True),
            :
        ]

    if as_dataframe:
        # 7) Return only the columns we care about (plus any extras)
        keep_cols = ["ligand", "receptor", "receptor_clean", "weight"] + [
            c for c in df.columns
            if c not in ("source", "target", "source_std", "target_std", "ligand", "receptor", "receptor_clean", "weight")
        ]
        return df[keep_cols].reset_index(drop=True)

    # Otherwise, build a NetworkX graph
    # Decide directed vs. undirected based on graph_type[0]
    graph_types = mux_info.get("graph_type", [])
    gt = graph_types[0] if graph_types else "00"
    is_directed = gt in ("01", "10")

    G = nx.DiGraph() if is_directed else nx.Graph()
    # Expect columns: ["source", "target", maybe "weight", maybe others...]
    for _, row in df.iterrows():
        u = row["source"]
        v = row["target"]
        w = row.get("weight", 1.0)
        G.add_edge(u, v, weight=w)

    return G


def get_top_tfs(
    results_df: pd.DataFrame,
    cell_type,
    n: int = 5,
) -> pd.DataFrame:
    """
    From a `results_df` of node‐scores, identify the top‐n TF nodes
    for a given cell type.

    Parameters
    ----------
    results_df : pandas.DataFrame
        Must have columns:
          - "multiplex" (e.g. "Endothelial cell_grn")
          - "node"      (e.g. "GENE_TF::Endothelial cell")
          - "score"     (numeric)
    cell_type : str
        The cell type of interest, e.g. "Endothelial cell" or "Fibroblast".
    n : int, default=5
        The number of top TFs to return.

    Returns
    -------
    pandas.DataFrame
        A subset of `results_df` (only rows where `multiplex` corresponds to the
        given cell type and whose `node` ends with `_TF::{cell_type}`),
        sorted by descending score, limited to the top‐n rows.
    """

    results = results_df.copy()
    results.loc[:, "celltype"] = results["multiplex"].str.split("_", expand=True)[0]
    results = results[results["celltype"]==cell_type]

    results = results[results["node"].str.endswith(f"_TF::{cell_type}")]
    return results.iloc[:n, :]


def get_top_receptors(
    results_df: pd.DataFrame,
    cell_type,
    n: int = 5,
) -> pd.DataFrame:
    """
    From a `results_df` of node‐scores, identify the top‐n receptor nodes
    for a given cell type.

    Parameters
    ----------
    results_df : pandas.DataFrame
        Must have columns:
          - "multiplex" (e.g. "cell_communication")
          - "node"      (e.g. "GENE_receptor::Endothelial cell")
          - "score"     (numeric)
    cell_type : str
        The cell type of interest, e.g. "Endothelial cell" or "Fibroblast".
    n : int, default=5
        The number of top receptors to return.

    Returns
    -------
    pandas.DataFrame
        A subset of `results_df` (only rows where `multiplex` corresponds to "cell_communication"
        and whose `node` ends with `_receptor::{cell_type}`), sorted by descending score,
        limited to the top‐n rows.
    """

    results = results_df.copy()
    results.loc[:, "celltype"] = results.loc[:, "multiplex"].str.split("_", expand=True)[0]
    results = results[results["celltype"]==cell_type]
    
    # remove fake node
    results = results[results["node"]!=f"fake_receptor::{cell_type}"]

    results = results[results["node"].str.endswith(f"_receptor::{cell_type}")]
    return results.iloc[:n, :]


def get_top_ligands(
    results_df: pd.DataFrame,
    receptor_ligand_df: pd.DataFrame,
    n: int = 5,
    per_celltype: bool = False
) -> pd.DataFrame:
    """
    From a `results_df` of node‐scores and a `receptor_ligand_df`
    (as returned by extract_receptor_ligand_pairs), identify the top‐n ligand nodes
    (globally or per ligand‐celltype) based on their score in the “cell_communication” multiplex.

    This version will first convert any "GENE-CellType" strings in results_df
    (for multiplex=="cell_communication") into "GENE::CellType" so that they line up
    with the ligand names in `receptor_ligand_df`.

    Parameters
    ----------
    results_df : pandas.DataFrame
        Must have columns:
          - "multiplex" (e.g. "cell_communication")
          - "node"      (e.g. "VEGFA‐Endothelial cell"  or  "VEGFA::Endothelial cell")
          - "score"     (numeric)
    receptor_ligand_df : pandas.DataFrame
        Output of `extract_receptor_ligand_pairs(...)`, containing at least:
          - "ligand"   (already in "GENE::CellType" form)
          - "receptor" (already in "GENE::CellType" form)
    n : int, default=5
        If per_celltype=False: pick the overall top‐n ligands by score.
        If per_celltype=True: pick the top‐n ligands within each distinct ligand‐celltype.
    per_celltype : bool, default=False
        If False, do a global ranking. If True, do one ranking per source cell type.

    Returns
    -------
    pandas.DataFrame
        A subset of `results_df` (only rows where `multiplex=="cell_communication"`
        and whose standardized‐node is in the chosen top ligand set). If `per_celltype=True`,
        this DataFrame will include a new column "ligand_celltype" indicating each
        ligand’s source cell type.

    Raises
    ------
    KeyError
        - If `results_df` is missing any of {"multiplex", "node", "score"}.
        - If `receptor_ligand_df` is missing "ligand".
    """

    # 1) Validate required columns
    missing_res = {"multiplex", "node", "score"}.difference(results_df.columns)
    if missing_res:
        raise KeyError(f"`results_df` is missing required columns: {missing_res}")
    if "ligand" not in receptor_ligand_df.columns:
        raise KeyError("`receptor_ligand_df` must contain a column named 'ligand'.")

    # 2) Restrict to just the "cell_communication" rows
    cc_scores = results_df.loc[results_df["multiplex"] == "cell_communication", :].copy()
    if cc_scores.empty:
        return cc_scores  # nothing to do if that multiplex is absent

    # 3) STANDARDIZE any "GENE-CellType" → "GENE::CellType" in cc_scores["node"]
    def _hyphen_to_doublecolon(x: str) -> str:
        # If there's a hyphen but no "::", split on the last hyphen.
        # E.g. "VEGFA-Endothelial cell" → "VEGFA::Endothelial cell"
        if "-" in x and "::" not in x:
            gene_part, cell_part = x.rsplit("-", 1)
            return f"{gene_part}::{cell_part}"
        else:
            return x

    # Create a new column "node_std" that always uses "::" form
    cc_scores.loc[:, "node_std"] = cc_scores.loc[:, "node"].astype(str).apply(_hyphen_to_doublecolon)

    # 4) Find the set of all ligands that actually appear in receptor_ligand_df
    ligands_in_pairs = set(receptor_ligand_df.loc[:, "ligand"].astype(str).unique())

    # 5) Keep only those rows whose standardized node is an actual ligand
    cc_scores = cc_scores.loc[cc_scores.loc[:,"node_std"].isin(ligands_in_pairs), :].copy()
    if cc_scores.empty:
        return cc_scores  # no overlap → nothing to return

    # 6) Extract the ligand_celltype from the "::" form.  E.g. "VEGFA::Endothelial cell" → "Endothelial cell"
    cc_scores.loc[:, "ligand_celltype"] = cc_scores.loc[:, "node_std"].str.split("::", n=1).str[-1]

    # 7) Decide which top‐n nodes to pick
    if per_celltype:
        #   Within each ligand_celltype, pick the top‐n by "score"
        top_ligand_names = (
            cc_scores
            .sort_values(["ligand_celltype", "score"], ascending=[True, False])
            .groupby("ligand_celltype", group_keys=False)["node_std"]
            .head(n)
            .tolist()
        )
    else:
        #   Globally (ignore celltype): pick the top‐n by "score"
        top_ligand_names = (
            cc_scores
            .sort_values("score", ascending=False)
            .head(n)["node_std"]
            .tolist()
        )

    # 8) Filter cc_scores to only those whose node_std is in our top‐ligand set
    filtered = cc_scores.loc[cc_scores.loc[:, "node_std"].isin(top_ligand_names), :].copy()

    # 9) Return the original columns (plus ligand_celltype).  Drop node_std
    return filtered.drop(columns=["node_std"]).reset_index(drop=True)


def extract_gene_tf_pairs(
    tf_gene_layer, top_tfs, seeds
) -> pd.DataFrame:
    """
    Given a tf_gene_layer DataFrame (from get_celltype_grn_receptor_bipartite or
    get_celltype_gene_layer) and a top_tfs DataFrame (from get_top_tfs),
    return only those rows of tf_gene_layer whose tf and gene both appear in the
    respective top lists.
    This function assumes:
        • tf_gene_layer has columns "source" (tf) and "target" (gene).
        • top_tfs has a column "node" containing tf names in the form "GENE_TF::CellType".
        • seeds is a pd.Index or list of gene names (in the form "GENE::CellType").
    We therefore standardize all three sources to "GENE::CellType" before filtering.
    Parameters
    ----------
    tf_gene_layer : pandas.DataFrame
        Output of get_celltype_grn_receptor_bipartite() or get_celltype
        _gene_layer(), with at least:
            • "source" (tf) e.g. "GENE_TF::Endothelial cell"
            • "target" (gene) e.g. "GENE::Endothelial cell"
    top_tfs : pandas.DataFrame
        Subset of results_df for multiplex=="{cell_type}_grn", with column:
            • "node" e.g. "GENE_TF::Endothelial cell"
    seeds : list or pd.Index
        List of gene names (in the form "GENE::CellType") to filter for.
    Returns
    -------
    pandas.DataFrame
        A filtered DataFrame containing only those rows of tf_gene_layer where:
            tf ∈ standardized(top_tfs["node"])  AND
            gene ∈ seeds.
    """

    sources_list = seeds.tolist()
    targets_list = list(top_tfs["node"].values)

    print(targets_list)
    # (2) Now filter tf_gene_layer for rows where both ‘source’ ∈ sources_list AND ‘target’ ∈ targets_list:
    filtered_df = tf_gene_layer[
        tf_gene_layer["source"].isin(sources_list) &
        tf_gene_layer["target"].isin(targets_list)
    ]
    filtered_df.columns = ["gene", "tf", "weight", "network_key"]
    filtered_df.loc[:, 'tf_clean'] = filtered_df['tf'].str.replace('_TF', '', regex=False)


    return filtered_df


def extract_receptor_tf_pairs(
    receptor_gene_layer,
    top_tfs,
    top_receptors
) -> pd.DataFrame:

    """
    Given a receptor_gene_layer DataFrame (from get_celltype_grn_receptor_bip
    artite or get_celltype_gene_layer) and two “top‐N” DataFrames
    (one for tfs, one for receptors, each coming from results_df), return only
    those rows of receptor_gene_layer whose tf and receptor both appear in the
    respective top‐N lists.
    This function assumes:
        • receptor_gene_layer has columns "col1" (receptor) and "col2
            (tf).
        • top_tfs has a column "node" containing tf names in the form
            "GENE_TF::CellType".
        • top_receptors has a column "node" containing receptor names in
            either "GENE_receptor::CellType" or "GENE-receptor-CellType" or
            "GENE::CellType" form.
    We therefore standardize all three sources to "GENE::CellType" before filtering.

    Parameters
    ----------
    receptor_gene_layer : pandas.DataFrame
        Output of get_celltype_grn_receptor_bipartite() or
        get_celltype_gene_layer(), with at least:
            • "col1" (receptor) e.g. "GENE_receptor::Endothelial cell"
            • "col2" (tf) e.g. "GENE_TF::Endothelial cell"
    top_tfs : pandas.DataFrame
        Subset of results_df for multiplex=="{cell_type}_grn", with column:
            • "node" e.g. "GENE_TF::Endothelial cell"
    top_receptors : pandas.DataFrame
        Subset of results_df for multiplex=="cell_communication", with column:
            • "node" e.g. "GENE_receptor::Endothelial cell" or "GENE-receptor-Endothelial cell"
            or "GENE::Endothelial cell"

    Returns
    -------
    pandas.DataFrame
        A filtered DataFrame containing only those rows of receptor_gene_layer where:
            tf ∈ standardized(top_tfs["node"])  AND
            receptor ∈ standardized(top_receptors["node"]).
    """

    sources_list = list(top_tfs["node"].values)
    sources_list = ["".join(e.split("_TF")) for e in sources_list]

    targets_list = list(top_receptors["node"].values)

    print(sources_list)
    print(targets_list)
    # (2) Now filter tf_gene_layer for rows where both ‘source’ ∈ sources_list AND ‘target’ ∈ targets_list:
    filtered_df = receptor_gene_layer[
        receptor_gene_layer.loc[:, "col1"].isin(sources_list) &
        receptor_gene_layer.loc[:, "col2"].isin(targets_list)
    ]

    filtered_df.columns = ["receptor", "tf", "weight", "network_key"]
    return filtered_df


def extract_receptor_ligand_pairs(
    receptor_ligand_df: pd.DataFrame,
    top_ligands_df: pd.DataFrame,
    top_receptors_df: pd.DataFrame
) -> pd.DataFrame:
    """
    Given a receptor–ligand DataFrame (ligand_receptor_df) and two “top‐N” DataFrames
    (one for ligands, one for receptors, each coming from results_df), return only
    those rows of receptor_ligand_df whose ligand and receptor both appear in the
    respective top‐N lists.

    This function assumes:
      • receptor_ligand_
      df has columns “ligand” and “receptor” in the form “GENE::CellType”.
      • top_ligands_df has a column “node” containing ligand names in either
        “GENE-CellType” or “GENE::CellType” form.
      • top_receptors_df has a column “node” containing receptor names in either
        “GENE_receptor::CellType” or “GENE-receptor-CellType” or “GENE::CellType” form.

    We therefore standardize all three sources to “GENE::CellType” before filtering.

    Parameters
    ----------
    receptor_ligand_df : pandas.DataFrame
        Output of extract_receptor_ligand_pairs(), with at least:
          • "ligand"   e.g. "VEGFA::Endothelial cell"
          • "receptor" e.g. "FLT1::Fibroblast"
    top_ligands_df : pandas.DataFrame
        Subset of results_df for multiplex=="cell_communication", with column:
          • "node"   e.g. "VEGFA-Endothelial cell"  or  "VEGFA::Endothelial cell"
    top_receptors_df : pandas.DataFrame
        Subset of results_df for, e.g., "Fibroblast_receptor", with column:
          • "node"   e.g. "FLT1_receptor::Fibroblast" or "FLT1-Fibroblast" or "FLT1::Fibroblast"

    Returns
    -------
    pandas.DataFrame
        A filtered DataFrame containing only those rows of receptor_ligand_df where:
          ligand ∈ standardized(top_ligands_df["node"])  AND
          receptor ∈ standardized(top_receptors_df["node"]).
    """

    # Helper 1: standardize any "GENE-CellType" → "GENE::CellType"
    def _hyphen_to_doublecolon(x: str) -> str:
        if "-" in x and "::" not in x:
            gene_part, cell_part = x.rsplit("-", 1)
            return f"{gene_part}::{cell_part}"
        else:
            return x

    # Helper 2: strip off any "_receptor" suffix and standardize hyphens
    def _std_receptor(x: str) -> str:
        # Case A: "GENE_receptor::CellType"
        if "_receptor::" in x:
            gene, cell = x.split("_receptor::", 1)
            return f"{gene}::{cell}"
        # Case B: "GENE-receptor-CellType" (rare)
        if "-receptor-" in x and "::" not in x:
            gene_part, cell_part = x.rsplit("-receptor-", 1)
            return f"{gene_part}::{cell_part}"
        # Otherwise, just apply the hyphen-to-doublecolon rule
        return _hyphen_to_doublecolon(x)

    # 1) Build set of standardized ligand‐names from top_ligands_df["node"]
    raw_lig_nodes = top_ligands_df.loc[:, "node"].astype(str)
    standardized_ligs = set(raw_lig_nodes.apply(_hyphen_to_doublecolon))

    # 2) Build set of standardized receptor‐names from top_receptors_df["node"]
    raw_rec_nodes = top_receptors_df.loc[:, "node"].astype(str)
    standardized_recs = set(raw_rec_nodes.apply(_std_receptor))

    # 3) Filter receptor_ligand_df:
    #    Keep only rows where ligand ∈ standardized_ligs AND receptor ∈ standardized_recs
    mask = (
        receptor_ligand_df.loc[:, "ligand"].isin(standardized_ligs) &
        receptor_ligand_df.loc[:, "receptor"].isin(standardized_recs)
    )
    
    return receptor_ligand_df.loc[mask].reset_index(drop=True)


def build_partial_networks(
    multicell_obj,
    results,
    cell_type: str,
    seeds: Union[List[str], pd.Series],
    ligand_cells: List[str] = None,
    top_ligand_n: int = 100,
    top_receptor_n: int = 30,
    top_tf_n: int = 10,
    before_top_n: int = 5,
    per_celltype: bool = True,
    include_before_cells: bool = False
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """
    Construct partial or full network layers needed for Sankey plots.

    Parameters
    ----------
    include_before_cells : bool
        If True, include the before-cell layers (receptor → TF and TF → ligand).
        If False, returns empty DataFrames for those layers.

    Returns
    -------
    Tuple of 5 DataFrames:
      (before_receptor_tf_df, before_tf_ligand_df,
       receptor_ligand_df, receptor_tf_df, gene_tf_df)
    """
    cc_df = get_cell_communication_layer(
        multicell_obj,
        as_dataframe=True,
        ligand_cells=ligand_cells,
        receptor_cells=[cell_type])

    top_ligands = get_top_ligands(results, cc_df, n=top_ligand_n, per_celltype=per_celltype)
    top_receptors = get_top_receptors(results, cell_type=cell_type, n=top_receptor_n)
    top_tfs = get_top_tfs(results, cell_type=cell_type, n=top_tf_n)

    receptor_ligand_top = extract_receptor_ligand_pairs(
        receptor_ligand_df=cc_df,
        top_ligands_df=top_ligands,
        top_receptors_df=top_receptors
    )

    tf_gene_df = get_celltype_gene_layer(
        multicell_obj=multicell_obj,
        cell_type=cell_type,
        layer_name="grn",
        as_dataframe=True
    )

    receptor_tf_df = get_celltype_grn_receptor_bipartite(
        multicell_obj=multicell_obj,
        cell_type=cell_type,
        as_dataframe=True
    )

    seeds_prefixed = pd.Series([f"{gene}::{cell_type}" for gene in seeds])
    gene_tf_pairs = extract_gene_tf_pairs(tf_gene_df, top_tfs, seeds_prefixed)
    receptor_tf_pairs = extract_receptor_tf_pairs(receptor_tf_df, top_tfs, top_receptors)

    if not include_before_cells:
        # Return empty frames for before-layers
        return (
            pd.DataFrame(columns=["receptor", "tf", "weight"]),
            pd.DataFrame(columns=["tf_clean", "gene", "weight"]),
            receptor_ligand_top,
            receptor_tf_pairs,
            gene_tf_pairs
        )

    # Otherwise build full before-cell layers
    before_cell_types = top_ligands["ligand_celltype"].unique()
    all_before_receptor_tf_pairs = []
    all_before_gene_tf_pairs = []

    for before_cell_type in before_cell_types:
        before_top_receptors = get_top_receptors(results, cell_type=before_cell_type, n=before_top_n)
        before_top_tfs = get_top_tfs(results, cell_type=before_cell_type, n=before_top_n)

        before_tf_gene_df = get_celltype_gene_layer(multicell_obj, before_cell_type, "grn", as_dataframe=True)
        before_receptor_tf_df = get_celltype_grn_receptor_bipartite(multicell_obj, before_cell_type, as_dataframe=True)

        before_gene_tf_pairs = extract_gene_tf_pairs(
            tf_gene_layer=before_tf_gene_df,
            top_tfs=before_top_tfs,
            seeds=receptor_ligand_top[
                receptor_ligand_top["celltype_source"] == before_cell_type
            ]["ligand"].values
        )

        before_receptor_tf_pairs = extract_receptor_tf_pairs(
            receptor_gene_layer=before_receptor_tf_df,
            top_tfs=before_top_tfs,
            top_receptors=before_top_receptors
        )

        all_before_receptor_tf_pairs.append(before_receptor_tf_pairs)
        all_before_gene_tf_pairs.append(before_gene_tf_pairs)

    all_before_receptor_tf_df = pd.concat(all_before_receptor_tf_pairs, ignore_index=True)
    all_before_gene_tf_df = pd.concat(all_before_gene_tf_pairs, ignore_index=True)

    return (
        all_before_receptor_tf_df,
        all_before_gene_tf_df,
        receptor_ligand_top,
        receptor_tf_pairs,
        gene_tf_pairs
    )


def plot_3layer_sankey(
    receptor_tf_df: pd.DataFrame,
    gene_tf_df: pd.DataFrame,
    cell_type: Union[str, None] = None,
    flow: str = "upstream",
    color: str = "rgba(160, 160, 160, 0.4)",
    save_path=None
):
    def format_links(df, source_col, target_col):
        return df.loc[:, [source_col, target_col, "weight"]].rename(columns={
            source_col: "source",
            target_col: "target",
            "weight": "value"
        })

    r_t = format_links(receptor_tf_df, "receptor", "tf")
    t_g = format_links(gene_tf_df, "tf_clean", "gene")

    r_t = r_t[r_t["target"].isin(t_g["source"])]
    t_g = t_g[t_g["source"].isin(r_t["target"])]

    if flow.lower() == "downstream":
        for df in [r_t, t_g]:
            df[["source", "target"]] = df[["target", "source"]]

    r_t["color"] = color
    t_g["color"] = color

    for df in [r_t, t_g]:
        total = df["value"].sum()
        if total > 0:
            df["value"] /= total

    links = pd.concat([r_t, t_g], ignore_index=True)

    all_nodes = pd.unique(links[["source", "target"]].values.ravel())
    node_idx = {name: i for i, name in enumerate(all_nodes)}
    links["source_idx"] = links["source"].map(node_idx)
    links["target_idx"] = links["target"].map(node_idx)

    def _format_label(x: str) -> str:
        parts = x.split("::", 1)
        return parts[0] if len(parts) == 2 else x

    labels = [_format_label(n) for n in all_nodes]

    sankey_data = go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=labels
        ),
        link=dict(
            source=links["source_idx"],
            target=links["target_idx"],
            value=links["value"],
            color=links["color"]
        ),
        orientation="h"
    )

    title_text = (
        f"Top regulators in {cell_type}: Receptor → TF → Gene"
        if flow.lower() == "upstream"
        else f"Top regulators in {cell_type}: Gene → TF → Receptor"
    )
    layer_names = (
        ["Receptors", "TFs", "Genes"]
        if flow.lower() == "upstream"
        else ["Genes", "TFs", "Receptors"]
    )
    x_positions = [0.0, 0.5, 1.0]

    fig = go.Figure(data=[sankey_data])

    fig.update_layout(title_text=title_text, font_size=14, font_color="black")

    for x, name in zip(x_positions, layer_names):
        fig.add_annotation(
            x=x, y=-0.15,
            text=f"<b>{name}</b>",
            showarrow=False,
            font=dict(size=16)
        )

    if save_path:
        fig.write_html(save_path)

    fig.show(renderer="notebook")


def plot_4layer_sankey(
    ligand_receptor_df: pd.DataFrame,
    receptor_tf_df: pd.DataFrame,
    gene_tf_df: pd.DataFrame,
    flow: str = "upstream",
    save_path: Union[str, None] = None
):
    """
    Plot a 4-layer Sankey diagram showing ligand → receptor → TF → gene
    or the reverse flow.
    
    Parameters
    ----------
    ligand_receptor_df : pandas.DataFrame
        Output of extract_receptor_ligand_pairs(), with at least:
          - "ligand"   (e.g. "VEGFA::Endothelial cell")
          - "receptor" (e.g. "FLT1::Fibroblast")
          - "receptor_clean" (e.g. "FLT1_receptor::Fibroblast")
          - "weight"   (numeric)
    receptor_tf_df : pandas.DataFrame
        Output of extract_receptor_tf_pairs(), with at least:
          - "receptor" (e.g. "FLT1::Fibroblast")
          - "tf"       (e.g. "STAT3_TF::Fibroblast")
          - "weight"   (numeric)
    gene_tf_df : pandas.DataFrame
        Output of extract_gene_tf_pairs(), with at least:
          - "tf_clean" (e.g. "STAT3::Fibroblast")
          - "gene"     (e.g. "JUN::Fibroblast")
          - "weight"   (numeric)
    flow : str, default="upstream"
        If "upstream", plot ligand → receptor → TF → gene.
        If "downstream", plot gene → TF → receptor → ligand.
    save_path : str or None, default=None
        If provided, save the figure as an HTML file to this path.

    Returns
    -------
    None
    """
    def format_links(df, source_col, target_col):
        return df.loc[:, [source_col, target_col, "weight"]].rename(columns={
            source_col: "source",
            target_col: "target",
            "weight": "value"
        })

    def hex_to_rgba(hex_color, alpha=0.6):
        h = hex_color.lstrip("#")
        return f"rgba({int(h[0:2], 16)}, {int(h[2:4], 16)}, {int(h[4:6], 16)}, {alpha})"

    def string_to_color(string):
        h = hashlib.md5(string.encode()).hexdigest()
        return "#" + h[:6]

    def assign_group_colors(df, column):
        unique_types = df[column].str.extract(r"::(.+)$")[0].fillna("Unknown")
        return unique_types.apply(lambda ct: hex_to_rgba(string_to_color(ct)))

    l_r = format_links(ligand_receptor_df, "ligand", "receptor_clean")
    r_t = format_links(receptor_tf_df, "receptor", "tf")
    t_g = format_links(gene_tf_df, "tf_clean", "gene")

    l_r = l_r[l_r["target"].isin(r_t["source"])]
    r_t = r_t[r_t["source"].isin(l_r["target"])]
    t_g = t_g[t_g["source"].isin(r_t["target"])]

    if flow.lower() == "downstream":
        for df in [l_r, r_t, t_g]:
            df[["source", "target"]] = df[["target", "source"]]

    l_r["color"] = assign_group_colors(l_r, "source")
    r_t["color"] = "rgba(100,200,100,1)"
    t_g["color"] = "rgba(100,200,100,1)"

    for df in [l_r, r_t, t_g]:
        total = df["value"].sum()
        if total > 0:
            df["value"] /= total

    links = pd.concat([l_r, r_t, t_g], ignore_index=True)

    all_nodes = pd.unique(links[["source", "target"]].values.ravel())
    node_idx = {name: i for i, name in enumerate(all_nodes)}
    links["source_idx"] = links["source"].map(node_idx)
    links["target_idx"] = links["target"].map(node_idx)

    def _format_label(x: str) -> str:
        parts = x.split("::", 1)
        return parts[0] if len(parts) == 2 else x

    labels = [_format_label(n) for n in all_nodes]

    sankey_data = go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=labels
        ),
        link=dict(
            source=links["source_idx"],
            target=links["target_idx"],
            value=links["value"],
            color=links["color"]
        ),
        orientation="h"
    )

    title_text = (
        #l_r["source"].str.extract(r"::(.+)$")[0].unique()
        f"Top regulators from other cell types:\n Ligand → Receptor → TF → Gene"
        if flow.lower() == "upstream"
        else f"Top regulators from other cell types:\n Gene → TF → Receptor → Ligand"
    )
    layer_names = (
        ["Ligands", "Receptors", "TFs", "Genes"]
        if flow.lower() == "upstream"
        else ["Genes", "TFs", "Receptors", "Ligands"]
    )
    x_positions = [0.0, 0.33, 0.66, 1.0]

    fig = go.Figure(data=[sankey_data])
    fig.update_layout(title_text=title_text, font_size=14, font_color="black")

    color_map = (
        pd.concat([l_r, r_t, t_g])[["source", "color"]]
        .dropna()
        .drop_duplicates()
    )
    color_map["celltype"] = color_map["source"].str.extract(r"::(.+)$")[0]
    color_map = color_map.dropna(subset=["celltype"])

    # Final dictionary: celltype → color
    color_map = dict(zip(color_map["celltype"], color_map["color"]))

    for i, (ct, color) in enumerate(sorted(color_map.items())):
        fig.add_annotation(
            x=1.02, y=1.0 - i * 0.05,
            text=f"<b>{ct}</b>",
            showarrow=False,
            font=dict(size=12),
            bgcolor=color,
            bordercolor="black",
            borderwidth=0.5,
            align="left",
            xanchor="left"
        )

    for x, name in zip(x_positions, layer_names):
        fig.add_annotation(
            x=x, y=-0.15,
            text=f"<b>{name}</b>",
            showarrow=False,
            font=dict(size=16)
        )

    if save_path:
        fig.write_html(save_path)

    fig.show(renderer="notebook")


def plot_6layer_sankey(
    before_receptor_tf_df: pd.DataFrame,
    before_tf_ligand_df: pd.DataFrame,
    ligand_receptor_df: pd.DataFrame,
    receptor_tf_df: pd.DataFrame,
    gene_tf_df: pd.DataFrame,
    flow: str = "upstream",
    save_path: Union[str, None] = None
):
    def format_links(df, source_col, target_col):
        return df.loc[:, [source_col, target_col, "weight"]].rename(columns={
            source_col: "source",
            target_col: "target",
            "weight": "value"
        })

    def hex_to_rgba(hex_color, alpha=0.6):
        h = hex_color.lstrip("#")
        return f"rgba({int(h[0:2], 16)}, {int(h[2:4], 16)}, {int(h[4:6], 16)}, {alpha})"

    def string_to_color(string):
        h = hashlib.md5(string.encode()).hexdigest()
        return "#" + h[:6]

    # Format all links
    br_bt = format_links(before_receptor_tf_df, "receptor", "tf")
    bt_l = format_links(before_tf_ligand_df, "tf_clean", "gene")
    l_r = format_links(ligand_receptor_df, "ligand", "receptor_clean")
    r_t = format_links(receptor_tf_df, "receptor", "tf")
    t_g = format_links(gene_tf_df, "tf_clean", "gene")

    # Filter unconnected
    l_r = l_r[l_r["target"].isin(r_t["source"])]
    l_r = l_r[l_r["source"].isin(bt_l["target"])]
    r_t = r_t[r_t["source"].isin(l_r["target"])]
    t_g = t_g[t_g["source"].isin(r_t["target"])]
    bt_l = bt_l[bt_l["target"].isin(l_r["source"])]
    bt_l = bt_l[bt_l["source"].isin(br_bt["target"])]
    l_r = l_r[l_r["source"].isin(bt_l["target"])]
    br_bt = br_bt[br_bt["target"].isin(bt_l["source"])]

    if flow.lower() == "downstream":
        for df in [br_bt, bt_l, l_r, r_t, t_g]:
            df.loc[:, ["source", "target"]] = df.loc[:, ["target", "source"]]

    # Assign colors by celltype
    def assign_group_colors(df, column):
        unique_types = df[column].str.extract(r"::(.+)$")[0].fillna("Unknown")
        return unique_types.apply(lambda ct: hex_to_rgba(string_to_color(ct)))

    br_bt.loc[:, "color"] = assign_group_colors(br_bt, "source")
    bt_l.loc[:, "color"] = assign_group_colors(bt_l, "source")
    l_r.loc[:, "color"] = "rgba(160,160,160,0.4)"
    r_t.loc[:, "color"] = "rgba(100,200,100,1)"
    t_g.loc[:, "color"] = "rgba(100,200,100,1)"

    for df in [br_bt, bt_l, l_r, r_t, t_g]:
        total = df["value"].sum()
        if total > 0:
            df["value"] /= total

    links = pd.concat([br_bt, bt_l, l_r, r_t, t_g], ignore_index=True)

    all_nodes = pd.unique(links[["source", "target"]].values.ravel())
    node_idx = {name: i for i, name in enumerate(all_nodes)}
    links.loc[:, "source_idx"] = links.loc[:, "source"].map(node_idx)
    links.loc[:, "target_idx"] = links.loc[:, "target"].map(node_idx)

    def _format_label(x: str) -> str:
        parts = x.split("::", 1)
        return parts[0].split("_")[0] if len(parts) == 2 else x

    labels = [_format_label(n) for n in all_nodes]

    sankey_data = go.Sankey(
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=labels,
        ),
        link=dict(
            source=links.loc[:, "source_idx"],
            target=links.loc[:, "target_idx"],
            value=links.loc[:, "value"],
            color=links.loc[:, "color"]
        ),
        orientation="h"
    )

    fig = go.Figure(data=[sankey_data])

    color_map = (
        pd.concat([br_bt, bt_l, r_t, t_g]).loc[:, ["source", "color"]]
        .dropna()
        .drop_duplicates()
    )

    color_map["celltype"] = color_map.loc[:, "source"].str.extract(r"::(.+)$")[0]
    color_map = color_map.dropna(subset=["celltype"])

    # Final dictionary: celltype → color
    color_map = dict(zip(color_map.loc[:, "celltype"], color_map.loc[:, "color"]))

    for i, (ct, color) in enumerate(sorted(color_map.items())):
        fig.add_annotation(
            x=1.02, y=1.0 - i * 0.05,
            text=f"<b>{ct}</b>",
            showarrow=False,
            font=dict(size=12),
            bgcolor=color,
            bordercolor="black",
            borderwidth=0.5,
            align="left",
            xanchor="left"
        )

    if flow.lower() == "upstream":
        title_text = "Upstream Receptor → Upstream TF → Ligand → Downstream Receptor → Downstream TF → Gene"
        layer_names = ["Upstream Receptors", "Upstream TFs", "Ligands", "Receptors", "TFs", "Genes"]
    else:
        title_text = "Gene → First TF → First Receptor → Ligand → Upstream TF → Upstream Receptor"
        layer_names = ["Genes", "TFs", "Receptors", "Ligands", "Upstream TFs", "Upstream Receptors"]

    fig.update_layout(title_text=title_text, font_size=14, font_color="black")

    x_positions = [0.0, 0.2, 0.4, 0.6, 0.8, 1.0]

    for x, name in zip(x_positions, layer_names):
        fig.add_annotation(
            x=x, y=-0.15,
            text=f"<b>{name}</b>",
            showarrow=False,
            font=dict(size=16)
        )

    if save_path:
        fig.write_html(save_path)

    fig.show(renderer="notebook")


def plot_intracell_sankey(
    multicell_obj,
    results,
    cell_type,
    seeds,
    top_receptor_n: int = 30,
    top_tf_n: int = 10,
    flow="upstream",
    save_path=None
):
    """
    Plot the intracellular regulatory Sankey diagram for a given cell type.
    This includes only the receptor → TF → gene layers, without any ligand or
    before-cell layers.

    Parameters
    ----------
    multicell_obj : object
        The multicell object containing the data.
    results : pandas.DataFrame
        The results DataFrame containing scores for nodes.
    cell_type : str
        The cell type of interest.
    seeds : list or pd.Series
        List of gene names (in the form "GENE::CellType") to filter for.
    top_receptor_n : int, default=30
        Number of top receptors to include.
    top_tf_n : int, default=10
        Number of top transcription factors to include.
    flow : str, default="upstream"
        If "upstream", plot receptor → TF → gene.
        If "downstream", plot gene → TF → receptor.
    save_path : str or None, default=None
        If provided, save the figure as an HTML file to this path.

    Returns
    -------
    None
    """
    networks = build_partial_networks(
        multicell_obj=multicell_obj,
        results=results,
        cell_type=cell_type,
        seeds=seeds,
        ligand_cells=[],
        top_receptor_n=top_receptor_n,
        top_tf_n=top_tf_n,
        include_before_cells=False
    )

    networks = {
        "before_receptor_tf": networks[0],
        "before_tf_ligand": networks[1],
        "ligand_receptor": networks[2],
        "receptor_tf": networks[3],
        "gene_tf": networks[4],
    }

    plot_3layer_sankey(
        receptor_tf_df=networks["receptor_tf"],
        gene_tf_df=networks["gene_tf"],
        cell_type=cell_type,
        flow=flow,
        save_path=save_path
    )


def plot_ligand_sankey(
    multicell_obj,
    results,
    cell_type,
    seeds,
    ligand_cells,
    top_ligand_n: int = 100,
    top_receptor_n: int = 30,
    top_tf_n: int = 10,
    per_celltype: bool = True,
    flow="upstream",
    save_path=None
):
    """
    Plot the ligand regulatory Sankey diagram for a given cell type.
    This includes the ligand → receptor → TF → gene layers, without any
    before-cell layers.

    Parameters
    ----------
    multicell_obj : object
        The multicell object containing the data.
    results : pandas.DataFrame
        The results DataFrame containing scores for nodes.
    cell_type : str
        The cell type of interest.
    seeds : list or pd.Series
        List of gene names (in the form "GENE::CellType") to filter for.
    ligand_cells : list of str
        List of cell types to consider as ligand sources.
    top_ligand_n : int, default=100
        Number of top ligands to include.
    top_receptor_n : int, default=30
        Number of top receptors to include.
    top_tf_n : int, default=10
        Number of top transcription factors to include.
    per_celltype : bool, default=True
        If True, select top ligands per ligand cell type.
        If False, select top ligands globally.
    flow : str, default="upstream"
        If "upstream", plot ligand → receptor → TF → gene.
        If "downstream", plot gene → TF → receptor → ligand.
    save_path : str or None, default=None
        If provided, save the figure as an HTML file to this path.

    Returns
    -------
    None
    """
    networks = build_partial_networks(
        multicell_obj=multicell_obj,
        results=results,
        cell_type=cell_type,
        seeds=seeds,
        ligand_cells=ligand_cells,
        top_ligand_n=top_ligand_n,
        top_receptor_n=top_receptor_n,
        top_tf_n=top_tf_n,
        per_celltype=per_celltype,
        include_before_cells=False
    )

    networks = {
        "before_receptor_tf": networks[0],
        "before_tf_ligand": networks[1],
        "ligand_receptor": networks[2],
        "receptor_tf": networks[3],
        "gene_tf": networks[4],
    }

    ligand_receptor_df = networks["ligand_receptor"]
    receptor_tf_df = networks["receptor_tf"]
    gene_tf_df = networks["gene_tf"]

    plot_4layer_sankey(
        ligand_receptor_df=ligand_receptor_df,
        receptor_tf_df=receptor_tf_df,
        gene_tf_df=gene_tf_df,
        flow=flow,
        save_path=save_path
    )


def plot_intercell_sankey(
    multicell_obj,
    results,
    cell_type,
    seeds,
    ligand_cells,
    top_ligand_n: int = 100,
    top_receptor_n: int = 30,
    top_tf_n: int = 10,
    before_top_n: int = 5,
    per_celltype: bool = True,
    flow="upstream",
    save_path=None
):
    """
    Plot the full intercellular regulatory Sankey diagram for a given cell type.
    This includes the before-cell layers (receptor → TF and TF → ligand) as well
    as the ligand → receptor → TF → gene layers.

    Parameters
    ----------
    multicell_obj : object
        The multicell object containing the data.
    results : pandas.DataFrame
        The results DataFrame containing scores for nodes.
    cell_type : str
        The cell type of interest.
    seeds : list or pd.Series
        List of gene names (in the form "GENE::CellType") to filter for.
    ligand_cells : list of str
        List of cell types to consider as ligand sources.
    top_ligand_n : int, default=100
        Number of top ligands to include.
    top_receptor_n : int, default=30
        Number of top receptors to include.
    top_tf_n : int, default=10
        Number of top transcription factors to include.
    before_top_n : int, default=5
        Number of top receptors and TFs to include in the before-cell layers.
    per_celltype : bool, default=True
        If True, select top ligands per ligand cell type.
        If False, select top ligands globally.
    flow : str, default="upstream"
        If "upstream", plot upstream receptor → upstream TF → ligand → receptor → TF → gene.
        If "downstream", plot gene → TF → receptor → ligand → upstream TF → upstream receptor.
    save_path : str or None, default=None
        If provided, save the figure as an HTML file to this path.

    Returns
    -------
    None
    """

    networks = build_partial_networks(
        multicell_obj=multicell_obj,
        results=results,
        cell_type=cell_type,
        seeds=seeds,
        ligand_cells=ligand_cells,
        top_ligand_n=top_ligand_n,
        top_receptor_n=top_receptor_n,
        top_tf_n=top_tf_n,
        before_top_n=before_top_n,
        per_celltype=per_celltype,
        include_before_cells=True
    )

    networks = {
        "before_receptor_tf": networks[0],
        "before_tf_ligand": networks[1],
        "ligand_receptor": networks[2],
        "receptor_tf": networks[3],
        "gene_tf": networks[4],
    }

    plot_6layer_sankey(
        before_receptor_tf_df=networks["before_receptor_tf"],
        before_tf_ligand_df=networks["before_tf_ligand"],
        ligand_receptor_df=networks["ligand_receptor"],
        receptor_tf_df=networks["receptor_tf"],
        gene_tf_df=networks["gene_tf"],
        flow=flow,
        save_path=save_path
    )
