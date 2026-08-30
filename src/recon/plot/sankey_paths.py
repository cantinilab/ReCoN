import pandas as pd
import networkx as nx
from typing import Union, List, Tuple
import plotly.graph_objects as go
import hashlib


def _normalize_seed_nodes(seeds: Union[List[str], pd.Series], cell_type: str) -> pd.Series:
    """Return seed genes in the canonical ``GENE::CellType`` form."""
    seed_values = seeds.tolist() if hasattr(seeds, "tolist") else list(seeds)
    normalized = []

    for seed in seed_values:
        seed = str(seed)
        if "::" in seed:
            normalized.append(seed)
        else:
            normalized.append(f"{seed}::{cell_type}")

    return pd.Series(normalized)


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
    layer_name: str = "gene",
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

    # Step 1: build a dict mapping each name to its DataFrame.
    name_to_layer = dict(zip(mux_info["names"], mux_info["layers"]))

    # ``gene`` is the semantic intracellular GRN layer expected by the cascade
    # plotting code. Some older ReCoN objects store the same layer as ``grn``.
    layer_lookup_order = [layer_name]
    if layer_name == "gene":
        layer_lookup_order.extend(["grn", "tf_gene"])

    layer = None
    for candidate in layer_lookup_order:
        if candidate in name_to_layer:
            layer = name_to_layer[candidate]
            break
    if layer is None:
        available = list(name_to_layer)
        raise KeyError(
            f"Multiplex '{multiplex_key}' has no layer {layer_name} listed. "
            f"Available layers: {available}"
        )

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

    # 2) By convention:  source_std → ligand (sender),  target_std → receptor (receiver)
    df.loc[:, "ligand"] = df.loc[:, "source_std"]
    df.loc[:, "receptor"] = df.loc[:, "target_std"]

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

    results = results_df.sort_values(by="score", ascending=False)
    if results.empty:
        return results
    # Strip the known suffix instead of splitting on "_"
    results.loc[:, "celltype"] = results["multiplex"].str.replace(r"_(grn|receptor)$", "", regex=True)
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

    results = results_df.sort_values(by="score", ascending=False)
    if results.empty:
        return results
    # Strip the known suffix instead of splitting on "_"
    results.loc[:, "celltype"] = results["multiplex"].str.replace(r"_(grn|receptor)$", "", regex=True)
    results = results[results["celltype"]==cell_type]
    
    # remove fake node
    results = results[results["node"]!=f"fake_receptor::{cell_type}"]

    results = results[results["node"].str.endswith(f"_receptor::{cell_type}")]
    return results.iloc[:n, :]


def _seed_names_for_cell_type(seeds, cell_type: str) -> pd.Series:
    """Return seed names in GENE::CellType form."""
    seed_index = pd.Index(seeds.keys() if isinstance(seeds, dict) else seeds)
    seed_values = seed_index.astype(str)
    return pd.Series([
        seed if "::" in seed else f"{seed}::{cell_type}"
        for seed in seed_values
    ])


def get_top_genes(
    results_df: pd.DataFrame,
    cell_type,
    n: int = 5,
    exclude_nodes=None,
) -> pd.DataFrame:
    """
    From a `results_df` of node-scores, identify top gene nodes for a cell type.

    Gene nodes are plain ``GENE::CellType`` nodes in the GRN multiplex; TF and
    receptor helper nodes keep their ``_TF``/``_receptor`` layer suffixes.
    """
    results = results_df.sort_values(by="score", ascending=False)
    if results.empty:
        return results

    results = results.copy()
    results.loc[:, "celltype"] = results["multiplex"].str.replace(r"_(grn|receptor)$", "", regex=True)
    results = results[results["celltype"] == cell_type]
    results = results[results["node"].str.endswith(f"::{cell_type}")]
    results = results[
        ~results["node"].str.contains(r"_(?:TF|receptor)::", regex=True)
    ]

    if exclude_nodes is not None:
        results = results[~results["node"].isin(set(exclude_nodes))]

    return results.iloc[:n, :]


def _normalize_flow(flow: str) -> str:
    flow = flow.lower()
    if flow not in {"upstream", "downstream"}:
        raise ValueError("flow must be 'upstream' or 'downstream'")
    return flow


def _top_seed_connected_tfs(results, tf_gene_layer, seeds, cell_type, n):
    seeds = set(pd.Index(seeds).astype(str))
    source = tf_gene_layer["source"].astype(str)
    target = tf_gene_layer["target"].astype(str)

    # Support both TF -> gene and older gene -> TF GRN orientations.
    connected_tfs = set(source[target.isin(seeds)])
    connected_tfs.update(target[source.isin(seeds) & _is_tf_node_series(target)])

    top_tfs = get_top_tfs(results, cell_type=cell_type, n=len(results))
    top_tfs = top_tfs[top_tfs["node"].isin(connected_tfs)]
    return top_tfs.iloc[:n, :]


def _top_tf_connected_receptors(
    results, receptor_tf_layer, top_tfs, cell_type, n, receptor_candidates=None
):
    tf_nodes = set(top_tfs["node"]) if "node" in top_tfs.columns else set()
    tf_nodes_clean = {tf.replace("_TF::", "::") for tf in tf_nodes}
    connected_receptors = set(receptor_tf_layer.loc[
        receptor_tf_layer["col1"].isin(tf_nodes) |
        receptor_tf_layer["col1"].isin(tf_nodes_clean),
        "col2"
    ])
    if receptor_candidates is not None:
        connected_receptors &= set(receptor_candidates)
    top_receptors = get_top_receptors(results, cell_type=cell_type, n=len(results))
    top_receptors = top_receptors[top_receptors["node"].isin(connected_receptors)]
    return top_receptors.iloc[:n, :]


def _ligand_source_celltypes(top_ligands, receiver_cell_type):
    return [
        ct for ct in top_ligands["ligand_celltype"].unique()
        if ct != receiver_cell_type
    ]


def _filter_connected_sankey_layers(br_bt, bt_l, l_r, r_t, t_g):
    """Keep only links connected to adjacent layers in the cascade."""
    layers = [br_bt, bt_l, l_r, r_t, t_g]

    def _connect(left, right):
        if len(left) == 0 or len(right) == 0:
            return left, right
        left = left[left["target"].isin(right["source"])]
        right = right[right["source"].isin(left["target"])]
        return left, right

    changed = True
    while changed:
        sizes_before = [len(layer) for layer in layers]
        layers[0], layers[1] = _connect(layers[0], layers[1])
        layers[1], layers[2] = _connect(layers[1], layers[2])
        layers[2], layers[3] = _connect(layers[2], layers[3])
        layers[3], layers[4] = _connect(layers[3], layers[4])

        has_left = [
            False,
            len(layers[0]) > 0 and len(layers[1]) > 0,
            len(layers[1]) > 0 and len(layers[2]) > 0,
            len(layers[2]) > 0 and len(layers[3]) > 0,
            len(layers[3]) > 0 and len(layers[4]) > 0,
        ]
        has_right = [
            len(layers[0]) > 0 and len(layers[1]) > 0,
            len(layers[1]) > 0 and len(layers[2]) > 0,
            len(layers[2]) > 0 and len(layers[3]) > 0,
            len(layers[3]) > 0 and len(layers[4]) > 0,
            False,
        ]

        for idx, layer in enumerate(layers):
            if len(layer) > 0 and not has_left[idx] and not has_right[idx]:
                layers[idx] = layer.iloc[0:0].copy()

        changed = sizes_before != [len(layer) for layer in layers]

    return tuple(layers)


def _node_celltype(node: str) -> str:
    return node.rsplit("::", 1)[-1] if "::" in node else ""


def _as_tf_node(node: str) -> str:
    if "_TF::" in node:
        return node
    if "::" in node:
        gene, cell = node.split("::", 1)
        return f"{gene}_TF::{cell}"
    return node


def _as_receptor_node(node: str) -> str:
    if "_receptor::" in node:
        return node
    if "::" in node:
        gene, cell = node.split("::", 1)
        return f"{gene}_receptor::{cell}"
    return node


def _clean_tf_node(node: str) -> str:
    return str(node).replace("_TF::", "::")


def _rank_nodes(results, nodes, n):
    if len(nodes) == 0:
        return pd.DataFrame(columns=results.columns)
    ranked = results[results["node"].isin(set(nodes))].sort_values("score", ascending=False)
    return ranked.iloc[:n, :]


def _rank_downstream_receptors(results, receptors, n, per_celltype):
    ranked = _rank_nodes(results, receptors, len(results))
    if not per_celltype or ranked.empty:
        return ranked.iloc[:n, :]

    ranked = ranked.copy()
    ranked.loc[:, "celltype"] = ranked["node"].astype(str).apply(_node_celltype)
    return (
        ranked
        .sort_values(["celltype", "score"], ascending=[True, False])
        .groupby("celltype", group_keys=False)
        .head(n)
        .drop(columns=["celltype"])
        .reset_index(drop=True)
    )


def _standardize_node_name(node: str) -> str:
    node = str(node)
    if "-" in node and "::" not in node:
        gene, cell = node.rsplit("-", 1)
        return f"{gene}::{cell}"
    return node


def _rank_ligands_for_downstream(results, ligands, n, per_celltype):
    if len(ligands) == 0:
        return pd.DataFrame(columns=["multiplex", "node", "score", "node_std", "ligand_celltype"])

    ligand_set = set(ligands)
    ranked = results[results["multiplex"] == "cell_communication"].copy()
    if ranked.empty:
        ranked = pd.DataFrame(columns=list(results.columns) + ["node_std", "ligand_celltype"])
    else:
        ranked.loc[:, "node_std"] = ranked["node"].astype(str).apply(_standardize_node_name)
        ranked = ranked[ranked["node_std"].isin(ligand_set)].sort_values("score", ascending=False)

    seen = set(ranked["node_std"]) if "node_std" in ranked.columns else set()
    missing = [lig for lig in ligands if lig not in seen]
    if missing:
        ranked = pd.concat([
            ranked,
            pd.DataFrame({
                "multiplex": "cell_communication",
                "node": missing,
                "score": 0.0,
                "node_std": missing,
            }),
        ], ignore_index=True, sort=False)

    ranked.loc[:, "ligand_celltype"] = ranked["node_std"].str.split("::", n=1).str[-1]
    if per_celltype:
        return (
            ranked
            .sort_values(["ligand_celltype", "score"], ascending=[True, False])
            .groupby("ligand_celltype", group_keys=False)
            .head(n)
            .reset_index(drop=True)
        )
    return ranked.sort_values("score", ascending=False).head(n).reset_index(drop=True)


def _is_tf_node_series(values):
    return values.astype(str).str.contains(r"_TF::", regex=True)


def _tf_name_sets(top_tfs):
    if "node" not in top_tfs.columns:
        return set(), set()
    tf_nodes = set(top_tfs["node"].astype(str))
    tf_clean = {tf.replace("_TF::", "::") for tf in tf_nodes}
    return tf_nodes, tf_clean


def _filter_receptor_targets_to_tfs(receptor_tf_df, top_tfs):
    tf_nodes, tf_clean = _tf_name_sets(top_tfs)
    return receptor_tf_df[
        receptor_tf_df["col1"].isin(tf_nodes) |
        receptor_tf_df["col1"].isin(tf_clean)
    ].copy()


def _normalize_downstream_plot_links(l_r, r_t, t_g):
    l_r = l_r.copy()
    r_t = r_t.copy()
    t_g = t_g.copy()

    if len(l_r) > 0:
        l_r.loc[:, "target"] = l_r["target"].astype(str).apply(_as_receptor_node)
    if len(r_t) > 0:
        r_t.loc[:, "source"] = r_t["source"].astype(str).apply(_as_receptor_node)
        r_t.loc[:, "target"] = r_t["target"].astype(str).apply(_clean_tf_node)
    if len(t_g) > 0:
        t_g.loc[:, "source"] = t_g["source"].astype(str).apply(_clean_tf_node)

    return l_r, r_t, t_g


def _format_sankey_links(df, source_col, target_col):
    return df.loc[:, [source_col, target_col, "weight"]].rename(columns={
        source_col: "source",
        target_col: "target",
        "weight": "value"
    })


def _normalize_link_values(link_dfs, flow):
    for df in link_dfs:
        total = df["value"].sum()
        if total > 0:
            df["value"] /= total


def _build_layered_links(layer_edges, layer_x):
    links_list = []
    node_ids = []
    labels_by_id = {}
    x_by_id = {}
    layer_by_id = {}

    for link_df, source_layer, target_layer in layer_edges:
        if len(link_df) == 0:
            continue
        link_df = link_df.copy()
        link_df.loc[:, "source_id"] = source_layer + "::" + link_df["source"].astype(str)
        link_df.loc[:, "target_id"] = target_layer + "::" + link_df["target"].astype(str)
        links_list.append(link_df)

    links = pd.concat(links_list, ignore_index=True)

    for col in ["source_id", "target_id"]:
        for node_id in links[col].dropna().astype(str):
            if node_id in labels_by_id:
                continue
            layer, label = node_id.split("::", 1)
            node_ids.append(node_id)
            labels_by_id[node_id] = label
            x_by_id[node_id] = layer_x[layer]
            layer_by_id[node_id] = layer

    node_idx = {node_id: i for i, node_id in enumerate(node_ids)}
    links.loc[:, "source_idx"] = links["source_id"].map(node_idx)
    links.loc[:, "target_idx"] = links["target_id"].map(node_idx)

    y_by_id = {}
    for layer in layer_x:
        layer_nodes = [node_id for node_id in node_ids if layer_by_id[node_id] == layer]
        if len(layer_nodes) == 1:
            y_by_id[layer_nodes[0]] = 0.5
            continue
        for idx, node_id in enumerate(layer_nodes):
            y_by_id[node_id] = 0.02 + (0.96 * idx / max(len(layer_nodes) - 1, 1))

    return (
        links,
        [labels_by_id[node_id] for node_id in node_ids],
        [x_by_id[node_id] for node_id in node_ids],
        [y_by_id[node_id] for node_id in node_ids],
    )


def _prepare_6layer_links(
    before_receptor_tf_df,
    before_tf_ligand_df,
    ligand_receptor_df,
    receptor_tf_df,
    gene_tf_df,
    flow,
):
    br_bt = _format_sankey_links(before_receptor_tf_df, "receptor", "tf")
    bt_l = _format_sankey_links(before_tf_ligand_df, "tf_clean", "gene")
    l_r = _format_sankey_links(ligand_receptor_df, "ligand", "receptor_clean")
    receptor_tf_col = "tf_clean" if flow == "downstream" and "tf_clean" in receptor_tf_df.columns else "tf"
    r_t = _format_sankey_links(receptor_tf_df, "receptor", receptor_tf_col)
    t_g = _format_sankey_links(gene_tf_df, "tf_clean", "gene")

    if flow == "downstream":
        l_r, r_t, t_g = _normalize_downstream_plot_links(l_r, r_t, t_g)

    return _filter_connected_sankey_layers(br_bt, bt_l, l_r, r_t, t_g)


def _extract_downstream_gene_tf_pairs(tf_gene_layer, top_tfs, top_genes):
    tfs, tfs_clean = _tf_name_sets(top_tfs)
    genes = set(top_genes["node"]) if "node" in top_genes.columns else set()
    filtered_df = tf_gene_layer[
        (
            tf_gene_layer["source"].isin(tfs) |
            tf_gene_layer["source"].isin(tfs_clean)
        ) &
        tf_gene_layer["target"].isin(genes)
    ].copy()
    filtered_df = filtered_df.rename(columns={"source": "tf", "target": "gene"})
    filtered_df.loc[:, "tf_clean"] = filtered_df["tf"].str.replace("_TF", "", regex=False)
    return filtered_df


def _build_downstream_networks(
    multicell_obj,
    results,
    cell_type,
    seeds,
    ligand_cells,
    top_ligand_n,
    top_receptor_n,
    top_tf_n,
    per_celltype,
    include_before_cells,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    seeds_prefixed = _seed_names_for_cell_type(seeds, cell_type)
    seed_tfs = set(seeds_prefixed[seeds_prefixed.str.contains(r"_TF::", regex=True)])
    seed_receptors = set(seeds_prefixed[seeds_prefixed.str.contains(r"_receptor::", regex=True)])
    plain_gene_seeds = set(seeds_prefixed) - seed_tfs - seed_receptors

    source_tf_gene_df = get_celltype_gene_layer(
        multicell_obj=multicell_obj,
        cell_type=cell_type,
        layer_name="gene",
        as_dataframe=True,
    )
    source_receptor_tf_df = get_celltype_grn_receptor_bipartite(
        multicell_obj=multicell_obj,
        cell_type=cell_type,
        as_dataframe=True,
    )

    if seed_receptors:
        source_receptor_tf_pairs = source_receptor_tf_df[
            source_receptor_tf_df["col2"].isin(seed_receptors)
        ].copy()
        downstream_source_tfs = set(source_receptor_tf_pairs["col1"])
        source_receptor_tf_pairs = source_receptor_tf_pairs.rename(
            columns={"col2": "receptor", "col1": "tf"}
        )
    else:
        source_receptor_tf_pairs = pd.DataFrame(columns=["receptor", "tf", "weight"])
        downstream_source_tfs = set()

    source_tfs = seed_tfs | downstream_source_tfs
    source_tf_gene_pairs = source_tf_gene_df[
        source_tf_gene_df["source"].isin(source_tfs)
    ].copy()
    source_ligand_candidates = set(source_tf_gene_pairs["target"]) | plain_gene_seeds
    source_tf_gene_pairs = source_tf_gene_pairs.rename(columns={"source": "tf", "target": "gene"})
    if len(source_tf_gene_pairs):
        source_tf_gene_pairs.loc[:, "tf_clean"] = source_tf_gene_pairs["tf"].str.replace(
            "_TF", "", regex=False
        )
    else:
        source_tf_gene_pairs = pd.DataFrame(columns=["tf_clean", "gene", "weight", "tf"])

    if not include_before_cells:
        return (
            pd.DataFrame(columns=["receptor", "tf", "weight"]),
            pd.DataFrame(columns=["tf_clean", "gene", "weight"]),
            pd.DataFrame(columns=["ligand", "receptor", "receptor_clean", "weight"]),
            source_receptor_tf_pairs,
            source_tf_gene_pairs,
        )

    receptor_cells = ligand_cells if ligand_cells else None
    cc_df = get_cell_communication_layer(
        multicell_obj,
        as_dataframe=True,
        ligand_cells=[cell_type],
        receptor_cells=receptor_cells,
    )
    cc_df = cc_df[cc_df["ligand"].isin(source_ligand_candidates)].copy()
    if cc_df.empty:
        return (
            source_receptor_tf_pairs,
            source_tf_gene_pairs,
            pd.DataFrame(columns=["ligand", "receptor", "receptor_clean", "weight"]),
            pd.DataFrame(columns=["receptor", "tf", "weight"]),
            pd.DataFrame(columns=["tf_clean", "gene", "weight"]),
        )

    top_ligands = _rank_ligands_for_downstream(
        results,
        list(pd.unique(cc_df["ligand"])),
        top_ligand_n,
        per_celltype,
    )
    top_ligand_nodes = set(top_ligands["node_std"].astype(str)) if "node_std" in top_ligands.columns else set()
    cc_df = cc_df[cc_df["ligand"].isin(top_ligand_nodes)]

    connected_receptors = set(cc_df["receptor_clean"])
    receptors_with_tf_targets = set()
    receptor_tf_by_cell = {}
    top_tfs_by_cell = {}
    for downstream_cell in sorted({_node_celltype(r) for r in connected_receptors}):
        if not downstream_cell:
            continue
        receptor_tf_df = get_celltype_grn_receptor_bipartite(
            multicell_obj=multicell_obj,
            cell_type=downstream_cell,
            as_dataframe=True,
        )
        candidate_tfs = set(receptor_tf_df.loc[
            receptor_tf_df["col2"].isin(connected_receptors),
            "col1"
        ])
        top_tfs_for_cell = _rank_nodes(
            results,
            {tf for tf in candidate_tfs if "_TF::" in str(tf)} |
            {_as_tf_node(tf) for tf in candidate_tfs},
            len(results),
        )
        receptor_tf_df = _filter_receptor_targets_to_tfs(receptor_tf_df, top_tfs_for_cell)
        receptor_tf_df = receptor_tf_df[
            receptor_tf_df["col2"].isin(connected_receptors)
        ]
        receptors_with_tf_targets.update(receptor_tf_df["col2"])
        receptor_tf_by_cell[downstream_cell] = receptor_tf_df
        top_tfs_by_cell[downstream_cell] = top_tfs_for_cell

    top_receptors = _rank_downstream_receptors(
        results,
        receptors_with_tf_targets,
        top_receptor_n,
        per_celltype,
    )
    ligand_receptor_top = extract_receptor_ligand_pairs(
        receptor_ligand_df=cc_df,
        top_ligands_df=top_ligands,
        top_receptors_df=top_receptors,
    )

    all_receptor_tf_pairs = []
    all_gene_tf_pairs = []
    for downstream_cell in sorted({_node_celltype(r) for r in top_receptors["node"]}):
        if not downstream_cell:
            continue
        receptor_tf_df = receptor_tf_by_cell.get(downstream_cell)
        if receptor_tf_df is None:
            receptor_tf_df = get_celltype_grn_receptor_bipartite(
                multicell_obj=multicell_obj,
                cell_type=downstream_cell,
                as_dataframe=True,
            )
            receptor_tf_df = _filter_receptor_targets_to_tfs(
                receptor_tf_df,
                top_tfs_by_cell.get(downstream_cell, pd.DataFrame(columns=results.columns)),
            )
        connected_tfs = set(receptor_tf_df.loc[
            receptor_tf_df["col2"].isin(set(top_receptors["node"])),
            "col1"
        ])
        top_tfs = _rank_nodes(
            results,
            {tf for tf in connected_tfs if "_TF::" in str(tf)} |
            {_as_tf_node(tf) for tf in connected_tfs},
            top_tf_n,
        )
        receptor_tf_pairs = extract_receptor_tf_pairs(
            receptor_gene_layer=receptor_tf_df,
            top_tfs=top_tfs,
            top_receptors=top_receptors[top_receptors["node"].str.endswith(f"::{downstream_cell}")],
        )
        all_receptor_tf_pairs.append(receptor_tf_pairs)

        tf_gene_df = get_celltype_gene_layer(
            multicell_obj=multicell_obj,
            cell_type=downstream_cell,
            layer_name="gene",
            as_dataframe=True,
        )
        tf_nodes, tf_nodes_clean = _tf_name_sets(top_tfs)
        connected_genes = set(tf_gene_df.loc[
            tf_gene_df["source"].isin(tf_nodes) |
            tf_gene_df["source"].isin(tf_nodes_clean),
            "target"
        ])
        top_genes = _rank_nodes(results, connected_genes, top_tf_n)
        all_gene_tf_pairs.append(_extract_downstream_gene_tf_pairs(tf_gene_df, top_tfs, top_genes))

    receptor_tf_pairs = (
        pd.concat(all_receptor_tf_pairs, ignore_index=True)
        if all_receptor_tf_pairs else pd.DataFrame(columns=["receptor", "tf", "weight"])
    )
    gene_tf_pairs = (
        pd.concat(all_gene_tf_pairs, ignore_index=True)
        if all_gene_tf_pairs else pd.DataFrame(columns=["tf_clean", "gene", "weight"])
    )
    return (
        source_receptor_tf_pairs,
        source_tf_gene_pairs,
        ligand_receptor_top,
        receptor_tf_pairs,
        gene_tf_pairs,
    )


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
    cc_scores = cc_scores.sort_values(by="score", ascending=False)
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
    tf_gene_layer, top_tfs, seeds, verbose: bool = False
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
    verbose : bool, default=False
        If True, print debugging information about filtering.
    
    Returns
    -------
    pandas.DataFrame
        A filtered DataFrame containing only those rows of tf_gene_layer where:
            tf ∈ standardized(top_tfs["node"])  AND
            gene ∈ seeds.
    """

    genes_list = pd.Index(seeds).astype(str).tolist()
    tfs_list = list(top_tfs["node"].values)

    if verbose:
        print(f"\n[extract_gene_tf_pairs] === INPUT ===")
        print(f"  Input genes: {len(genes_list)} - {genes_list[:3]}")
        print(f"  Input TFs: {len(tfs_list)} - {tfs_list[:3]}")
        print(f"\n[extract_gene_tf_pairs] === TF-GENE LAYER ===")
        print(f"  Shape: {tf_gene_layer.shape}")
        print(f"  Columns: {tf_gene_layer.columns.tolist()}")
        print(f"  Unique sources (TFs): {tf_gene_layer['source'].nunique()}")
        print(f"  Unique targets (genes): {tf_gene_layer['target'].nunique()}")
        print(f"  Sample rows (first 3):")
        print(tf_gene_layer.head(3).to_string(index=False))

    # Support both TF -> gene and older gene -> TF GRN orientations.
    source = tf_gene_layer["source"].astype(str)
    target = tf_gene_layer["target"].astype(str)
    tf_to_gene = tf_gene_layer[
        source.isin(tfs_list) & target.isin(genes_list)
    ].copy()
    gene_to_tf = tf_gene_layer[
        source.isin(genes_list) & target.isin(tfs_list)
    ].copy()
    if len(tf_to_gene):
        filtered_df = tf_to_gene.rename(columns={"source": "tf", "target": "gene"})
    elif len(gene_to_tf):
        filtered_df = gene_to_tf.rename(columns={"source": "gene", "target": "tf"})
    else:
        filtered_df = tf_to_gene.rename(columns={"source": "tf", "target": "gene"})
    
    if verbose:
        print(f"\n[extract_gene_tf_pairs] === FILTERING RESULT ===")
        print(f"  Filtered pairs: {len(filtered_df)} rows")
        if len(filtered_df) == 0:
            print(f"  ⚠️  WARNING: No TF-gene pairs found!")
            # Check overlap
            tfs_in_layer = set(tf_gene_layer['source'].unique())
            genes_in_layer = set(tf_gene_layer['target'].unique())
            tf_overlap = set(tfs_list).intersection(tfs_in_layer)
            gene_overlap = set(genes_list).intersection(genes_in_layer)
            print(f"  TFs found in layer: {len(tf_overlap)}/{len(tfs_list)}")
            if tf_overlap:
                print(f"    Examples: {list(tf_overlap)[:3]}")
            print(f"  Genes found in layer: {len(gene_overlap)}/{len(genes_list)}")
            if gene_overlap:
                print(f"    Examples: {list(gene_overlap)[:3]}")
        else:
            print(f"  Sample filtered pairs (first 3):")
            print(filtered_df[['source', 'target', 'weight']].head(3).to_string(index=False))
    
    if "tf" not in filtered_df.columns and "source" in filtered_df.columns:
        filtered_df = filtered_df.rename(columns={"source": "tf", "target": "gene"})
    filtered_df.loc[:, 'tf_clean'] = filtered_df['tf'].str.replace('_TF', '', regex=False)

    return filtered_df


def extract_receptor_tf_pairs(
    receptor_gene_layer,
    top_tfs,
    top_receptors,
    verbose: bool = False
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
    verbose : bool, default=False
        If True, print debugging information about filtering.

    Returns
    -------
    pandas.DataFrame
        A filtered DataFrame containing only those rows of receptor_gene_layer where:
            tf ∈ standardized(top_tfs["node"])  AND
            receptor ∈ standardized(top_receptors["node"]).
    """
    tfs_list = list(top_tfs["node"].values)
    tfs_list_clean = ["".join(e.split("_TF")) for e in tfs_list]
    
    receptors_list = list(top_receptors["node"].values)

    if verbose:
        print(f"\n[extract_receptor_tf_pairs] === INPUT ===")
        print(f"  Input TFs: {len(tfs_list)} - {tfs_list[:3]}")
        print(f"  Input TFs (cleaned): {tfs_list_clean[:3]}")
        print(f"  Input receptors: {len(receptors_list)} - {receptors_list[:3]}")
        print(f"\n[extract_receptor_tf_pairs] === BIPARTITE LAYER ===")
        print(f"  Shape: {receptor_gene_layer.shape}")
        print(f"  Columns: {receptor_gene_layer.columns.tolist()}")
        print(f"  Unique col1 (genes/TFs): {receptor_gene_layer['col1'].nunique()}")
        print(f"  Unique col2 (receptors): {receptor_gene_layer['col2'].nunique()}")
        print(f"  Sample rows (first 3):")
        print(receptor_gene_layer.head(3).to_string(index=False))

    # FIXED: col2=receptors, col1=genes/TFs (actual data structure!)
    filtered_df = receptor_gene_layer[
        receptor_gene_layer.loc[:, "col2"].isin(receptors_list) &
        (
            receptor_gene_layer.loc[:, "col1"].isin(tfs_list_clean) |
            receptor_gene_layer.loc[:, "col1"].isin(tfs_list)
        )
    ].copy()

    if verbose:
        print(f"\n[extract_receptor_tf_pairs] === FILTERING RESULT ===")
        print(f"  Filtered pairs: {len(filtered_df)} rows")
        if len(filtered_df) == 0:
            print(f"  ⚠️  WARNING: No receptor-TF pairs found!")
            # Check overlap
            genes_in_layer = set(receptor_gene_layer['col1'].unique())
            receptors_in_layer = set(receptor_gene_layer['col2'].unique())
            tf_overlap = set(tfs_list_clean).intersection(genes_in_layer)
            receptor_overlap = set(receptors_list).intersection(receptors_in_layer)
            print(f"  TFs (genes) found in layer: {len(tf_overlap)}/{len(tfs_list_clean)}")
            if tf_overlap:
                print(f"    Examples: {list(tf_overlap)[:3]}")
            print(f"  Receptors found in layer: {len(receptor_overlap)}/{len(receptors_list)}")
            if receptor_overlap:
                print(f"    Examples: {list(receptor_overlap)[:3]}")
        else:
            print(f"  Sample filtered pairs (first 3):")
            print(filtered_df[['col1', 'col2', 'weight']].head(3).to_string(index=False))

    filtered_df = filtered_df.rename(columns={"col2": "receptor", "col1": "tf"})
    filtered_df.loc[:, "tf_clean"] = filtered_df["tf"].str.replace("_TF", "", regex=False)
    return filtered_df


def extract_receptor_ligand_pairs(
    receptor_ligand_df: pd.DataFrame,
    top_ligands_df: pd.DataFrame,
    top_receptors_df: pd.DataFrame,
    verbose: bool = False
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
    verbose : bool, default=False
        If True, print debugging information about filtering.

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

    if verbose:
        print(f"\n[extract_receptor_ligand_pairs] === INPUT ===")
        print(f"  Input ligands: {len(top_ligands_df)}")
        if len(top_ligands_df) > 0:
            print(f"    Examples: {top_ligands_df['node'].head(3).tolist()}")
        print(f"  Input receptors: {len(top_receptors_df)}")
        if len(top_receptors_df) > 0:
            print(f"    Examples: {top_receptors_df['node'].head(3).tolist()}")
        print(f"\n[extract_receptor_ligand_pairs] === RECEPTOR-LIGAND LAYER ===")
        print(f"  Shape: {receptor_ligand_df.shape}")
        print(f"  Columns: {receptor_ligand_df.columns.tolist()}")
        print(f"  Sample rows (first 3):")
        print(receptor_ligand_df.head(3).to_string(index=False))

    # 3) Filter receptor_ligand_df:
    #    Keep only rows where ligand ∈ standardized_ligs AND receptor ∈ standardized_recs
    mask = (
        receptor_ligand_df.loc[:, "ligand"].isin(standardized_ligs) &
        receptor_ligand_df.loc[:, "receptor"].isin(standardized_recs)
    )
    
    filtered = receptor_ligand_df.loc[mask].reset_index(drop=True)
    
    if verbose:
        print(f"\n[extract_receptor_ligand_pairs] === FILTERING RESULT ===")
        print(f"  Filtered pairs: {len(filtered)} rows")
        if len(filtered) == 0:
            print(f"  ⚠️  WARNING: No receptor-ligand pairs found!")
            # Check overlap
            ligands_in_layer = set(receptor_ligand_df['ligand'].unique())
            receptors_in_layer = set(receptor_ligand_df['receptor'].unique())
            ligand_overlap = standardized_ligs.intersection(ligands_in_layer)
            receptor_overlap = standardized_recs.intersection(receptors_in_layer)
            print(f"  Ligands found in layer: {len(ligand_overlap)}/{len(standardized_ligs)}")
            if ligand_overlap:
                print(f"    Examples: {list(ligand_overlap)[:3]}")
            print(f"  Receptors found in layer: {len(receptor_overlap)}/{len(standardized_recs)}")
            if receptor_overlap:
                print(f"    Examples: {list(receptor_overlap)[:3]}")
        else:
            print(f"  Sample filtered pairs (first 3):")
            print(filtered[['ligand', 'receptor', 'weight']].head(3).to_string(index=False))
    
    return filtered


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
    include_before_cells: bool = False,
    verbose: bool = False,
    direction: str = "upstream"
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
    direction = _normalize_flow(direction)
    if direction == "downstream":
        return _build_downstream_networks(
            multicell_obj=multicell_obj,
            results=results,
            cell_type=cell_type,
            seeds=seeds,
            ligand_cells=ligand_cells,
            top_ligand_n=top_ligand_n,
            top_receptor_n=top_receptor_n,
            top_tf_n=top_tf_n,
            per_celltype=per_celltype,
            include_before_cells=include_before_cells,
        )

    seeds_prefixed = _seed_names_for_cell_type(seeds, cell_type)
    cc_df = get_cell_communication_layer(
        multicell_obj,
        as_dataframe=True,
        ligand_cells=ligand_cells,
        receptor_cells=[cell_type])

    top_ligands = get_top_ligands(results, cc_df, n=top_ligand_n, per_celltype=per_celltype)

    tf_gene_df = get_celltype_gene_layer(
        multicell_obj=multicell_obj,
        cell_type=cell_type,
        layer_name="gene",
        as_dataframe=True
    )

    receptor_tf_df = get_celltype_grn_receptor_bipartite(
        multicell_obj=multicell_obj,
        cell_type=cell_type,
        as_dataframe=True
    )

    top_tfs = _top_seed_connected_tfs(
        results, tf_gene_df, seeds_prefixed, cell_type, top_tf_n
    )
    top_receptors = _top_tf_connected_receptors(
        results,
        receptor_tf_df,
        top_tfs,
        cell_type,
        top_receptor_n,
        receptor_candidates=set(cc_df["receptor_clean"]),
    )

    receptor_ligand_top = extract_receptor_ligand_pairs(
        receptor_ligand_df=cc_df,
        top_ligands_df=top_ligands,
        top_receptors_df=top_receptors,
        verbose=verbose
    )

    gene_tf_pairs = extract_gene_tf_pairs(
        tf_gene_df,
        top_tfs,
        seeds_prefixed,
        verbose=verbose,
    )
    receptor_tf_pairs = extract_receptor_tf_pairs(receptor_tf_df, top_tfs, top_receptors, verbose=verbose)

    if not include_before_cells:
        # Return empty frames for before-layers
        return (
            pd.DataFrame(columns=["receptor", "tf", "weight"]),
            pd.DataFrame(columns=["tf_clean", "gene", "weight"]),
            receptor_ligand_top,
            receptor_tf_pairs,
            gene_tf_pairs
        )

    # Otherwise build full before-cell layers. Infer sender cells from the
    # ligand node suffix; some saved objects have stale/inverted metadata in
    # ``celltype_source`` / ``celltype_target``.
    receptor_ligand_top = receptor_ligand_top.copy()
    if len(receptor_ligand_top):
        receptor_ligand_top.loc[:, "ligand_celltype"] = (
            receptor_ligand_top["ligand"].astype(str).apply(_node_celltype)
        )
    before_cell_types = [
        ct for ct in pd.unique(receptor_ligand_top.get("ligand_celltype", pd.Series(dtype=object)))
        if ct and ct != cell_type
    ]
    all_before_receptor_tf_pairs = []
    all_before_gene_tf_pairs = []

    for before_cell_type in before_cell_types:
        before_ligands = receptor_ligand_top.loc[
            receptor_ligand_top["ligand_celltype"] == before_cell_type, "ligand"
        ].values
        if len(before_ligands) == 0:
            continue

        before_tf_gene_df = get_celltype_gene_layer(multicell_obj, before_cell_type, "gene", as_dataframe=True)
        before_receptor_tf_df = get_celltype_grn_receptor_bipartite(multicell_obj, before_cell_type, as_dataframe=True)

        before_top_tfs = _top_seed_connected_tfs(
            results, before_tf_gene_df, before_ligands, before_cell_type, before_top_n
        )
        before_gene_tf_pairs = extract_gene_tf_pairs(
            tf_gene_layer=before_tf_gene_df,
            top_tfs=before_top_tfs,
            seeds=before_ligands,
        )

        before_top_receptors = _top_tf_connected_receptors(
            results, before_receptor_tf_df, before_top_tfs, before_cell_type, before_top_n
        )
        before_receptor_tf_pairs = extract_receptor_tf_pairs(
            receptor_gene_layer=before_receptor_tf_df,
            top_tfs=before_top_tfs,
            top_receptors=before_top_receptors,
        )

        all_before_receptor_tf_pairs.append(before_receptor_tf_pairs)
        all_before_gene_tf_pairs.append(before_gene_tf_pairs)

    all_before_receptor_tf_df = (
        pd.concat(all_before_receptor_tf_pairs, ignore_index=True)
        if all_before_receptor_tf_pairs
        else pd.DataFrame(columns=["receptor", "tf", "weight"])
    )
    all_before_gene_tf_df = (
        pd.concat(all_before_gene_tf_pairs, ignore_index=True)
        if all_before_gene_tf_pairs
        else pd.DataFrame(columns=["tf_clean", "gene", "weight"])
    )

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
    flow = _normalize_flow(flow)

    def format_links(df, source_col, target_col):
        return df.loc[:, [source_col, target_col, "weight"]].rename(columns={
            source_col: "source",
            target_col: "target",
            "weight": "value"
        })

    receptor_tf_col = "tf_clean" if flow == "downstream" and "tf_clean" in receptor_tf_df.columns else "tf"
    r_t = format_links(receptor_tf_df, "receptor", receptor_tf_col)
    t_g = format_links(gene_tf_df, "tf_clean", "gene")
    if flow == "downstream":
        _, r_t, t_g = _normalize_downstream_plot_links(
            pd.DataFrame(columns=["source", "target", "value"]), r_t, t_g
        )

    if len(r_t) > 0 and len(t_g) > 0:
        r_t = r_t[r_t["target"].isin(t_g["source"])]
        t_g = t_g[t_g["source"].isin(r_t["target"])]

    if len(r_t) == 0 and len(t_g) == 0:
        print("Warning: All networks are empty. Cannot generate Sankey plot.")
        return

    r_t["color"] = color
    t_g["color"] = color

    for df in [r_t, t_g]:
        total = df["value"].sum()
        if total > 0:
            df["value"] /= total

    links = pd.concat([df for df in [r_t, t_g] if len(df) > 0], ignore_index=True)

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
        if flow == "upstream"
        else f"Top regulators of downstream genes in {cell_type}: Receptor → TF → Gene"
    )
    layer_names = ["Receptors", "TFs", "Genes"]
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

    fig.show()


def plot_4layer_sankey(
    ligand_receptor_df: pd.DataFrame,
    receptor_tf_df: pd.DataFrame,
    gene_tf_df: pd.DataFrame,
    flow: str = "upstream",
    save_path: Union[str, None] = None
):
    """
    Plot a 4-layer Sankey diagram showing ligand → receptor → TF → gene.
    
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
        Direction of the result scores. The diagram keeps the same
        ligand → receptor → TF → gene layout for both directions.
    save_path : str or None, default=None
        If provided, save the figure as an HTML file to this path.

    Returns
    -------
    None
    """
    flow = _normalize_flow(flow)

    def hex_to_rgba(hex_color, alpha=0.6):
        h = hex_color.lstrip("#")
        return f"rgba({int(h[0:2], 16)}, {int(h[2:4], 16)}, {int(h[4:6], 16)}, {alpha})"

    def string_to_color(string):
        h = hashlib.md5(string.encode()).hexdigest()
        return "#" + h[:6]

    def assign_group_colors(df, column):
        unique_types = df[column].str.extract(r"::(.+)$")[0].fillna("Unknown")
        return unique_types.apply(lambda ct: hex_to_rgba(string_to_color(ct)))

    l_r = _format_sankey_links(ligand_receptor_df, "ligand", "receptor_clean")
    receptor_tf_col = "tf_clean" if flow == "downstream" and "tf_clean" in receptor_tf_df.columns else "tf"
    r_t = _format_sankey_links(receptor_tf_df, "receptor", receptor_tf_col)
    t_g = _format_sankey_links(gene_tf_df, "tf_clean", "gene")
    if flow == "downstream":
        l_r, r_t, t_g = _normalize_downstream_plot_links(l_r, r_t, t_g)

    l_r = l_r[l_r["target"].isin(r_t["source"])]
    r_t = r_t[r_t["source"].isin(l_r["target"])]
    t_g = t_g[t_g["source"].isin(r_t["target"])]

    if flow == "downstream":
        if len(l_r) > 0:
            l_r.loc[:, "color"] = assign_group_colors(l_r, "source")
        if len(r_t) > 0:
            r_t.loc[:, "color"] = assign_group_colors(r_t, "source")
        if len(t_g) > 0:
            t_g.loc[:, "color"] = assign_group_colors(t_g, "source")
    else:
        if len(l_r) > 0:
            l_r["color"] = assign_group_colors(l_r, "source")
        if len(r_t) > 0:
            r_t["color"] = "rgba(100,200,100,1)"
        if len(t_g) > 0:
            t_g["color"] = "rgba(100,200,100,1)"

    _normalize_link_values([l_r, r_t, t_g], flow)

    # Filter out empty DataFrames
    non_empty_dfs = [df for df in [l_r, r_t, t_g] if len(df) > 0]
    
    if len(non_empty_dfs) == 0:
        print("Warning: All networks are empty. Cannot generate Sankey plot.")
        return
    
    def _format_label(x: str) -> str:
        parts = x.split("::", 1)
        return parts[0] if len(parts) == 2 else x

    if flow == "downstream":
        links, raw_labels, node_x, node_y = _build_layered_links(
            [
                (l_r, "ligand", "receptor"),
                (r_t, "receptor", "tf"),
                (t_g, "tf", "gene"),
            ],
            {"ligand": 0.0, "receptor": 0.33, "tf": 0.66, "gene": 1.0},
        )
        labels = [_format_label(n) for n in raw_labels]
        node_extra = {"x": node_x, "y": node_y}
    else:
        links = pd.concat(non_empty_dfs, ignore_index=True)
        all_nodes = pd.unique(links[["source", "target"]].values.ravel())
        node_idx = {name: i for i, name in enumerate(all_nodes)}
        links["source_idx"] = links["source"].map(node_idx)
        links["target_idx"] = links["target"].map(node_idx)
        labels = [_format_label(n) for n in all_nodes]
        node_extra = {}

    sankey_data = go.Sankey(
        arrangement="fixed" if flow == "downstream" else "snap",
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=labels,
            **node_extra,
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
        if flow == "upstream"
        else f"Top regulators of downstream genes from other cell types:\n Ligand → Receptor → TF → Gene"
    )
    layer_names = ["Ligands", "Receptors", "TFs", "Genes"]
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

    fig.show()


def plot_6layer_sankey(
    before_receptor_tf_df: pd.DataFrame,
    before_tf_ligand_df: pd.DataFrame,
    ligand_receptor_df: pd.DataFrame,
    receptor_tf_df: pd.DataFrame,
    gene_tf_df: pd.DataFrame,
    flow: str = "upstream",
    save_path: Union[str, None] = None
):
    flow = _normalize_flow(flow)

    def hex_to_rgba(hex_color, alpha=0.6):
        h = hex_color.lstrip("#")
        return f"rgba({int(h[0:2], 16)}, {int(h[2:4], 16)}, {int(h[4:6], 16)}, {alpha})"

    def string_to_color(string):
        h = hashlib.md5(string.encode()).hexdigest()
        return "#" + h[:6]

    br_bt, bt_l, l_r, r_t, t_g = _prepare_6layer_links(
        before_receptor_tf_df,
        before_tf_ligand_df,
        ligand_receptor_df,
        receptor_tf_df,
        gene_tf_df,
        flow,
    )

    # Assign colors by celltype
    def assign_group_colors(df, column):
        unique_types = df[column].str.extract(r"::(.+)$")[0].fillna("Unknown")
        return unique_types.apply(lambda ct: hex_to_rgba(string_to_color(ct)))

    if flow == "downstream":
        if len(br_bt) > 0:
            br_bt.loc[:, "color"] = assign_group_colors(br_bt, "source")
        if len(bt_l) > 0:
            bt_l.loc[:, "color"] = assign_group_colors(bt_l, "target")
        if len(l_r) > 0:
            l_r.loc[:, "color"] = assign_group_colors(l_r, "source")
        if len(r_t) > 0:
            r_t.loc[:, "color"] = assign_group_colors(r_t, "source")
        if len(t_g) > 0:
            t_g.loc[:, "color"] = assign_group_colors(t_g, "source")
    else:
        if len(br_bt) > 0:
            br_bt.loc[:, "color"] = assign_group_colors(br_bt, "source")
        if len(bt_l) > 0:
            bt_l.loc[:, "color"] = assign_group_colors(bt_l, "source")
        if len(l_r) > 0:
            l_r.loc[:, "color"] = "rgba(160,160,160,0.4)"
        if len(r_t) > 0:
            r_t.loc[:, "color"] = "rgba(100,200,100,1)"
        if len(t_g) > 0:
            t_g.loc[:, "color"] = "rgba(100,200,100,1)"

    _normalize_link_values([br_bt, bt_l, l_r, r_t, t_g], flow)

    # Filter out empty DataFrames before concatenation
    non_empty_dfs = [df for df in [br_bt, bt_l, l_r, r_t, t_g] if len(df) > 0]
    
    if len(non_empty_dfs) == 0:
        print("Warning: All networks are empty. Cannot generate Sankey plot.")
        return
    
    def _format_label(x: str) -> str:
        parts = x.split("::", 1)
        return parts[0].split("_")[0] if len(parts) == 2 else x

    if flow == "downstream":
        links, raw_labels, node_x, node_y = _build_layered_links(
            [
                (br_bt, "sender_receptor", "sender_tf"),
                (bt_l, "sender_tf", "ligand"),
                (l_r, "ligand", "receiver_receptor"),
                (r_t, "receiver_receptor", "receiver_tf"),
                (t_g, "receiver_tf", "gene"),
            ],
            {
                "sender_receptor": 0.0,
                "sender_tf": 0.2,
                "ligand": 0.4,
                "receiver_receptor": 0.6,
                "receiver_tf": 0.8,
                "gene": 1.0,
            },
        )
        labels = [_format_label(n) for n in raw_labels]
        node_extra = {"x": node_x, "y": node_y}
    else:
        links = pd.concat(non_empty_dfs, ignore_index=True)
        all_nodes = pd.unique(links[["source", "target"]].values.ravel())
        node_idx = {name: i for i, name in enumerate(all_nodes)}
        links.loc[:, "source_idx"] = links.loc[:, "source"].map(node_idx)
        links.loc[:, "target_idx"] = links.loc[:, "target"].map(node_idx)
        labels = [_format_label(n) for n in all_nodes]
        node_extra = {}

    sankey_data = go.Sankey(
        arrangement="fixed" if flow == "downstream" else "snap",
        node=dict(
            pad=15,
            thickness=20,
            line=dict(color="black", width=0.5),
            label=labels,
            **node_extra,
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
        pd.concat([br_bt, bt_l, l_r, r_t, t_g]).loc[:, ["source", "color"]]
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

    if flow == "upstream":
        title_text = "Upstream Receptor → Upstream TF → Ligand → Downstream Receptor → Downstream TF → Gene"
        layer_names = ["Sender Receptors", "Sender TFs", "Ligands", "Receiver Receptors", "Receiver TFs", "Genes"]
    else:
        title_text = "Regulators of downstream genes: Upstream Receptor → Upstream TF → Ligand → Receiver Receptor → Receiver TF → Gene"
        layer_names = ["Sender Receptors", "Sender TFs", "Ligands", "Receiver Receptors", "Receiver TFs", "Genes"]

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

    fig.show()


def plot_intracell_sankey(
    multicell_obj,
    results,
    cell_type,
    seeds,
    top_receptor_n: int = 30,
    top_tf_n: int = 10,
    flow="upstream",
    verbose: bool = False,
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
        Direction of the result scores. The diagram keeps the same
        receptor → TF → gene layout for both directions.
    verbose : bool, default=False
        If True, print detailed debugging information about network construction.
    save_path : str or None, default=None
        If provided, save the figure as an HTML file to this path.

    Returns
    -------
    None
    """
    flow = _normalize_flow(flow)
    networks = build_partial_networks(
        multicell_obj=multicell_obj,
        results=results,
        cell_type=cell_type,
        seeds=seeds,
        ligand_cells=[],
        top_receptor_n=top_receptor_n,
        top_tf_n=top_tf_n,
        include_before_cells=False,
        verbose=verbose,
        direction=flow
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
    verbose: bool = False,
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
        Direction of the result scores. The diagram keeps the same
        ligand → receptor → TF → gene layout for both directions.
    verbose : bool, default=False
        If True, print detailed debugging information about network construction.
    save_path : str or None, default=None
        If provided, save the figure as an HTML file to this path.

    Returns
    -------
    None
    """
    flow = _normalize_flow(flow)
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
        include_before_cells=False,
        verbose=verbose,
        direction=flow
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
    verbose: bool = False,
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
        Direction of the result scores. The diagram keeps the same upstream
        receptor → upstream TF → ligand → receiver receptor → TF → gene layout.
    verbose : bool, default=False
        If True, print detailed debugging information about network construction.
    save_path : str or None, default=None
        If provided, save the figure as an HTML file to this path.

    Returns
    -------
    None
    """
    flow = _normalize_flow(flow)

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
        include_before_cells=True,
        verbose=verbose,
        direction=flow
    )

    networks = {
        "before_receptor_tf": networks[0],
        "before_tf_ligand": networks[1],
        "ligand_receptor": networks[2],
        "receptor_tf": networks[3],
        "gene_tf": networks[4],
    }

    if (
        flow == "downstream"
        and len(networks["before_receptor_tf"]) == 0
        and len(networks["before_tf_ligand"]) == 0
    ):
        plot_4layer_sankey(
            ligand_receptor_df=networks["ligand_receptor"],
            receptor_tf_df=networks["receptor_tf"],
            gene_tf_df=networks["gene_tf"],
            flow=flow,
            save_path=save_path
        )
    else:
        plot_6layer_sankey(
            before_receptor_tf_df=networks["before_receptor_tf"],
            before_tf_ligand_df=networks["before_tf_ligand"],
            ligand_receptor_df=networks["ligand_receptor"],
            receptor_tf_df=networks["receptor_tf"],
            gene_tf_df=networks["gene_tf"],
            flow=flow,
            save_path=save_path
        )
