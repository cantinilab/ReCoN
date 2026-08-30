"""
Biological cell cascade diagram — public API.

Two entry points
────────────────
  cascade_plot            — single-run diagram, score-based coloring
  contrast_cascade_plot   — two-run contrast, diverging delta coloring
"""

import warnings

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from .cascade_core import (
    # Network & layout
    build_networks,
    build_networks_contrast,
    trim_cascade_edges,
    empty_requested_cascade_edges,
    collect_node_sets,
    compute_global_geometry,
    setup_figure,
    draw_cells_and_edges,
    add_section_labels,
    normalize_flow,
    # Node helpers
    celltype_of,
    compute_node_radii,
    # Color helpers
    resolve_node_type_colors,
    score_to_rgba,
    score_to_grey,
    # Legend
    draw_type_legend,
)


_EMPTY_LAYER_HELP = {
    "upstream_r_tf": "sender receptor → TF (increase before_top_n)",
    "upstream_tf_lig": "sender TF → ligand (increase before_top_n or top_ligand_n)",
    "lig_rec": "ligand → receiver receptor (increase top_ligand_n or top_receptor_n)",
    "rec_tf": "receiver receptor → TF (increase top_receptor_n or top_tf_n)",
    "tf_gene": "receiver TF → gene (increase top_tf_n)",
}


def _warn_empty_requested_layers(edges, *, start_layer=None, stop_layer=None):
    empty = empty_requested_cascade_edges(
        edges, start_layer=start_layer, stop_layer=stop_layer,
    )
    if not empty:
        return
    transitions = "; ".join(_EMPTY_LAYER_HELP[key] for key in empty)
    warnings.warn(
        "No connections were retained for the requested cascade transition(s): "
        f"{transitions}. If increasing the suggested selection limit does not "
        "help, check gene identifiers and confirm that the required edges exist "
        "in the GRN, receptor prior, or communication network.",
        UserWarning,
        stacklevel=3,
    )


def _limit_displayed_seed_gene_edges(edges, seeds, cell_type, max_seed_nodes):
    """Limit displayed seed-gene targets without changing upstream selection."""
    if seeds is None or max_seed_nodes is None or "tf_gene" not in edges:
        return edges
    seed_values = list(seeds.keys()) if hasattr(seeds, "keys") else list(seeds)
    seed_nodes = []
    for seed in seed_values:
        seed = str(seed)
        seed_nodes.append(seed if "::" in seed else f"{seed}::{cell_type}")
    present_seed_nodes = [
        seed for seed in seed_nodes
        if seed in set(edges["tf_gene"].get("target", pd.Series(dtype=object)).astype(str))
    ]
    if len(present_seed_nodes) <= max_seed_nodes:
        return edges

    warnings.warn(
        f"Showing only the first {max_seed_nodes} seed genes out of {len(present_seed_nodes)} "
        "in the cascade plot. All seeds were still used for network selection; "
        "increase max_seed_nodes to show more.",
        UserWarning,
    )
    keep_seed_nodes = set(present_seed_nodes[:max_seed_nodes])
    limited = {key: value.copy() for key, value in edges.items()}
    tf_gene = limited["tf_gene"]
    seed_target_mask = tf_gene["target"].astype(str).isin(set(seed_nodes))
    limited["tf_gene"] = tf_gene[~seed_target_mask | tf_gene["target"].astype(str).isin(keep_seed_nodes)].copy()
    return limited

def _score_to_rgba_with_floor(score, vmax, base_rgb, alpha=0.85, min_fraction=0.22):
    if not np.isfinite(score) or vmax <= 0:
        return (1.0, 1.0, 1.0, alpha)
    t = min(max(score / vmax, 0.0), 1.0)
    if t > 0:
        t = min_fraction + t * (1.0 - min_fraction)
    r = 1.0 + t * (base_rgb[0] - 1.0)
    g = 1.0 + t * (base_rgb[1] - 1.0)
    b = 1.0 + t * (base_rgb[2] - 1.0)
    return (r, g, b, alpha)


def _score_to_grey_with_floor(score, vmax, alpha=0.65, min_fraction=0.25):
    if not np.isfinite(score) or vmax <= 0:
        return (1.0, 1.0, 1.0, alpha)
    t = min(max(score / vmax, 0.0), 1.0)
    if t > 0:
        t = min_fraction + t * (1.0 - min_fraction)
    grey = 1.0 - t * 0.55
    return (grey, grey, grey, alpha)


# ═══════════════════════════════════════════════════════════════════════════
# cascade_plot — single run, score-based coloring
# ═══════════════════════════════════════════════════════════════════════════

def cascade_plot(
    multicell_obj,
    results: pd.DataFrame,
    *,
    cell_type: str = None,
    seeds: list = None,
    ligand_cells: list = None,
    top_ligand_n: int = 20,
    top_receptor_n: int = 10,
    top_tf_n: int = 10,
    before_top_n: int = 5,
    per_celltype: bool = True,
    show_seeds: bool = False,
    max_seed_nodes: int = 20,
    celltype_colors: dict = None,
    node_type_colors: dict = None,
    celltype_display_names: dict = None,
    show_labels: bool = True,
    label_fontsize: int = 14,
    node_size_by_weight: bool = True,
    node_type_halo: bool = False,
    node_alpha: float = 0.85,
    edge_alpha: float = 0.65,
    seed_label: str = None,
    seed_label_fontsize: int = 10,
    figsize: tuple = None,
    save_path: str = None,
    verbose: bool = False,
    title: str = None,
    flow: str = "upstream",
    start_layer: str = None,
    stop_layer: str = None,
    return_edges: bool = False,
):
    """Biological cell cascade diagram for a single run.

    Node color: white → type color gradient (score-based).
    Edge color: white → grey gradient (source score-based),
                with black outline underneath.
    Node radius: proportional to score.

    Parameters
    ----------
    multicell_obj : single multicell object
    results : DataFrame with columns ``node``, ``score``.
    seeds : list of str, optional
        Raw seed gene names (e.g. ['Nfkb1', 'Tnf']).
    show_seeds : bool
        If True, seed genes are shown as individual nodes in the nucleus.
    flow : {'upstream', 'downstream'}, default='upstream'
        Controls how the cascade network is selected from ReCoN scores.
    start_layer, stop_layer : str, optional
        Inclusive biological layer window. Examples: ``stop_layer="ligands"``
        keeps sender regulation up to ligands; ``start_layer="ligands"``
        keeps ligand-to-receiver regulation.
    return_edges : bool, default=False
        If True, return ``(fig, ax, edges)`` for testing/debugging.
    node_type_colors : dict, optional
        Partial overrides for ``receptor``, ``tf``, ``ligand``, ``gene``,
        and ``seed``. Values may use any Matplotlib color format.
    """
    flow = normalize_flow(flow)
    node_type_colors = resolve_node_type_colors(node_type_colors)

    # ── 1. Build & filter networks ──────────────────────────────────────
    edges = build_networks(
        multicell_obj, results,
        cell_type=cell_type, seeds=seeds, ligand_cells=ligand_cells,
        top_ligand_n=top_ligand_n, top_receptor_n=top_receptor_n,
        top_tf_n=top_tf_n, before_top_n=before_top_n,
        per_celltype=per_celltype, verbose=verbose,
        flow=flow,
    )
    edges = trim_cascade_edges(edges, start_layer=start_layer, stop_layer=stop_layer)
    _warn_empty_requested_layers(
        edges, start_layer=start_layer, stop_layer=stop_layer,
    )
    edges = _limit_displayed_seed_gene_edges(edges, seeds, cell_type, max_seed_nodes)

    # ── 2. Collect node sets ────────────────────────────────────────────
    seed_display = seeds if show_seeds else None
    nodes = collect_node_sets(
        edges,
        seeds=seed_display,
        cell_type=cell_type,
    )

    if verbose:
        print(f"Sender CTs   : {list(dict.fromkeys(celltype_of(n) for n in nodes['ligands']))}")
        print(f"Sender recs  : {len(nodes['sender_recs'])}, TFs: {len(nodes['sender_tfs'])}")
        print(f"Ligands      : {len(nodes['ligands'])}")
        print(f"Recv recs    : {len(nodes['recv_recs'])}, TFs: {len(nodes['recv_tfs'])}, "
              f"genes: {len(nodes['recv_genes'])}")
        if show_seeds:
            print(f"Seeds        : {len(nodes['seed_nodes'])}")

    # ── 3. Score dict for coloring & sizing ─────────────────────────────
    score_dict = results.set_index("node")["score"].to_dict()

    # ── 4. Global geometry ──────────────────────────────────────────────
    geo = compute_global_geometry(nodes)

    # ── 5. Node radii ───────────────────────────────────────────────────
    all_drawn = (nodes["sender_recs"] + nodes["sender_tfs"] + nodes["ligands"]
                 + nodes["recv_recs"] + nodes["recv_tfs"] + nodes["recv_genes"]
                 + nodes["seed_nodes"])
    NODE_R = max(0.22, geo["RECV_R"] * 0.055)
    node_radii = compute_node_radii(
        all_drawn, score_dict, NODE_R, node_size_by_weight,
        group_lists=[nodes["sender_recs"], nodes["sender_tfs"],
                     nodes["ligands"], nodes["recv_recs"],
                     nodes["recv_tfs"], nodes["recv_genes"],
                     nodes["seed_nodes"]],
    )
    TF_SCALE = 0.75
    for n in nodes["sender_tfs"] + nodes["recv_tfs"] + nodes["sender_recs"]:
        if n in node_radii:
            node_radii[n] *= TF_SCALE

    # Cap seed & gene radii to fit inside the nucleus (with 10% margin)
    RECV_RY_RATIO = 0.85
    nuc_a = 0.40 * geo["RECV_R"] * 0.90
    nuc_b = 0.40 * geo["RECV_R"] * RECV_RY_RATIO * 0.90
    nuc_min_dim = min(nuc_a, nuc_b)
    nucleus_nodes = nodes["recv_genes"] + nodes["seed_nodes"]
    n_nucleus = max(len(nucleus_nodes), 1)
    max_nucleus_r = (nuc_min_dim * 0.9) / max(2.0, np.sqrt(n_nucleus) * 1.2)
    for n in nucleus_nodes:
        if n in node_radii and node_radii[n] > max_nucleus_r:
            node_radii[n] = max_nucleus_r

    # ── 6. Per-group score vmax for color gradient ──────────────────────
    _group_map = {}
    for n in nodes["sender_recs"]:
        _group_map[n] = "sender_receptor"
    for n in nodes["recv_recs"]:
        _group_map[n] = "receiver_receptor"
    for n in nodes["sender_tfs"]:
        _group_map[n] = "sender_tf"
    for n in nodes["recv_tfs"]:
        _group_map[n] = "receiver_tf"
    for n in nodes["ligands"]:
        _group_map[n] = "ligand"
    for n in nodes["recv_genes"]:
        _group_map[n] = "gene"
    for n in nodes["seed_nodes"]:
        _group_map[n] = "seed"

    _group_vmax = {}
    for group_name in ("sender_receptor", "receiver_receptor", "sender_tf", "receiver_tf", "ligand", "gene", "seed"):
        group_nodes = [n for n, g in _group_map.items() if g == group_name]
        group_scores = [score_dict.get(n, 0.0) for n in group_nodes]
        _group_vmax[group_name] = max(group_scores) if group_scores else 1.0
        if _group_vmax[group_name] <= 0:
            _group_vmax[group_name] = 1.0

    # ── 7. Color functions ──────────────────────────────────────────────
    _ROLE_BASE_RGB = {k: v[:3] for k, v in node_type_colors.items()}
    _DISPLAY_ROLE = {
        "sender_receptor": "receptor",
        "receiver_receptor": "receptor",
        "sender_tf": "tf",
        "receiver_tf": "tf",
    }

    def _node_fill(node_name, role):
        score_group = _group_map.get(node_name, role)
        display_role = _DISPLAY_ROLE.get(score_group, role)
        base_rgb = _ROLE_BASE_RGB.get(display_role, (0.5, 0.5, 0.5))
        s = score_dict.get(node_name, 0.0)
        vmax = _group_vmax.get(score_group, 1.0)
        return _score_to_rgba_with_floor(s, vmax, base_rgb, alpha=node_alpha)

    def _edge_color(node_ref, layer_key):
        s = score_dict.get(node_ref, 0.0)
        score_group = _group_map.get(node_ref, "gene")
        vmax = _group_vmax.get(score_group, 1.0)
        return _score_to_grey_with_floor(s, vmax, alpha=edge_alpha)

    # ── 8. Figure setup & draw ──────────────────────────────────────────
    fig, ax, y_max = setup_figure(geo, figsize)

    draw_cells_and_edges(
        ax, nodes, edges, geo,
        cell_type=cell_type,
        celltype_colors=celltype_colors,
        celltype_display_names=celltype_display_names,
        label_fontsize=label_fontsize,
        node_radii=node_radii,
        node_type_halo=node_type_halo,
        show_labels=show_labels,
        show_seeds=show_seeds,
        fill_fn=_node_fill,
        edge_color_fn=_edge_color,
        seed_label=seed_label,
        seed_label_fontsize=seed_label_fontsize,
        flow=flow,
        node_type_colors=node_type_colors,
    )

    # ── 9. Labels & legend ──────────────────────────────────────────────
    add_section_labels(ax, geo, y_max, label_fontsize)
    draw_type_legend(
        ax, label_fontsize, show_seeds=show_seeds,
        node_type_colors=node_type_colors,
    )

    # ── 10. Title ───────────────────────────────────────────────────────
    fig_title = title or f"Cell signaling cascade ({flow}) — {cell_type or '?'}"
    plt.title(fig_title, fontsize=label_fontsize + 6, pad=12)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, bbox_inches="tight", dpi=150)
        print(f"Saved: {save_path}")
    plt.show()
    if return_edges:
        return fig, ax, edges
    return fig, ax


# ═══════════════════════════════════════════════════════════════════════════
# contrast_cascade_plot — two-run comparison
# ═══════════════════════════════════════════════════════════════════════════

def contrast_cascade_plot(
    multicell_objs: dict,
    results_a: pd.DataFrame,
    results_b: pd.DataFrame,
    *,
    cell_type: str = None,
    seeds: list = None,
    ligand_cells: list = None,
    top_ligand_n: int = 20,
    top_receptor_n: int = 10,
    top_tf_n: int = 10,
    before_top_n: int = 5,
    per_celltype: bool = True,
    show_seeds: bool = False,
    max_seed_nodes: int = 20,
    delta_vmax: float = 2.0,
    delta_min_quantile: float = 0.3,
    delta_min_color_fraction: float = 0.0,
    contrast_scheme: str = "temperature",
    celltype_colors: dict = None,
    node_type_colors: dict = None,
    celltype_display_names: dict = None,
    normalize_receiver_scores: bool = True,
    show_labels: bool = True,
    label_fontsize: int = 14,
    node_size_by_weight: bool = True,
    node_type_halo: bool = False,
    node_alpha: float = 0.85,
    edge_alpha: float = 0.65,
    seed_label: str = None,
    seed_label_fontsize: int = 10,
    figsize: tuple = None,
    save_path: str = None,
    verbose: bool = False,
    title: str = None,
    flow: str = "upstream",
    start_layer: str = None,
    stop_layer: str = None,
    return_edges: bool = False,
):
    """Biological cell cascade diagram comparing two runs.

    Node and edge color uses a diverging scale based on the score
    difference between run A and run B.

    Parameters
    ----------
    multicell_objs : dict
        Condition-keyed dict of multicell objects.
    results_a : DataFrame
        Run A results with columns ``node``, ``score``.
    results_b : DataFrame
        Run B results with columns ``node``, ``score``.
    contrast_scheme : {'temperature', 'sex'}
        'temperature': blue (B-enriched) ↔ white ↔ red (A-enriched).
        'sex': blue (A-enriched) ↔ white ↔ orange (B-enriched).
    normalize_receiver_scores : bool
        Z-score color values within each receiver layer.
    node_type_colors : dict, optional
        Partial overrides for the molecular-type halo and legend colors.
        Contrast fill colors continue to represent the score difference.
    flow : {'upstream', 'downstream'}, default='upstream'
        If 'upstream', draw ligand → receptor → TF → gene arrows.
        If 'downstream', draw gene → TF → receptor → ligand arrows.
    """
    from .cascade_core import (
        delta_to_rgba,
        parse_rgba,
        boost_values_for_color,
        zscore_layer,
        draw_contrast_legend,
    )
    flow = normalize_flow(flow)
    node_type_colors = resolve_node_type_colors(node_type_colors)

    # ── 1. Compute delta metrics from two score vectors ─────────────────
    scores_a = results_a.set_index("node")["score"]
    scores_b = results_b.set_index("node")["score"]
    all_nodes_idx = scores_a.index.union(scores_b.index)
    scores_a = scores_a.reindex(all_nodes_idx, fill_value=0.0)
    scores_b = scores_b.reindex(all_nodes_idx, fill_value=0.0)

    signed_delta = (scores_a - scores_b).to_dict()
    abs_delta    = (scores_b.abs() - scores_a.abs()).to_dict()
    max_abs      = np.maximum(scores_a.abs(), scores_b.abs()).to_dict()

    # Use one of the results as scaffold
    scaffold = results_a.copy()

    # ── 2. Build & filter networks ──────────────────────────────────────
    edges = build_networks_contrast(
        multicell_objs, scaffold, signed_delta,
        cell_type=cell_type, seeds=seeds, ligand_cells=ligand_cells,
        top_ligand_n=top_ligand_n, top_receptor_n=top_receptor_n,
        top_tf_n=top_tf_n, before_top_n=before_top_n,
        per_celltype=per_celltype, verbose=verbose, flow=flow,
    )
    edges = trim_cascade_edges(edges, start_layer=start_layer, stop_layer=stop_layer)
    _warn_empty_requested_layers(
        edges, start_layer=start_layer, stop_layer=stop_layer,
    )
    edges = _limit_displayed_seed_gene_edges(edges, seeds, cell_type, max_seed_nodes)

    # ── 3. Collect node sets ────────────────────────────────────────────
    nodes = collect_node_sets(
        edges,
        seeds=seeds if show_seeds else None,
        cell_type=cell_type,
    )

    if verbose:
        print(f"Sender CTs   : {list(dict.fromkeys(celltype_of(n) for n in nodes['ligands']))}")
        print(f"Sender recs  : {len(nodes['sender_recs'])}, TFs: {len(nodes['sender_tfs'])}")
        print(f"Ligands      : {len(nodes['ligands'])}")
        print(f"Recv recs    : {len(nodes['recv_recs'])}, TFs: {len(nodes['recv_tfs'])}, "
              f"genes: {len(nodes['recv_genes'])}")
        if show_seeds:
            print(f"Seeds        : {len(nodes['seed_nodes'])}")

    # ── 4. Global geometry ──────────────────────────────────────────────
    geo = compute_global_geometry(nodes)

    # ── 5. Prepare color dict (abs_delta) ───────────────────────────────
    color_dict = abs_delta.copy()
    for node in list(color_dict.keys()):
        node_to_replace = "::".join(node.rsplit("-", 1))
        if node_to_replace not in color_dict:
            color_dict[node_to_replace] = color_dict[node]
        elif color_dict[node_to_replace] == 0:
            color_dict[node_to_replace] = color_dict[node]

    if normalize_receiver_scores:
        zscore_layer(color_dict, nodes["recv_recs"])
        zscore_layer(color_dict, nodes["recv_tfs"])
        zscore_layer(color_dict, nodes["recv_genes"])

    if delta_min_color_fraction > 0:
        color_dict = boost_values_for_color(
            color_dict, delta_vmax, delta_min_color_fraction)

    # ── 6. Node radii ───────────────────────────────────────────────────
    all_drawn = (nodes["sender_recs"] + nodes["sender_tfs"] + nodes["ligands"]
                 + nodes["recv_recs"] + nodes["recv_tfs"] + nodes["recv_genes"]
                 + nodes["seed_nodes"])
    NODE_R = max(0.22, geo["RECV_R"] * 0.055)
    node_radii = compute_node_radii(
        all_drawn, max_abs, NODE_R, node_size_by_weight,
        group_lists=[nodes["sender_recs"], nodes["sender_tfs"],
                     nodes["ligands"], nodes["recv_recs"],
                     nodes["recv_tfs"], nodes["recv_genes"],
                     nodes["seed_nodes"]],
    )
    TF_SCALE = 0.75
    for n in nodes["sender_tfs"] + nodes["recv_tfs"] + nodes["sender_recs"]:
        if n in node_radii:
            node_radii[n] *= TF_SCALE

    # ── 7. Color functions ──────────────────────────────────────────────
    def _node_fill(node_name, _role):
        val = color_dict.get(node_name, 0.0)
        return parse_rgba(
            delta_to_rgba(val, vmax=delta_vmax, alpha=node_alpha,
                          scheme=contrast_scheme)
        )

    def _edge_color(node_ref, layer_key):
        val = color_dict.get(node_ref, 0.0)
        return parse_rgba(
            delta_to_rgba(val, vmax=delta_vmax, alpha=edge_alpha,
                          scheme=contrast_scheme)
        )

    # ── 8. Delta thresholds for edge filtering ──────────────────────────
    delta_thresholds = {}
    if delta_min_quantile > 0:
        for key in ("upstream_r_tf", "upstream_tf_lig", "rec_tf", "tf_gene"):
            df = edges[key]
            if len(df):
                abs_d = df["source"].map(lambda n: abs(signed_delta.get(n, 0.0)))
                delta_thresholds[key] = float(abs_d.quantile(delta_min_quantile))

    # ── 9. Figure setup & draw ──────────────────────────────────────────
    fig, ax, y_max = setup_figure(geo, figsize)

    draw_cells_and_edges(
        ax, nodes, edges, geo,
        cell_type=cell_type,
        celltype_colors=celltype_colors,
        celltype_display_names=celltype_display_names,
        label_fontsize=label_fontsize,
        node_radii=node_radii,
        node_type_halo=node_type_halo,
        show_labels=show_labels,
        show_seeds=show_seeds,
        fill_fn=_node_fill,
        edge_color_fn=_edge_color,
        seed_label=seed_label,
        seed_label_fontsize=seed_label_fontsize,
        score_filter_dict=signed_delta,
        lfc_thresholds=delta_thresholds,
        flow=flow,
        node_type_colors=node_type_colors,
    )

    # ── 10. Labels & legend ─────────────────────────────────────────────
    add_section_labels(ax, geo, y_max, label_fontsize)
    draw_contrast_legend(ax, contrast_scheme, delta_vmax,
                         normalize_receiver_scores, label_fontsize,
                         show_seeds=show_seeds,
                         node_type_colors=node_type_colors)

    # ── 11. Title ───────────────────────────────────────────────────────
    if title is not None:
        fig_title = title
        if normalize_receiver_scores:
            fig_title += "  †norm"
    else:
        if contrast_scheme == "sex":
            suffix = "  [A (hi) vs B (lo)]"
        else:
            suffix = "  [warm (hi) vs cold (lo)]"
        if normalize_receiver_scores:
            suffix += "  †norm"
        fig_title = f"Cell signaling diagram ({flow}) — {cell_type or '?'}{suffix}"

    plt.title(fig_title, fontsize=label_fontsize + 6, pad=12)

    plt.tight_layout()
    if save_path:
        plt.savefig(save_path, bbox_inches="tight", dpi=150)
        print(f"Saved: {save_path}")
    plt.show()
    if return_edges:
        return fig, ax, edges
    return fig, ax
