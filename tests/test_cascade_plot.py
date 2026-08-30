"""Tests for recon.plot cascade plots."""
import importlib
import sys
import types
import pytest
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyArrowPatch, FancyBboxPatch

from recon.plot import cascade_core
cascade_plot_module = importlib.import_module("recon.plot.cascade_plot")

cascade_plot_fn = cascade_plot_module.cascade_plot
contrast_cascade_plot_fn = cascade_plot_module.contrast_cascade_plot


@pytest.fixture
def cascade_edges():
    return {
        "upstream_r_tf": pd.DataFrame({
            "source": ["SREC_receptor::SenderA", "DUP_receptor::Receiver"],
            "target": ["STF_TF::SenderA", "RTF_TF::Receiver"],
            "value": [0.5, 0.7],
        }),
        "upstream_tf_lig": pd.DataFrame({
            "source": ["STF_TF::SenderA", "STF_TF::SenderB"],
            "target": ["LIG1::SenderA", "LIG2::SenderB"],
            "value": [0.8, 0.9],
        }),
        "lig_rec": pd.DataFrame({
            "source": ["LIG1::SenderA", "LIG2::SenderB", "DUP_receptor::Receiver"],
            "target": ["RREC_receptor::Receiver", "RREC_receptor::Receiver", "DUP_receptor::Receiver"],
            "value": [1.0, 2.0, 3.0],
        }),
        "rec_tf": pd.DataFrame({
            "source": ["RREC_receptor::Receiver", "DUP_receptor::Receiver"],
            "target": ["RTF_TF::Receiver", "GENE_COLLIDE_TF::Receiver"],
            "value": [1.1, 1.2],
        }),
        "tf_gene": pd.DataFrame({
            "source": ["RTF_TF::Receiver", "GENE_COLLIDE_TF::Receiver", "RTF_TF::Receiver", "RTF_TF::Receiver"],
            "target": ["GENE1::Receiver", "GENE_COLLIDE::Receiver", "SEED1::Receiver", "SEED2::Receiver"],
            "value": [1.3, 1.4, 1.5, 1.6],
        }),
    }


@pytest.fixture
def cascade_results():
    return pd.DataFrame({
        "node": [
            "SREC_receptor::SenderA", "STF_TF::SenderA", "STF_TF::SenderB",
            "LIG1::SenderA", "LIG2::SenderB", "RREC_receptor::Receiver",
            "DUP_receptor::Receiver", "RTF_TF::Receiver",
            "GENE_COLLIDE_TF::Receiver", "GENE1::Receiver",
            "GENE_COLLIDE::Receiver", "SEED1::Receiver", "SEED2::Receiver",
        ],
        "score": [1.0, 2.0, 1.5, 3.0, 2.5, 4.0, 1.0, 5.0, 1.0, 2.0, 0.5, 6.0, 5.5],
    })


def _empty_edges():
    return {
        key: pd.DataFrame(columns=["source", "target", "value"])
        for key in ("upstream_r_tf", "upstream_tf_lig", "lig_rec", "rec_tf", "tf_gene")
    }


class TestCascadeFlow:
    """Test cascade plot flow direction handling."""

    def test_draw_single_edge_upstream_direction(self, monkeypatch):
        captured = []

        def mock_add_patch(patch):
            if isinstance(patch, FancyArrowPatch):
                captured.append(patch._posA_posB)

        _, ax = plt.subplots()
        monkeypatch.setattr(ax, "add_patch", mock_add_patch)

        cascade_core.draw_single_edge(
            ax=ax,
            src="LIGAND::Sender",
            tgt="RECEPTOR::Receiver",
            value=1.0,
            layer_key="lig_rec",
            rad=0.1,
            positions={
                "LIGAND::Sender": (0.0, 0.0),
                "RECEPTOR::Receiver": (1.0, 0.0),
            },
            lw_ranges={"lig_rec": (0.0, 1.0)},
            color_fn=lambda _node, _layer: (0.5, 0.5, 0.5, 0.5),
            flow="upstream",
        )

        assert captured[0] == [(0.0, 0.0), (1.0, 0.0)]
        plt.close("all")

    def test_draw_single_edge_downstream_keeps_topology_direction(self, monkeypatch):
        captured = []

        def mock_add_patch(patch):
            if isinstance(patch, FancyArrowPatch):
                captured.append(patch._posA_posB)

        _, ax = plt.subplots()
        monkeypatch.setattr(ax, "add_patch", mock_add_patch)

        cascade_core.draw_single_edge(
            ax=ax,
            src="LIGAND::Sender",
            tgt="RECEPTOR::Receiver",
            value=1.0,
            layer_key="lig_rec",
            rad=0.1,
            positions={
                "LIGAND::Sender": (0.0, 0.0),
                "RECEPTOR::Receiver": (1.0, 0.0),
            },
            lw_ranges={"lig_rec": (0.0, 1.0)},
            color_fn=lambda _node, _layer: (0.5, 0.5, 0.5, 0.5),
            flow="downstream",
        )

        assert captured[0] == [(0.0, 0.0), (1.0, 0.0)]
        plt.close("all")

    def test_invalid_flow_raises(self):
        with pytest.raises(ValueError, match="flow must be either"):
            cascade_core.normalize_flow("sideways")

    def test_flow_is_case_insensitive(self):
        assert cascade_core.normalize_flow("DownStream") == "downstream"


class TestCascadeLayerWindow:
    """Test layer-window controls for cascade plots."""

    def test_trim_cascade_edges_keeps_all_keys_and_empties_outside_window(self, cascade_edges):
        trimmed = cascade_core.trim_cascade_edges(cascade_edges, stop_layer="ligands")

        assert list(trimmed) == cascade_core.CASCADE_EDGE_ORDER
        assert len(trimmed["upstream_r_tf"]) == len(cascade_edges["upstream_r_tf"])
        assert len(trimmed["upstream_tf_lig"]) == len(cascade_edges["upstream_tf_lig"])
        assert trimmed["lig_rec"].empty
        assert trimmed["rec_tf"].empty
        assert trimmed["tf_gene"].empty
        assert trimmed["lig_rec"].columns.tolist() == ["source", "target", "value"]

    def test_trim_cascade_edges_can_start_at_ligands(self, cascade_edges):
        trimmed = cascade_core.trim_cascade_edges(cascade_edges, start_layer="ligands")

        assert trimmed["upstream_r_tf"].empty
        assert trimmed["upstream_tf_lig"].empty
        assert len(trimmed["lig_rec"]) == len(cascade_edges["lig_rec"])
        assert len(trimmed["rec_tf"]) == len(cascade_edges["rec_tf"])
        assert len(trimmed["tf_gene"]) == len(cascade_edges["tf_gene"])

    def test_trim_cascade_edges_validates_layer_order(self, cascade_edges):
        with pytest.raises(ValueError, match="start_layer must"):
            cascade_core.trim_cascade_edges(cascade_edges, start_layer="genes", stop_layer="ligands")
        with pytest.raises(ValueError, match="Unknown cascade layer"):
            cascade_core.trim_cascade_edges(cascade_edges, stop_layer="banana")

    def test_empty_requested_edges_ignore_transitions_outside_window(self, cascade_edges):
        trimmed = cascade_core.trim_cascade_edges(
            cascade_edges, start_layer="receptors", stop_layer="genes",
        )
        assert cascade_core.empty_requested_cascade_edges(
            trimmed, start_layer="receptors", stop_layer="genes",
        ) == []

        trimmed["tf_gene"] = trimmed["tf_gene"].iloc[0:0]
        assert cascade_core.empty_requested_cascade_edges(
            trimmed, start_layer="receptors", stop_layer="genes",
        ) == ["tf_gene"]

    def test_cascade_plot_warns_for_empty_requested_transition(
            self, monkeypatch, cascade_edges, cascade_results):
        empty_tf_gene = {key: value.copy() for key, value in cascade_edges.items()}
        empty_tf_gene["tf_gene"] = empty_tf_gene["tf_gene"].iloc[0:0]
        monkeypatch.setattr(
            cascade_plot_module, "build_networks",
            lambda *args, **kwargs: empty_tf_gene,
        )

        with pytest.warns(UserWarning, match="receiver TF → gene.*top_tf_n"):
            fig, ax = cascade_plot_fn(
                object(), cascade_results,
                cell_type="Receiver",
                start_layer="receptors", stop_layer="genes",
                show_labels=False,
            )
        plt.close(fig)

    def test_collect_node_sets_preserves_trimmed_boundary_ligands(self, cascade_edges):
        trimmed = cascade_core.trim_cascade_edges(cascade_edges, stop_layer="ligands")
        nodes = cascade_core.collect_node_sets(trimmed, seeds=None, cell_type="Receiver")

        assert nodes["ligands"] == ["LIG1::SenderA", "LIG2::SenderB"]
        assert nodes["recv_recs"] == []

    def test_cascade_plot_applies_layer_window_and_returns_edges(self, monkeypatch, cascade_edges, cascade_results):
        monkeypatch.setattr(cascade_plot_module, "build_networks", lambda *args, **kwargs: cascade_edges)

        fig, ax, edges = cascade_plot_fn(
            object(),
            cascade_results,
            cell_type="Receiver",
            stop_layer="ligands",
            show_labels=False,
            return_edges=True,
        )

        assert not edges["upstream_tf_lig"].empty
        assert edges["lig_rec"].empty
        plt.close(fig)


    def test_cascade_plot_warns_when_displayed_seeds_are_capped(self, monkeypatch, cascade_edges, cascade_results):
        monkeypatch.setattr(cascade_plot_module, "build_networks", lambda *args, **kwargs: cascade_edges)

        with pytest.warns(UserWarning, match="Showing only the first 1 seed genes"):
            fig, ax, edges = cascade_plot_fn(
                object(),
                cascade_results,
                cell_type="Receiver",
                seeds=["SEED1", "SEED2"],
                show_seeds=True,
                max_seed_nodes=1,
                show_labels=False,
                return_edges=True,
            )

        plt.close(fig)

    def test_contrast_cascade_plot_applies_layer_window_and_returns_edges(self, monkeypatch, cascade_edges, cascade_results):
        monkeypatch.setattr(cascade_plot_module, "build_networks_contrast", lambda *args, **kwargs: cascade_edges)
        result_b = cascade_results.copy()
        result_b["score"] = result_b["score"] * 0.5

        fig, ax, edges = contrast_cascade_plot_fn(
            {"a": object(), "b": object()},
            cascade_results,
            result_b,
            cell_type="Receiver",
            start_layer="ligands",
            show_labels=False,
            return_edges=True,
        )

        assert edges["upstream_r_tf"].empty
        assert not edges["lig_rec"].empty
        plt.close(fig)


class TestCascadeCoreHelpers:
    """Test pure helper behavior used by the cascade plots."""

    def test_node_type_color_overrides_are_partial_and_validated(self):
        palette = cascade_core.resolve_node_type_colors({
            "ligand": "#ff0000",
            "tf": (0.1, 0.2, 0.3),
        })

        assert palette["ligand"] == pytest.approx((1.0, 0.0, 0.0, 1.0))
        assert palette["tf"] == pytest.approx((0.1, 0.2, 0.3, 1.0))
        assert palette["receptor"] == cascade_core.TYPE_RGBA["receptor"]
        with pytest.raises(ValueError, match="Unknown node type color"):
            cascade_core.resolve_node_type_colors({"protein": "red"})
        with pytest.raises(ValueError, match="Invalid color"):
            cascade_core.resolve_node_type_colors({"gene": "not-a-color"})
        with pytest.raises(TypeError, match="must be a dictionary"):
            cascade_core.resolve_node_type_colors(["red"])

    def test_collect_node_sets_deduplicates_and_tracks_label_collisions(self, cascade_edges):
        nodes = cascade_core.collect_node_sets(cascade_edges, seeds=["SEED1"], cell_type="Receiver")

        assert nodes["sender_recs"] == ["SREC_receptor::SenderA"]
        assert nodes["sender_tfs"] == ["STF_TF::SenderA"]
        assert nodes["ligands"] == ["LIG1::SenderA", "LIG2::SenderB"]
        assert nodes["recv_recs"] == ["RREC_receptor::Receiver", "DUP_receptor::Receiver"]
        assert nodes["recv_tfs"] == ["RTF_TF::Receiver", "GENE_COLLIDE_TF::Receiver"]
        assert nodes["recv_genes"] == ["GENE1::Receiver", "GENE_COLLIDE::Receiver", "SEED2::Receiver"]
        assert nodes["seed_nodes"] == ["SEED1::Receiver"]
        assert nodes["name_collisions"] == {"GENE_COLLIDE"}

    def test_collect_node_sets_does_not_reintroduce_all_seed_targets_as_genes(
            self, cascade_edges):
        all_seed_edges = {key: value.copy() for key, value in cascade_edges.items()}
        all_seed_edges["tf_gene"] = pd.DataFrame({
            "source": ["RTF_TF::Receiver", "RTF_TF::Receiver"],
            "target": ["SEED1::Receiver", "SEED2::Receiver"],
            "value": [1.0, 0.8],
        })

        nodes = cascade_core.collect_node_sets(
            all_seed_edges, seeds=["SEED1", "SEED2"], cell_type="Receiver",
        )

        assert nodes["recv_genes"] == []
        assert nodes["seed_nodes"] == ["SEED1::Receiver", "SEED2::Receiver"]
        assert not (set(nodes["recv_genes"]) & set(nodes["seed_nodes"]))

    def test_compute_layout_positions_every_drawn_node_and_orders_receiver_receptors(self, cascade_edges):
        nodes = cascade_core.collect_node_sets(cascade_edges, seeds=["SEED1"], cell_type="Receiver")
        geo = cascade_core.compute_global_geometry(nodes)

        positions = cascade_core.compute_layout(
            nodes, cascade_edges, geo["sender_ct_order"], geo["sender_cy_by_ct"],
            geo["sender_cx"], geo["SENDER_R"], geo["extracell_x_lo"], geo["extracell_x_hi"],
            geo["recv_cx"], geo["recv_cy"], geo["RECV_R"], geo["NUC_CX"],
        )

        expected_nodes = (
            nodes["sender_recs"] + nodes["sender_tfs"] + nodes["ligands"]
            + nodes["recv_recs"] + nodes["recv_tfs"] + nodes["recv_genes"]
            + nodes["seed_nodes"]
        )
        assert set(expected_nodes) <= set(positions)
        assert nodes["recv_recs"] == ["RREC_receptor::Receiver", "DUP_receptor::Receiver"]

        nucleus_nodes = nodes["recv_genes"] + nodes["seed_nodes"]
        nucleus_a = 0.40 * geo["RECV_R"]
        nucleus_b = 0.40 * geo["RECV_R"] * geo["RECV_RY_RATIO"]
        for node in nucleus_nodes:
            x, y = positions[node]
            normalized_radius = np.sqrt(
                ((x - geo["NUC_CX"]) / nucleus_a) ** 2
                + ((y - geo["recv_cy"]) / nucleus_b) ** 2
            )
            assert normalized_radius <= cascade_core.NUCLEUS_LAYOUT_FRACTION

    def test_compute_node_radii_scales_by_group_and_can_be_disabled(self):
        nodes = ["low", "high", "missing"]
        radii = cascade_core.compute_node_radii(
            nodes, {"low": 1.0, "high": 9.0}, 1.0, True,
            group_lists=[["low", "high"], ["missing"]],
        )

        assert radii["high"] > radii["low"]
        assert radii["missing"] == pytest.approx(0.4)
        assert cascade_core.compute_node_radii(nodes, {}, 1.0, False, [nodes]) == {
            "low": 1.0, "high": 1.0, "missing": 1.0,
        }

    def test_color_helpers_handle_bounds_and_invalid_values(self):
        assert cascade_core.hex_to_rgb("#336699") == pytest.approx((0.2, 0.4, 0.6))
        assert cascade_core.lighten("#000000", 0.25) == pytest.approx((0.25, 0.25, 0.25))
        assert cascade_core.parse_rgba("rgba(255,128,0,0.5)") == pytest.approx((1.0, 128 / 255, 0.0, 0.5))
        assert cascade_core.score_to_rgba(np.nan, 1.0, (0.2, 0.4, 0.6)) == pytest.approx((1.0, 1.0, 1.0, 0.85))
        assert cascade_core.score_to_grey(1.0, 0.0, alpha=0.5) == pytest.approx((1.0, 1.0, 1.0, 0.5))
        assert cascade_core.score_to_grey(2.0, 1.0, alpha=0.5) == pytest.approx((0.45, 0.45, 0.45, 0.5))
        assert cascade_core.delta_to_rgba(np.nan) == pytest.approx((1, 1, 1, 0.5))

    def test_zscore_and_boost_preserve_empty_zero_and_sign_behavior(self):
        values = {"a": -2.0, "b": 4.0, "c": 0.0}
        cascade_core.zscore_layer(values, ["a", "b", "missing"])
        assert values["a"] == pytest.approx(-0.5)
        assert values["b"] == pytest.approx(1.0)
        assert values["c"] == 0.0

        boosted = cascade_core.boost_values_for_color({"a": -1.0, "b": 0.0, "c": 4.0}, 4.0, 0.25)
        assert boosted["a"] == pytest.approx(-1.75)
        assert boosted["b"] == 0.0
        assert boosted["c"] == pytest.approx(4.0)

    def test_draw_single_edge_filters_low_delta_and_uses_width_range(self, monkeypatch):
        added = []
        _, ax = plt.subplots()
        monkeypatch.setattr(ax, "add_patch", added.append)

        kwargs = dict(
            ax=ax, src="A", tgt="B", value=5.0, layer_key="rec_tf", rad=0.1,
            positions={"A": (0, 0), "B": (1, 0)},
            lw_ranges={"rec_tf": (0.0, 10.0)},
            color_fn=lambda _node, _layer: (0.2, 0.3, 0.4, 0.5),
            score_filter_dict={"A": 0.1},
            thresholds={"rec_tf": 0.2},
        )
        cascade_core.draw_single_edge(**kwargs)
        assert added == []

        kwargs["score_filter_dict"] = {"A": 0.3}
        cascade_core.draw_single_edge(**kwargs)
        assert len(added) == 2
        assert added[1].get_linewidth() == pytest.approx((cascade_core.EDGE_LW_MIN + cascade_core.EDGE_LW_MAX) / 2)
        plt.close("all")

    def test_geometry_helpers_cover_empty_single_and_grid_breaks(self):
        assert cascade_core.arc_positions(0, 0, 1, 0, 0, 90) == []
        assert cascade_core.arc_positions_ellipse(0, 0, 1, 2, 0, 0, 90) == []
        assert cascade_core.vertical_positions(1, 2, 0, 3) == []
        assert cascade_core.vertical_positions(1, 2, 1, 3) == [(1, 2)]
        assert cascade_core.vertical_positions(0, 0, 3, 1)[0] == pytest.approx((0, -1))
        assert cascade_core.ellipse_ring_positions(1, 2, 3, 4, 0) == []
        assert cascade_core.ellipse_ring_positions(1, 2, 3, 4, 1) == [(1, 2)]
        ring = cascade_core.ellipse_ring_positions(0, 0, 2, 1, 5)
        assert len(ring) == 5
        assert all(
            (x / 2) ** 2 + y ** 2 == pytest.approx(1.0)
            for x, y in ring
        )
        assert cascade_core.grid_positions_in_ellipse(0, 0, 1, 1, 0) == []
        assert len(cascade_core.grid_positions_in_ellipse(0, 0, 2, 1, 5)) == 5
        assert len(cascade_core.grid_positions_in_ellipse(0, 0, 2, 0, 3)) == 3
        seven_nodes = cascade_core.grid_positions_in_ellipse(2, 0, 1, 1, 7)
        assert seven_nodes[-1][0] == pytest.approx(2.0)

    def test_every_cascade_transition_has_a_drawing_configuration(self):
        assert list(cascade_core.CASCADE_EDGE_CURVATURE) == cascade_core.CASCADE_EDGE_ORDER

    def test_empty_collect_geometry_and_setup_figure_defaults(self):
        nodes = cascade_core.collect_node_sets(_empty_edges(), seeds=None, cell_type="Receiver")
        assert nodes["seed_nodes"] == []
        assert nodes["sender_recs"] == []

        geo = cascade_core.compute_global_geometry(nodes)
        assert geo["sender_ct_order"] == []
        assert geo["sender_cy_by_ct"] == {}
        assert geo["n_senders"] == 1
        assert geo["receiver_only"] is True
        assert geo["recv_cx"] == pytest.approx(0.0)
        assert geo["recv_cy"] == pytest.approx(0.0)

        fig, ax, y_max = cascade_core.setup_figure(geo)
        assert y_max == pytest.approx(
            geo["RECV_R"] * geo["RECV_RY_RATIO"] + 1.5
        )
        x_min, x_max = ax.get_xlim()
        assert (x_min + x_max) / 2 == pytest.approx(geo["recv_cx"])
        assert fig.get_size_inches()[0] < 10
        assert ax.axison is False
        plt.close(fig)

    def test_sender_receiver_colors_default_and_custom(self):
        sender_default = cascade_core.sender_colors(99, "SenderA", None)
        assert sender_default["edge"] == (0.35, 0.35, 0.35, 0.7)

        sender_custom = cascade_core.sender_colors(0, "SenderA", {"SenderA": "#336699"})
        assert sender_custom["label"] == pytest.approx((0.2, 0.4, 0.6))
        assert sender_custom["fill"][3] == pytest.approx(0.30)

        receiver_default = cascade_core.receiver_colors("Receiver", None)
        assert receiver_default["face"][3] == pytest.approx(0.30)

        receiver_custom = cascade_core.receiver_colors("Receiver", {"Receiver": "#336699"})
        assert receiver_custom["label"] == pytest.approx((0.2, 0.4, 0.6))
        assert receiver_custom["nuc_edge"][3] == pytest.approx(0.92)

    def test_edge_width_ranges_skip_empty_and_expand_constant_values(self):
        edges = _empty_edges()
        edges["rec_tf"] = pd.DataFrame({
            "source": ["A", "B"],
            "target": ["C", "D"],
            "value": [2.0, 2.0],
        })

        ranges = cascade_core.compute_edge_lw_ranges(edges)
        assert list(ranges) == ["rec_tf"]
        assert ranges["rec_tf"][0] == 2.0
        assert ranges["rec_tf"][1] > 2.0

    def test_zscore_layer_ignores_singletons_and_all_zero_values(self):
        one = {"a": 3.0}
        cascade_core.zscore_layer(one, ["a"])
        assert one == {"a": 3.0}

        zeros = {"a": 0.0, "b": 0.0}
        cascade_core.zscore_layer(zeros, ["a", "b"])
        assert zeros == {"a": 0.0, "b": 0.0}

    def test_node_radii_empty_and_default_for_ungrouped_nodes(self):
        assert cascade_core.compute_node_radii([], {}, 1.0, True, []) == {}
        radii = cascade_core.compute_node_radii(["grouped", "ungrouped"], {"grouped": 4.0}, 1.0, True, [["grouped"]])
        assert radii["ungrouped"] == 1.0

    def test_draw_helpers_create_expected_patches_and_text(self):
        _, ax = plt.subplots()

        cascade_core.draw_extracellular_region_per_ct(ax, 1, 2, [], {}, 1, 8, None)
        assert len(ax.patches) == 0

        cascade_core.draw_extracellular_region_per_ct(
            ax, 1, 2, ["SenderA", "SenderB"],
            {"SenderA": 0.0, "SenderB": 2.0}, 1.0, 8,
            {"SenderA": "#336699"},
        )
        boxes = [patch for patch in ax.patches if isinstance(patch, FancyBboxPatch)]
        assert len(boxes) == 2
        assert ax.texts[-1].get_text() == "extracellular"

        cascade_core.draw_sender_cell(
            ax, 0, 0, 1,
            {"fill": (1, 1, 1, 1), "edge": (0, 0, 0, 1), "arc": (0, 0, 0, 1), "label": (0, 0, 0)},
            "Sender label", 8,
        )
        cascade_core.draw_receiver_cell(
            ax, 4, 0, 1, 0.85, 4.1,
            {
                "face": (1, 1, 1, 1), "edge": (0, 0, 0, 1), "arc": (0, 0, 0, 1),
                "nuc_face": (1, 1, 1, 1), "nuc_edge": (0, 0, 0, 1), "label": (0, 0, 0),
            },
            "Receiver label", 8,
        )
        assert any(text.get_text() == "Sender label" for text in ax.texts)
        assert any(text.get_text() == "Receiver label" for text in ax.texts)
        plt.close("all")

    def test_draw_single_edge_and_node_early_returns_and_label_alignment(self, monkeypatch):
        added = []
        _, ax = plt.subplots()
        monkeypatch.setattr(ax, "add_patch", added.append)

        cascade_core.draw_single_edge(
            ax, "missing", "B", 1.0, "rec_tf", 0.1,
            {"B": (1, 0)}, {}, lambda _node, _layer: (1, 0, 0, 1),
        )
        assert added == []

        cascade_core.draw_single_edge(
            ax, "A", "B", 1.0, "unknown", 0.1,
            {"A": (0, 0), "B": (1, 0)}, {}, lambda _node, _layer: (1, 0, 0, 1),
        )
        assert len(added) == 2

        cascade_core.draw_single_node(
            ax, "missing", "gene", 0, 0, {}, {}, 0.2,
            lambda _node, _role: (1, 1, 1, 1), False, True, 10, set(),
        )
        assert len(added) == 2
        plt.close("all")

        _, ax = plt.subplots()
        for node, x, y, role in [
            ("LEFT_FAR::Receiver", -1.0, 1.0, "gene"),
            ("RIGHT::Receiver", 1.0, 0.0, "gene"),
            ("REC_receptor::Receiver", -1.0, 0.0, "recv_receptor"),
            ("LIG::Sender", 0.0, 0.0, "ligand"),
        ]:
            cascade_core.draw_single_node(
                ax, node, role, 0.0 if role != "ligand" else None, 0.0,
                {node: (x, y)}, {node: 0.2}, 0.2,
                lambda _node, _role: (1, 1, 1, 1), False, True, 10, set(),
            )

        assert [text.get_ha() for text in ax.texts] == ["right", "left", "left", "left"]
        assert ax.texts[2].get_fontsize() == pytest.approx(13)
        assert ax.texts[3].get_fontweight() == "bold"
        plt.close("all")

    def test_legends_include_expected_labels_and_normalization_note(self):
        _, ax = plt.subplots()
        cascade_core.draw_type_legend(ax, 10, show_seeds=False)
        assert [text.get_text() for text in ax.get_legend().texts] == ["Receptor", "TF", "Ligand", "Gene"]
        plt.close("all")

        _, ax = plt.subplots()
        cascade_core.draw_contrast_legend(ax, "temperature", 2.0, normalized=True, fontsize=10, show_seeds=True)
        labels = [text.get_text() for text in ax.get_legend().texts]
        assert "Seed" in labels
        assert "warm-enriched (hi)" in labels
        assert "receiver scores z-scored per layer" in [text.get_text() for text in ax.texts]
        plt.close("all")

        _, ax = plt.subplots()
        cascade_core.draw_contrast_legend(ax, "sex", 2.0, normalized=False, fontsize=10)
        labels = [text.get_text() for text in ax.get_legend().texts]
        assert "A-enriched (hi)" in labels
        assert len(ax.texts) == 0
        plt.close("all")

    def test_draw_cells_and_edges_orchestrates_edges_nodes_seed_label_and_section_labels(self, cascade_edges):
        nodes = cascade_core.collect_node_sets(cascade_edges, seeds=["SEED1"], cell_type="Receiver")
        geo = cascade_core.compute_global_geometry(nodes)
        fig, ax, y_max = cascade_core.setup_figure(geo, figsize=(8, 6))
        all_nodes = (
            nodes["sender_recs"] + nodes["sender_tfs"] + nodes["ligands"]
            + nodes["recv_recs"] + nodes["recv_tfs"] + nodes["recv_genes"]
            + nodes["seed_nodes"]
        )

        positions = cascade_core.draw_cells_and_edges(
            ax, nodes, cascade_edges, geo,
            cell_type="Receiver",
            celltype_colors={"SenderA": "#336699", "Receiver": "#663399"},
            celltype_display_names={"SenderA": "Sender A", "Receiver": "Receiver nice"},
            label_fontsize=7,
            node_radii={node: 0.15 for node in all_nodes},
            node_type_halo=True,
            show_labels=True,
            show_seeds=True,
            fill_fn=lambda _node, _role: (0.9, 0.9, 0.9, 1.0),
            edge_color_fn=lambda _node, _layer: (0.2, 0.2, 0.2, 1.0),
            seed_label="Seeds",
            seed_label_fontsize=6,
            score_filter_dict={node: 10.0 for node in all_nodes},
            lfc_thresholds={},
            flow="upstream",
        )
        cascade_core.add_section_labels(ax, geo, y_max, 7)

        assert set(all_nodes) <= set(positions)
        assert any(text.get_text() == "Seeds" for text in ax.texts)
        assert {"Sender cells", "Ligands", "Receiver"} <= {text.get_text() for text in ax.texts}
        expected_arrows = 2 * sum(
            (edges["source"].isin(positions) & edges["target"].isin(positions)).sum()
            for edges in cascade_edges.values()
        )
        assert len([patch for patch in ax.patches if isinstance(patch, FancyArrowPatch)]) == expected_arrows
        plt.close(fig)

    def test_build_networks_formats_maps_and_filters_sankey_layers(self, monkeypatch):
        nets = [
            pd.DataFrame({"receptor_clean": ["REC::Sender"], "tf_clean": ["TF::Sender"], "tf": ["TF_TF::Sender"], "weight": [1.0]}),
            pd.DataFrame({"tf_clean": ["TF::Sender"], "tf": ["TF_TF::Sender"], "gene": ["LIG::Sender"], "weight": [2.0]}),
            pd.DataFrame({"ligand": ["LIG::Sender"], "receptor_clean": ["RREC::Receiver"], "receptor": ["RREC_receptor::Receiver"], "weight": [3.0]}),
            pd.DataFrame({"receptor_clean": ["RREC::Receiver"], "tf_clean": ["RTF::Receiver"], "tf": ["RTF_TF::Receiver"], "weight": [4.0]}),
            pd.DataFrame({"tf_clean": ["RTF::Receiver"], "tf": ["RTF_TF::Receiver"], "gene": ["GENE::Receiver"], "weight": [5.0]}),
        ]
        calls = {}

        def fake_build_partial_networks(**kwargs):
            calls["kwargs"] = kwargs
            return nets

        fake_sankey_paths = types.SimpleNamespace(
            build_partial_networks=fake_build_partial_networks,
            _filter_connected_sankey_layers=lambda *dfs: dfs,
        )
        fake_plot = types.ModuleType("recon.plot")
        fake_plot.sankey_paths = fake_sankey_paths
        fake_recon = types.ModuleType("recon")
        monkeypatch.setitem(sys.modules, "recon", fake_recon)
        monkeypatch.setitem(sys.modules, "recon.plot", fake_plot)

        edges = cascade_core.build_networks(
            object(), pd.DataFrame({"node": [], "score": []}),
            cell_type="Receiver", seeds=None, ligand_cells=None,
            top_ligand_n=1, top_receptor_n=2, top_tf_n=3,
            before_top_n=4, per_celltype=False, verbose=True, flow="downstream",
        )

        assert calls["kwargs"]["seeds"] == []
        assert calls["kwargs"]["ligand_cells"] == []
        assert calls["kwargs"]["direction"] == "downstream"
        assert edges["upstream_r_tf"].to_dict("records") == [
            {"source": "REC::Sender", "target": "TF_TF::Sender", "value": 1.0}
        ]
        assert edges["rec_tf"].to_dict("records") == [
            {"source": "RREC::Receiver", "target": "RTF_TF::Receiver", "value": 4.0}
        ]
        assert edges["tf_gene"].loc[0, "source"] == "RTF_TF::Receiver"

    def test_build_networks_fallback_columns_and_contrast_score_replacement(self, monkeypatch):
        nets = [
            pd.DataFrame({"receptor": ["REC::Sender"], "tf_clean": ["TF::Sender"], "weight": [1.0]}),
            pd.DataFrame({"tf_clean": ["TF::Sender"], "gene": ["LIG::Sender"], "weight": [2.0]}),
            pd.DataFrame({"ligand": ["LIG::Sender"], "receptor_clean": ["RREC::Receiver"], "weight": [3.0]}),
            pd.DataFrame({"receptor": ["RREC::Receiver"], "tf": ["RTF::Receiver"], "weight": [4.0]}),
            pd.DataFrame({"tf_clean": ["RTF::Receiver"], "gene": ["GENE::Receiver"], "weight": [5.0]}),
        ]
        seen = {}

        def fake_build_networks(multicell_obj, results, **kwargs):
            seen["multicell_obj"] = multicell_obj
            seen["results"] = results.copy()
            seen["kwargs"] = kwargs
            return {"ok": True}

        fake_sankey_paths = types.SimpleNamespace(
            build_partial_networks=lambda **kwargs: nets,
            _filter_connected_sankey_layers=lambda *dfs: dfs,
        )
        fake_plot = types.ModuleType("recon.plot")
        fake_plot.sankey_paths = fake_sankey_paths
        monkeypatch.setitem(sys.modules, "recon", types.ModuleType("recon"))
        monkeypatch.setitem(sys.modules, "recon.plot", fake_plot)

        edges = cascade_core.build_networks(
            object(), pd.DataFrame({"node": [], "score": []}),
            cell_type="Receiver", seeds=["S"], ligand_cells=["L"],
            top_ligand_n=1, top_receptor_n=2, top_tf_n=3,
            before_top_n=4, per_celltype=True, verbose=False,
        )
        assert edges["upstream_r_tf"].to_dict("records") == [
            {"source": "REC::Sender", "target": "TF::Sender", "value": 1.0}
        ]
        assert edges["rec_tf"].to_dict("records") == [
            {"source": "RREC::Receiver", "target": "RTF::Receiver", "value": 4.0}
        ]

        monkeypatch.setattr(cascade_core, "build_networks", fake_build_networks)

        class MockMulti:
            def __init__(self, source):
                self.multiplexes = {
                    "cell_communication": {
                        "layers": [pd.DataFrame({
                            "source": [source, "DUP"],
                            "target": ["T", "T2"],
                            "weight": [1.0, 2.0],
                        })],
                    },
                }

        results = pd.DataFrame(
            {"node": ["LIG-CellA", "DIRECT::CellB"], "score": [99.0, 99.0]},
            index=[0, 0],
        )
        out = cascade_core.build_networks_contrast(
            {"A": MockMulti("DUP"), "B": MockMulti("UNIQ")},
            results,
            {"LIG::CellA": -2.0, "DIRECT::CellB": 3.0},
            cell_type="Receiver", seeds=["S"], ligand_cells=["L"],
            top_ligand_n=1, top_receptor_n=2, top_tf_n=3,
            before_top_n=4, per_celltype=True, verbose=False, flow="upstream",
        )

        assert out == {"ok": True}
        assert seen["results"]["score"].tolist() == [2.0, 3.0]
        assert seen["kwargs"]["flow"] == "upstream"
        merged_comm = seen["multicell_obj"].multiplexes["cell_communication"]["layers"][0]
        assert len(merged_comm) == 3

    def test_build_networks_oldest_upstream_column_fallback(self, monkeypatch):
        nets = [
            pd.DataFrame({"receptor": ["REC::Sender"], "tf": ["TF::Sender"], "weight": [1.0]}),
            pd.DataFrame({"tf_clean": ["TF::Sender"], "gene": ["LIG::Sender"], "weight": [2.0]}),
            pd.DataFrame({"ligand": ["LIG::Sender"], "receptor_clean": ["RREC::Receiver"], "weight": [3.0]}),
            pd.DataFrame({"receptor": ["RREC::Receiver"], "tf": ["RTF::Receiver"], "weight": [4.0]}),
            pd.DataFrame({"tf_clean": ["RTF::Receiver"], "gene": ["GENE::Receiver"], "weight": [5.0]}),
        ]

        fake_sankey_paths = types.SimpleNamespace(
            build_partial_networks=lambda **_kwargs: nets,
            _filter_connected_sankey_layers=lambda *dfs: dfs,
        )
        fake_plot = types.ModuleType("recon.plot")
        fake_plot.sankey_paths = fake_sankey_paths
        monkeypatch.setitem(sys.modules, "recon", types.ModuleType("recon"))
        monkeypatch.setitem(sys.modules, "recon.plot", fake_plot)

        edges = cascade_core.build_networks(
            object(), pd.DataFrame({"node": [], "score": []}),
            cell_type="Receiver", seeds=[], ligand_cells=[],
            top_ligand_n=1, top_receptor_n=1, top_tf_n=1,
            before_top_n=1, per_celltype=True, verbose=False,
        )

        assert edges["upstream_r_tf"].to_dict("records") == [
            {"source": "REC::Sender", "target": "TF::Sender", "value": 1.0}
        ]

    def test_draw_single_node_adds_halo_label_and_collision_suffix(self):
        _, ax = plt.subplots()
        cascade_core.draw_single_node(
            ax, "FOO_TF::Receiver", "tf", 0.0, 0.0,
            {"FOO_TF::Receiver": (-1.0, 0.0)},
            {"FOO_TF::Receiver": 0.3},
            0.2,
            lambda _node, _role: (1, 0, 0, 1),
            halo=True,
            show_labels=True,
            fontsize=10,
            name_collisions={"FOO"},
        )

        circles = [patch for patch in ax.patches if isinstance(patch, Circle)]
        assert len(circles) == 2
        assert ax.texts[0].get_text() == "FOO (TF)"
        plt.close("all")


class TestCascadePublicPlots:
    """Test public cascade plotting orchestration with deterministic networks."""

    def test_cascade_plot_wires_networks_nodes_labels_and_title(self, monkeypatch, cascade_edges, cascade_results):
        captured = {}

        def fake_draw(ax, nodes, edges, geo, **kwargs):
            captured.update(nodes=nodes, edges=edges, geo=geo, kwargs=kwargs)
            return {}

        monkeypatch.setattr(cascade_plot_module, "build_networks", lambda *args, **kwargs: cascade_edges)
        monkeypatch.setattr(cascade_plot_module, "draw_cells_and_edges", fake_draw)
        monkeypatch.setattr(plt, "show", lambda: None)

        fig, ax = cascade_plot_fn(
            object(), cascade_results, cell_type="Receiver", seeds=["SEED1"],
            show_seeds=True, flow="downstream", title=None, label_fontsize=8,
        )

        assert captured["nodes"]["seed_nodes"] == ["SEED1::Receiver"]
        assert captured["kwargs"]["show_seeds"] is True
        assert captured["kwargs"]["flow"] == "downstream"
        assert fig._suptitle is None
        assert ax.get_title() == "Cell signaling cascade (downstream) — Receiver"
        legend_labels = [text.get_text() for text in ax.get_legend().texts]
        assert "Seed" in legend_labels
        plt.close(fig)

    def test_cascade_plot_color_and_radius_options_reach_draw_layer(self, monkeypatch, cascade_edges, cascade_results):
        captured = {}

        def fake_draw(ax, nodes, edges, geo, **kwargs):
            captured.update(nodes=nodes, kwargs=kwargs)
            return {}

        monkeypatch.setattr(cascade_plot_module, "build_networks", lambda *args, **kwargs: cascade_edges)
        monkeypatch.setattr(cascade_plot_module, "draw_cells_and_edges", fake_draw)
        monkeypatch.setattr(plt, "show", lambda: None)

        fig, _ = cascade_plot_fn(
            object(), cascade_results, cell_type="Receiver", seeds=["SEED1"],
            show_seeds=False, node_size_by_weight=False,
            node_alpha=0.4, edge_alpha=0.3, label_fontsize=8,
        )

        radii = captured["kwargs"]["node_radii"]
        assert len(set(radii.values())) > 1  # TF/receptor scale-down is still applied.
        assert captured["kwargs"]["fill_fn"]("LIG1::SenderA", "ligand")[3] == pytest.approx(0.4)
        assert captured["kwargs"]["edge_color_fn"]("LIG1::SenderA", "lig_rec")[3] == pytest.approx(0.3)
        assert captured["nodes"]["seed_nodes"] == []
        plt.close(fig)

    def test_contrast_plot_builds_abs_delta_colors_thresholds_and_norm_title(
        self, monkeypatch, cascade_edges, cascade_results
    ):
        captured = {}

        def fake_draw(ax, nodes, edges, geo, **kwargs):
            captured.update(nodes=nodes, kwargs=kwargs)
            return {}

        results_b = cascade_results.copy()
        results_b["score"] = [0.5, 3.0, 1.5, 1.0, 5.0, 2.0, 2.0, 7.0, 0.0, 2.5, 1.5, 1.0, 1.2]

        monkeypatch.setattr(cascade_plot_module, "build_networks_contrast", lambda *args, **kwargs: cascade_edges)
        monkeypatch.setattr(cascade_plot_module, "draw_cells_and_edges", fake_draw)
        monkeypatch.setattr(plt, "show", lambda: None)

        fig, ax = contrast_cascade_plot_fn(
            {"A": object(), "B": object()}, cascade_results, results_b,
            cell_type="Receiver", seeds=["SEED1"], show_seeds=True,
            delta_vmax=4.0, delta_min_quantile=0.5,
            delta_min_color_fraction=0.25, normalize_receiver_scores=True,
            contrast_scheme="sex", flow="upstream", label_fontsize=8,
        )

        assert captured["kwargs"]["score_filter_dict"]["STF_TF::SenderA"] == pytest.approx(-1.0)
        assert "rec_tf" in captured["kwargs"]["lfc_thresholds"]
        assert captured["kwargs"]["fill_fn"]("LIG2::SenderB", "ligand")[0] < 1.0
        assert captured["kwargs"]["edge_color_fn"]("LIG1::SenderA", "lig_rec")[2] < 1.0
        assert "\u2020norm" in ax.get_title()
        assert "A (hi) vs B (lo)" in ax.get_title()
        plt.close(fig)

    def test_contrast_plot_custom_title_marks_normalization(self, monkeypatch, cascade_edges, cascade_results):
        monkeypatch.setattr(cascade_plot_module, "build_networks_contrast", lambda *args, **kwargs: cascade_edges)
        monkeypatch.setattr(cascade_plot_module, "draw_cells_and_edges", lambda *args, **kwargs: {})
        monkeypatch.setattr(plt, "show", lambda: None)

        fig, ax = contrast_cascade_plot_fn(
            {"A": object()}, cascade_results, cascade_results,
            cell_type="Receiver", title="My contrast", normalize_receiver_scores=True,
            label_fontsize=8,
        )

        assert ax.get_title() == "My contrast  \u2020norm"
        plt.close(fig)
