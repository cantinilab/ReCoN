"""Tests for recon.explore helper functions."""
import importlib
import warnings
import pytest
import pandas as pd
import numpy as np

from recon.explore.recon import (
    combine_effects,
    format_multicell_results,
    multicell_targets,
    set_lambda,
    summarize_indirect_effects,
)

recon_module = importlib.import_module("recon.explore.recon")


def test_format_multicell_results_filters_layers_and_pivots_profiles():
    """Role: verify Multicell RWR results become a gene x celltype score matrix."""
    results = pd.DataFrame({
        "node": ["G1::CellA", "G1::CellB", "G2::CellA", "LIG-CellA", "G3::Other"],
        "layer": ["gene", "gene", "gene", "cell_communication", "gene"],
        "score": [0.1, 0.2, 0.3, 0.9, 1.0],
    })

    formatted = format_multicell_results(
        results,
        celltypes=["CellA", "CellB"],
        keep_layers="gene",
    )

    assert formatted.index.tolist() == ["G1", "G2"]
    assert formatted.columns.tolist() == ["CellA", "CellB"]
    assert formatted.loc["G1", "CellA"] == pytest.approx(0.1)
    assert formatted.loc["G1", "CellB"] == pytest.approx(0.2)
    assert pd.isna(formatted.loc["G2", "CellB"])


def test_format_multicell_results_does_not_warn_on_slice_assignment():
    """Role: ensure formatting works on a copy instead of a pandas slice view."""
    results = pd.DataFrame({
        "node": ["G1::CellA", "G1::CellB", "LIG-CellA"],
        "layer": ["gene", "gene", "cell_communication"],
        "score": [0.1, 0.2, 0.9],
    })

    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        formatted = format_multicell_results(results, celltypes=["CellA"])

    assert not any(isinstance(w.message, pd.errors.SettingWithCopyWarning) for w in caught)
    assert formatted.loc["G1", "CellA"] == pytest.approx(0.1)


def test_set_lambda_celltypes_only_path_builds_expected_layers():
    """Role: the celltypes-only path builds one receptor/grn pair per celltype plus ccc."""
    lamb = set_lambda(multicell=None, celltypes=["CellA", "CellB"])

    expected = ["CellA_receptor", "CellA_grn", "CellB_receptor", "CellB_grn",
                "cell_communication"]
    assert list(lamb.index) == expected
    assert list(lamb.columns) == expected
    assert (lamb.sum(axis=1) - 1.0).abs().max() < 1e-10


@pytest.mark.parametrize("strategy", ["intracell", "intercell"])
@pytest.mark.parametrize("direction", ["downstream", "upstream"])
def test_set_lambda_keeps_overlapping_celltype_names_separate(direction, strategy):
    """Role: 'T' is a substring of 'T_reg'; their layers must never be linked directly."""
    lamb = set_lambda(celltypes=["T", "T_reg"], direction=direction, strategy=strategy)

    t_layers = ["T_receptor", "T_grn"]
    treg_layers = ["T_reg_receptor", "T_reg_grn"]
    assert (lamb.loc[t_layers, treg_layers] == 0).all().all()
    assert (lamb.loc[treg_layers, t_layers] == 0).all().all()

    # Same matrix as two non-overlapping celltypes
    reference = set_lambda(celltypes=["A", "B"], direction=direction, strategy=strategy)
    assert lamb.to_numpy().tolist() == reference.to_numpy().tolist()


def test_rename_celltype_replaces_exact_layer_names_only(simple_grn, simple_receptor_grn):
    """Role: renaming 're' must not touch the 're' inside '_receptor'."""
    from recon.explore.recon import Celltype

    ct = Celltype(
        celltype_name="re",
        grn_graph=simple_grn.copy(),
        receptor_grn_bipartite=simple_receptor_grn.copy(),
    )
    ct.rename_celltype("Mono")

    assert set(ct.multiplexes) == {"Mono_receptor", "Mono_grn"}
    assert set(ct.bipartites) == {"Mono_grn-Mono_receptor"}
    bip = ct.bipartites["Mono_grn-Mono_receptor"]
    assert (bip["source"], bip["target"]) == ("Mono_receptor", "Mono_grn")
    assert set(ct.lamb.index) == set(ct.lamb.columns) == {"Mono_receptor", "Mono_grn"}
    assert set(ct.eta.index) == {"Mono_receptor", "Mono_grn"}
    assert ct.celltype_name == "Mono"


def test_set_lambda_warns_when_both_multicell_and_celltypes_are_provided(
    simple_grn, simple_receptor_grn, simple_cell_communication
):
    """Role: verify the multicell object takes precedence when both inputs are passed."""
    from recon.explore.recon import Celltype, Multicell

    ct = Celltype(
        celltype_name="CellA",
        grn_graph=simple_grn.copy(),
        receptor_grn_bipartite=simple_receptor_grn.copy(),
    )
    multicell = Multicell(
        celltypes=[ct],
        cell_communication_graph=simple_cell_communication.copy(),
    )

    with pytest.warns(UserWarning, match="multicell will be used"):
        lamb = set_lambda(multicell=multicell, celltypes=["Ignored"])

    assert list(lamb.index) == list(multicell.multiplexes.keys())


def test_set_lambda_intracell_and_intercell_transition_roles(
    simple_grn, simple_receptor_grn, simple_cell_communication
):
    """Role: lock down strategy/direction semantics for transition matrices."""
    from recon.explore.recon import Celltype, Multicell

    ct = Celltype(
        celltype_name="CellA",
        grn_graph=simple_grn.copy(),
        receptor_grn_bipartite=simple_receptor_grn.copy(),
    )
    multicell = Multicell(
        celltypes=[ct],
        cell_communication_graph=simple_cell_communication.copy(),
    )

    intracell = set_lambda(multicell=multicell, direction="downstream", strategy="intracell")
    intercell = set_lambda(multicell=multicell, direction="downstream", strategy="intercell")
    upstream = set_lambda(multicell=multicell, direction="upstream", strategy="intercell")

    assert intracell.loc["CellA_grn", "cell_communication"] == 0
    assert intercell.loc["CellA_grn", "cell_communication"] > 0
    assert upstream.loc["cell_communication", "CellA_grn"] > 0
    assert (intercell.sum(axis=1) - 1.0).abs().max() < 1e-10


@pytest.mark.parametrize(
    "kwargs, message",
    [
        ({"direction": "sideways"}, "direction must be either"),
        ({"strategy": "extracell"}, "strategy must be either"),
    ],
)
def test_set_lambda_rejects_unknown_direction_or_strategy(kwargs, message):
    """Role: ensure transition matrices only accept documented modes."""
    with pytest.raises(ValueError, match=message):
        set_lambda(multicell=None, celltypes=None, **kwargs)


def test_summarize_indirect_effects_splits_direct_and_indirect_tables():
    """Role: verify combined contribution tables are split into public outputs."""
    predictions = pd.DataFrame(
        {
            ("CellA", "CellA_direct"): [0.7, 0.3],
            ("CellA", "CellA"): [0.1, 0.2],
            ("CellA", "CellB"): [0.4, 0.5],
            ("CellB", "CellB_direct"): [0.6, 0.4],
            ("CellB", "CellA"): [0.2, 0.1],
            ("CellB", "CellB"): [0.3, 0.2],
        },
        index=["G1", "G2"],
    )
    predictions.columns = pd.MultiIndex.from_tuples(predictions.columns)

    direct, indirect = summarize_indirect_effects(predictions)

    assert direct.index.name == "gene"
    assert direct.columns.name == "celltype_target"
    assert indirect.columns.names == ["celltype_target", "celltype_source"]
    assert direct.loc["G1", "CellA"] == pytest.approx(0.7)
    assert indirect.loc["G1", ("CellA", "CellB")] == pytest.approx(0.4)


def test_combine_effects_without_cell_comm_matrix_normalizes_and_blends():
    """Role: verify the default direct/indirect blending path."""
    direct = pd.DataFrame({"CellA": [2.0, 1.0], "CellB": [1.0, 3.0]}, index=["G1", "G2"])
    indirect = pd.DataFrame(
        {
            ("CellA", "CellA"): [1.0, 1.0],
            ("CellA", "CellB"): [3.0, 1.0],
            ("CellB", "CellA"): [1.0, 3.0],
            ("CellB", "CellB"): [1.0, 1.0],
        },
        index=["G1", "G2"],
    )
    indirect.columns = pd.MultiIndex.from_tuples(indirect.columns)

    combined = combine_effects(direct.copy(), indirect.copy(), alpha=0.5)

    assert combined.columns.tolist() == ["CellA", "CellB"]
    assert combined.loc["G1", "CellA"] == pytest.approx(2 / 3)
    assert combined.loc["G2", "CellB"] == pytest.approx(17 / 24)


def test_combine_effects_with_cell_comm_matrix_requires_known_celltypes():
    """Role: ensure weighted communication matrices fail loudly for missing targets."""
    direct = pd.DataFrame({"CellA": [1.0]}, index=["G1"])
    indirect = pd.DataFrame({("CellA", "CellA"): [1.0]}, index=["G1"])
    indirect.columns = pd.MultiIndex.from_tuples(indirect.columns)
    weights = pd.DataFrame({"CellA": [1.0]}, index=["Other"])

    with pytest.raises(ValueError, match="Cell type CellA not found"):
        combine_effects(direct, indirect, cell_comm_matrix=weights)


def test_combine_effects_with_cell_comm_matrix_uses_weighted_sources():
    """Role: verify source-specific communication weights contribute to the blend."""
    direct = pd.DataFrame({"CellA": [2.0, 1.0], "CellB": [1.0, 3.0]}, index=["G1", "G2"])
    indirect = pd.DataFrame(
        {
            ("CellA", "CellA"): [1.0, 1.0],
            ("CellA", "CellB"): [3.0, 1.0],
            ("CellB", "CellA"): [1.0, 3.0],
            ("CellB", "CellB"): [1.0, 1.0],
        },
        index=["G1", "G2"],
    )
    indirect.columns = pd.MultiIndex.from_tuples(indirect.columns)
    weights = pd.DataFrame(
        {"CellA": [0.25, 0.75], "CellB": [0.75, 0.25]},
        index=["CellA", "CellB"],
    )

    combined = combine_effects(direct.copy(), indirect.copy(), alpha=0.5, cell_comm_matrix=weights)

    assert combined.columns.tolist() == ["CellA", "CellB"]
    assert combined.notna().all().all()
    assert combined.loc["G1", "CellA"] > combined.loc["G2", "CellA"]


def test_multicell_targets_requires_non_empty_celltypes(simple_grn, simple_receptor_grn, simple_cell_communication):
    """Role: fail early when no target cell types are available."""
    with pytest.raises(ValueError, match="celltypes should be a non-empty"):
        multicell_targets(
            seeds=["LIG"],
            celltypes=[],
            ccc=simple_cell_communication.copy(),
            grn=simple_grn.copy(),
            receptor_grn=simple_receptor_grn.copy(),
            verbose=False,
        )


def test_multicell_targets_rejects_unknown_receptor_grn_path(simple_grn, simple_cell_communication):
    """Role: report invalid receptor prior resources/paths with a user-facing error."""
    with pytest.raises(ValueError, match="receptor_grn should be a valid resource"):
        multicell_targets(
            seeds=["LIG"],
            celltypes=["CellA"],
            ccc=simple_cell_communication.copy(),
            grn=simple_grn.copy(),
            receptor_grn="/definitely/not/a/receptor_prior.csv",
            verbose=False,
        )


def test_multicell_targets_orchestrates_direct_and_indirect_effects_with_mocked_rwr(
    simple_grn, simple_receptor_grn, simple_cell_communication, monkeypatch
):
    """Role: exercise multicell_targets control flow without running HuMMuS."""
    calls = []

    class FakeMultilayer:
        def __init__(self, multicell):
            self.multicell = multicell

        def random_walk_rank(self):
            calls.append(dict(seeds=self.multicell.seeds, lamb=self.multicell.lamb.copy()))
            if isinstance(self.multicell.seeds, dict):
                return pd.DataFrame({
                    "node": ["TARGET1::CellA", "TARGET2::CellA"],
                    "layer": ["gene", "gene"],
                    "multiplex": ["CellA_grn", "CellA_grn"],
                    "score": [0.5, 0.25],
                })
            return pd.DataFrame({
                "node": ["RECEPTOR2::CellA", "TARGET1::CellA", "TARGET2::CellA"],
                "layer": ["gene", "gene", "gene"],
                "multiplex": ["CellA_grn", "CellA_grn", "CellA_grn"],
                "score": [2.0, 1.0, 0.5],
            })

    def fake_multixrank(self, restart_proba=0.6, verbose=False, **kwargs):
        calls.append(dict(restart_proba=restart_proba, verbose=verbose))
        return FakeMultilayer(self)

    monkeypatch.setattr(recon_module.Multicell, "Multixrank", fake_multixrank)

    direct, indirect = multicell_targets(
        seeds=["LIGAND1"],
        celltypes=["CellA"],
        ccc=simple_cell_communication.copy(),
        grn=simple_grn.copy(),
        receptor_grn=simple_receptor_grn.copy(),
        restart_proba=0.4,
        extend_seeds=True,
        celltype_to_ccc_proba=pd.Series({"CellA": 0.75}),
        njobs=1,
        verbose=False,
    )

    assert calls[0]["restart_proba"] == 0.4
    assert calls[1]["seeds"] == ["LIGAND1-CellA"]
    assert calls[3]["seeds"] == {"RECEPTOR2-CellA": 2.0}
    assert direct.index.name == "gene"
    assert indirect.columns.names == ["celltype_target", "celltype_source"]
    assert direct.loc["RECEPTOR2", "CellA"] == pytest.approx(2.0)
    assert indirect.loc["TARGET1", ("CellA", "CellA")] == pytest.approx(0.75)


def test_multicell_targets_accepts_ccc_to_celltype_proba_after_celltype_conversion(
    simple_grn, simple_receptor_grn, simple_cell_communication, monkeypatch
):
    """Role: verify ccc_to_celltype_proba is indexed by celltype names, not the celltype dict."""
    monkeypatch.setattr(
        recon_module.Multicell,
        "Multixrank",
        lambda self, **_kwargs: type(
            "FakeMultilayer",
            (),
            {"random_walk_rank": lambda _self: pd.DataFrame({
                "node": ["RECEPTOR2::CellA"],
                "layer": ["gene"],
                "multiplex": ["CellA_grn"],
                "score": [1.0],
            })},
        )(),
    )

    multicell_targets(
        seeds=["LIGAND1"],
        celltypes=["CellA"],
        ccc=simple_cell_communication.copy(),
        grn=simple_grn.copy(),
        receptor_grn=simple_receptor_grn.copy(),
        ccc_to_celltype_proba=pd.Series({"CellA": 1.0}),
        njobs=1,
        verbose=False,
    )


def test_combine_effects_zero_sum_columns_currently_return_nan():
    """Role: document current zero-sum behavior so future handling can be explicit."""
    direct = pd.DataFrame({"CellA": [0.0, 0.0]}, index=["G1", "G2"])
    indirect = pd.DataFrame({("CellA", "CellA"): [0.0, 0.0]}, index=["G1", "G2"])
    indirect.columns = pd.MultiIndex.from_tuples(indirect.columns)

    combined = combine_effects(direct.copy(), indirect.copy(), alpha=0.5)

    assert np.isnan(combined.loc["G1", "CellA"])
