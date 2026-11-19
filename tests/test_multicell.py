"""Tests for recon.explore.Multicell class."""
import importlib
import pytest
import pandas as pd
import numpy as np
from recon.explore.recon import Celltype, Multicell

recon_module = importlib.import_module("recon.explore.recon")


class TestMulticellConstruction:
    """Test Multicell object creation and integration."""
    
    def test_create_multicell_from_celltypes(
        self, simple_grn, simple_receptor_grn, simple_cell_communication
    ):
        """Test Multicell creation from Celltype objects."""
        ct_a = Celltype(
            celltype_name="CellA",
            grn_graph=simple_grn.copy(),
            receptor_grn_bipartite=simple_receptor_grn.copy()
        )
        
        ct_b = Celltype(
            celltype_name="CellB",
            grn_graph=simple_grn.copy(),
            receptor_grn_bipartite=simple_receptor_grn.copy()
        )
        
        mc = Multicell(
            celltypes=[ct_a, ct_b],
            cell_communication_graph=simple_cell_communication
        )
        
        assert "cell_communication" in mc.multiplexes
        assert "CellA_grn" in mc.multiplexes
        assert "CellB_grn" in mc.multiplexes
        assert len(mc.celltypes_names) == 2
    
    def test_create_multicell_from_dicts(
        self, simple_grn, simple_receptor_grn, simple_cell_communication
    ):
        """Test Multicell creation from celltype dictionaries."""
        celltypes = [
            {
                "celltype_name": "CellA",
                "grn_graph": simple_grn.copy(),
                "receptor_grn_bipartite": simple_receptor_grn.copy()
            },
            {
                "celltype_name": "CellB",
                "grn_graph": simple_grn.copy(),
                "receptor_grn_bipartite": simple_receptor_grn.copy()
            }
        ]
        
        mc = Multicell(
            celltypes=celltypes,
            cell_communication_graph=simple_cell_communication
        )
        
        assert mc.celltypes_names == ["CellA", "CellB"]
    
    def test_node_suffix_addition(
        self, simple_grn, simple_receptor_grn, simple_cell_communication
    ):
        """Test that nodes get ::celltype suffix in Multicell."""
        ct_a = Celltype(
            celltype_name="CellA",
            grn_graph=simple_grn.copy(),
            receptor_grn_bipartite=simple_receptor_grn.copy()
        )
        
        ct_b = Celltype(
            celltype_name="CellB",
            grn_graph=simple_grn.copy(),
            receptor_grn_bipartite=simple_receptor_grn.copy()
        )
        
        mc = Multicell(
            celltypes=[ct_a, ct_b],
            cell_communication_graph=simple_cell_communication
        )
        
        # Check that GRN nodes have celltype suffix
        grn_layer = mc.multiplexes["CellA_grn"]["layers"][0]
        # At least one node should have ::CellA suffix
        assert any('::CellA' in str(node) for col in ['source', 'target'] 
                   for node in grn_layer[col].values if col in grn_layer.columns)
    
    def test_cell_communication_layer_format(
        self, simple_grn, simple_receptor_grn, simple_cell_communication
    ):
        """Test cell communication layer has correct format."""
        ct_a = Celltype(
            celltype_name="CellA",
            grn_graph=simple_grn.copy(),
            receptor_grn_bipartite=simple_receptor_grn.copy()
        )
        
        ct_b = Celltype(
            celltype_name="CellB",
            grn_graph=simple_grn.copy(),
            receptor_grn_bipartite=simple_receptor_grn.copy()
        )
        
        mc = Multicell(
            celltypes=[ct_a, ct_b],
            cell_communication_graph=simple_cell_communication
        )
        
        cc_layer = mc.multiplexes["cell_communication"]["layers"][0]
        
        # Should have ligand-celltype format: LIGAND1-CellA
        assert any('-' in str(node) for col in ['source', 'target'] 
                   for node in cc_layer[col].values if col in cc_layer.columns)
    
    def test_bipartite_ligand_connections(
        self, simple_grn, simple_receptor_grn, simple_cell_communication
    ):
        """Test bipartite connections from GRN to ligands."""
        ct_a = Celltype(
            celltype_name="CellA",
            grn_graph=simple_grn.copy(),
            receptor_grn_bipartite=simple_receptor_grn.copy()
        )
        
        mc = Multicell(
            celltypes=[ct_a],
            cell_communication_graph=simple_cell_communication
        )
        
        # Should have bipartite for ligands
        assert "CellA_to_ligands" in mc.bipartites
        ligand_bipartite = mc.bipartites["CellA_to_ligands"]
        assert ligand_bipartite["source"] == "cell_communication"
        assert ligand_bipartite["target"] == "CellA_grn"
    
    def test_bipartite_receptor_connections(
        self, simple_grn, simple_receptor_grn, simple_cell_communication
    ):
        """Test bipartite connections from receptors to cell communication."""
        ct_a = Celltype(
            celltype_name="CellA",
            grn_graph=simple_grn.copy(),
            receptor_grn_bipartite=simple_receptor_grn.copy()
        )
        
        mc = Multicell(
            celltypes=[ct_a],
            cell_communication_graph=simple_cell_communication
        )
        
        # Should have bipartite for receptors
        assert "CellA_to_receptors" in mc.bipartites
        receptor_bipartite = mc.bipartites["CellA_to_receptors"]
        assert receptor_bipartite["source"] == "cell_communication"
        assert receptor_bipartite["target"] == "CellA_receptor"

    @pytest.mark.parametrize(
        "kwargs, message",
        [
            ({"celltypes": "CellA"}, "celltypes must be a list or a dictionary"),
            ({"cell_communication_graph": "not a dataframe"}, "cell_communication_graph must be"),
            ({"bipartite_grn_cell_communication_directed": "bad"}, "bipartite_grn_cell_communication_directed"),
            ({"bipartite_cell_communication_receptor_directed": "bad"}, "bipartite_cell_communication_receptor_directed"),
            ({"bipartite_grn_cell_communication_weighted": "bad"}, "bipartite_grn_cell_communication_weighted"),
            ({"bipartite_cell_communication_receptor_weighted": "bad"}, "bipartite_cell_communication_receptor_weighted"),
            ({"cell_communication_graph_directed": "bad"}, "cell_communication_graph_directed"),
            ({"cell_communication_graph_weighted": "bad"}, "cell_communication_graph_weighted"),
        ],
    )
    def test_invalid_constructor_arguments_raise(
        self, kwargs, message, simple_grn, simple_receptor_grn, simple_cell_communication
    ):
        """Role: validate Multicell constructor rejects malformed inputs clearly."""
        ct = Celltype(
            celltype_name="CellA",
            grn_graph=simple_grn.copy(),
            receptor_grn_bipartite=simple_receptor_grn.copy(),
        )
        params = {
            "celltypes": [ct],
            "cell_communication_graph": simple_cell_communication.copy(),
        }
        params.update(kwargs)

        with pytest.raises(ValueError, match=message):
            Multicell(**params)

    def test_invalid_celltype_entry_raises(self, simple_cell_communication):
        """Role: document malformed list entries before constructor validation improves."""
        with pytest.raises(TypeError, match="object.*not subscriptable"):
            Multicell(
                celltypes=[object()],
                cell_communication_graph=simple_cell_communication.copy(),
            )


class TestMulticellLambMatrix:
    """Test transition matrix (lamb) for Multicell."""
    
    def test_default_lamb_shape(
        self, simple_grn, simple_receptor_grn, simple_cell_communication
    ):
        """Test that default lamb matrix has correct dimensions."""
        ct_a = Celltype(
            celltype_name="CellA",
            grn_graph=simple_grn.copy(),
            receptor_grn_bipartite=simple_receptor_grn.copy()
        )
        
        ct_b = Celltype(
            celltype_name="CellB",
            grn_graph=simple_grn.copy(),
            receptor_grn_bipartite=simple_receptor_grn.copy()
        )
        
        mc = Multicell(
            celltypes=[ct_a, ct_b],
            cell_communication_graph=simple_cell_communication
        )
        
        # Should have layers for: cell_comm, CellA_receptor, CellA_grn, CellB_receptor, CellB_grn
        expected_layers = 5
        assert mc.lamb.shape == (expected_layers, expected_layers)
    
    def test_lamb_allows_grn_transitions(
        self, simple_grn, simple_receptor_grn, simple_cell_communication
    ):
        """Test that lamb allows transitions within GRN layers."""
        ct_a = Celltype(
            celltype_name="CellA",
            grn_graph=simple_grn.copy(),
            receptor_grn_bipartite=simple_receptor_grn.copy()
        )
        
        mc = Multicell(
            celltypes=[ct_a],
            cell_communication_graph=simple_cell_communication
        )
        
        # GRN should allow self-transitions (normalized, so > 0)
        assert mc.lamb.loc["CellA_grn", "CellA_grn"] > 0


class TestMulticellRenaming:
    """Test celltype renaming in Multicell context."""
    
    def test_rename_via_dict(
        self, simple_grn, simple_receptor_grn, simple_cell_communication
    ):
        """Test renaming celltypes using dictionary keys."""
        ct_a = Celltype(
            celltype_name="OldA",
            grn_graph=simple_grn.copy(),
            receptor_grn_bipartite=simple_receptor_grn.copy()
        )
        
        ct_b = Celltype(
            celltype_name="OldB",
            grn_graph=simple_grn.copy(),
            receptor_grn_bipartite=simple_receptor_grn.copy()
        )
        
        mc = Multicell(
            celltypes={"NewA": ct_a, "NewB": ct_b},
            cell_communication_graph=simple_cell_communication
        )
        
        assert "NewA" in mc.celltypes_names
        assert "NewB" in mc.celltypes_names
        assert "OldA" not in mc.celltypes_names


class TestMulticellExplore:
    """Test Multicell.explore convenience method."""

    def test_explore_stores_multilayer_and_results(
        self, simple_grn, simple_receptor_grn, simple_cell_communication,
        monkeypatch
    ):
        """Test that Multicell.explore delegates to Multixrank then RWR."""
        expected = pd.DataFrame({"node": ["GENE1::CellA"], "score": [1.0]})

        class FakeMultilayer:
            def random_walk_rank(self, **kwargs):
                assert kwargs == {"tol": 1e-6}
                return expected

        def fake_multixrank(self, self_loops=True, restart_proba=0.7, verbose=True):
            assert self.seeds == ["GENE1::CellA"]
            assert self_loops is True
            assert restart_proba == 0.6
            assert verbose is False
            return FakeMultilayer()

        monkeypatch.setattr(Multicell, "Multixrank", fake_multixrank)

        ct_a = Celltype(
            celltype_name="CellA",
            grn_graph=simple_grn.copy(),
            receptor_grn_bipartite=simple_receptor_grn.copy()
        )

        mc = Multicell(
            celltypes=[ct_a],
            cell_communication_graph=simple_cell_communication
        )

        result = mc.explore(
            seeds=["GENE1::CellA"],
            restart_proba=0.6,
            verbose=False,
            tol=1e-6,
        )

        assert result is expected
        assert mc.results is expected
        assert hasattr(mc, "multilayer")


class TestMulticellIllustration:
    """Test plotting delegation methods on Multicell."""

    def test_illustrate_multicell_forwards_lamb_and_display_options(
        self, simple_grn, simple_receptor_grn, simple_cell_communication, monkeypatch
    ):
        """Role: verify Multicell.illustrate_multicell delegates to recon.plot correctly."""
        captured = {}

        def fake_illustrate_multicell(**kwargs):
            captured.update(kwargs)

        monkeypatch.setattr(recon_module.recon.plot, "illustrate_multicell", fake_illustrate_multicell)

        ct = Celltype(
            celltype_name="CellA",
            grn_graph=simple_grn.copy(),
            receptor_grn_bipartite=simple_receptor_grn.copy(),
        )
        mc = Multicell(
            celltypes=[ct],
            cell_communication_graph=simple_cell_communication.copy(),
        )

        mc.illustrate_multicell(
            figsize=(4, 3),
            azim=45,
            elev=10,
            display_layer_axis=False,
            display_self_proba=False,
            display_layer_names=True,
            alpha_layers=0.25,
            cell_communication_layer_name="ccc",
        )

        assert captured["lamb"] is mc.lamb
        assert captured["figsize"] == (4, 3)
        assert captured["azim"] == 45
        assert captured["elev"] == 10
        assert captured["display_layer_axis"] is False
        assert captured["display_self_proba"] is False
        assert captured["display_layer_names"] is True
        assert captured["alpha_layers"] == 0.25
        assert captured["cell_communication_layer_name"] == "ccc"


class TestMulticellOverlappingCelltypeNames:
    """Celltype names that are substrings of each other ('T' / 'T_reg') or of the
    ccc layer name ('cell' / 'cell_communication') must stay fully separate."""

    def test_overlapping_celltype_names_get_single_exact_suffix(
        self, simple_grn, simple_receptor_grn
    ):
        from recon.explore.recon import set_lambda

        names = ["T", "T_reg", "cell"]
        ccc = pd.DataFrame({
            "source": ["LIGAND1", "LIGAND2"],
            "target": ["RECEPTOR1", "RECEPTOR2"],
            "celltype_source": ["T", "T_reg"],
            "celltype_target": ["T_reg", "cell"],
            "lr_means": [0.5, 0.6],
        })
        multicell = Multicell(
            celltypes=[
                Celltype(
                    celltype_name=name,
                    grn_graph=simple_grn.copy(),
                    receptor_grn_bipartite=simple_receptor_grn.copy(),
                )
                for name in names
            ],
            cell_communication_graph=ccc,
        )

        # Explicit ownership records
        assert multicell.layer_celltypes == {
            f"{ct}_{kind}": ct for ct in names for kind in ("receptor", "grn")
        }
        assert multicell.bipartite_celltypes == {
            **{f"{ct}{suffix}": ct for ct in names
               for suffix in ("_to_ligands", "_to_receptors")},
            **{f"{ct}_grn-{ct}_receptor": ct for ct in names},
        }

        # Every node of a celltype layer carries exactly its own suffix
        for layer, ct in multicell.layer_celltypes.items():
            df = multicell.multiplexes[layer]["layers"][0]
            for node in pd.concat([df["source"], df["target"]]):
                assert node.count("::") == 1
                assert node.endswith(f"::{ct}")

        # The ccc layer itself is never suffixed
        ccc_layer = multicell.multiplexes["cell_communication"]["layers"][0]
        assert not pd.concat(
            [ccc_layer["source"], ccc_layer["target"]]).str.contains("::").any()

        # Bipartites: only the celltype-owned side is suffixed, exactly once
        for key, bip in multicell.bipartites.items():
            ct = multicell.bipartite_celltypes[key]
            edges = bip["edge_list_df"]
            if bip["target"] in multicell.layer_celltypes:
                assert edges["col1"].str.endswith(f"::{ct}").all()
                assert (edges["col1"].str.count("::") == 1).all()
            if bip["source"] in multicell.layer_celltypes:
                assert edges["col2"].str.endswith(f"::{ct}").all()
                assert (edges["col2"].str.count("::") == 1).all()
            else:
                assert not edges["col2"].str.contains("::").any()

        # Default and set_lambda transition matrices: no direct T <-> T_reg links
        t_layers = ["T_receptor", "T_grn"]
        treg_layers = ["T_reg_receptor", "T_reg_grn"]
        for lamb in (multicell.lamb, set_lambda(multicell=multicell)):
            assert (lamb.loc[t_layers, treg_layers] == 0).all().all()
            assert (lamb.loc[treg_layers, t_layers] == 0).all().all()
