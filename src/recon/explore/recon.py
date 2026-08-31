import hummuspy
from joblib import Parallel, delayed

import tqdm
import numpy as np
import pandas as pd

from typing import Union, List, Tuple, Dict, Any
import warnings

import copy
import recon.plot


class Celltype:
    """Single-cell-type multilayer network for ReCoN exploration.

    A ``Celltype`` contains two multiplexes, one for intracellular gene
    regulation and one for receptors, plus a bipartite receptor-to-gene layer.
    Use :meth:`Multixrank` when you need direct access to the underlying
    HuMMuS/Multixrank object, or :meth:`explore` for the common one-line
    workflow that builds the multilayer object and runs random walk with
    restart.
    """

    def __init__(
        self,
        celltype_name: str,
        grn_graph: Union[str, pd.DataFrame],
        receptor_grn_bipartite: Union[str, pd.DataFrame],
        receptor_graph: Union[str, pd.DataFrame, None] = None,
        lamb=None,
        eta=None,
        receptor_graph_directed: bool = True,
        receptor_graph_weighted: bool = False,
        grn_graph_directed: bool = False,
        grn_graph_weighted: bool = True,
        receptor_grn_bipartite_graph_directed: bool = False,
        receptor_grn_bipartite_graph_weighted: bool = True,
        copy_graphs=True,
        seeds: Union[List, Dict] = []
    ):
        """
        Create a Celltype object to explore the multilayer graph of a celltype.
        It can be used to predict the targets/regulators of a celltype based on
        the expression of its genes and the ligands it produces.

        It can be used in a Multicell object to explore the interactions between
        different celltypes.

        Parameters
        ----------
        receptor_graph : Union[str, pd.DataFrame, None], default None
            The graph of the receptors of the celltype.
            It can be a string with the path to a file containing the graph or a pandas DataFrame.
            If None, a fake receptor graph will be created with a fake receptor connecting all nodes
        grn_graph : Union[str, pd.DataFrame]
            The graph of the genes of the celltype.
            It can be a string with the path to a file containing the graph or a pandas DataFrame.
        receptor_grn_bipartite : Union[str, pd.DataFrame]
            The receptor_grn_bipartite graph between the receptors and the genes of the celltype.
            It can be a string with the path to a file containing the graph or a pandas DataFrame.
        celltype_name : str
            The name of the celltype.
        lamb : pd.DataFrame, optional
            The transition matrix between the layers of the celltype.
            If None, it will be set to a default one allowing transition between receptor,
            and receptor and grn layer exploration 
        eta : pd.Series, optional
            The probability of restarting at each layer.
            If None, it will be set to a vector with the same probability for each layer.
        receptor_graph_directed : bool, optional
            If True, the receptor graph is directed.
            The default is True.
        receptor_graph_weighted : bool, optional
            If True, the receptor graph is weighted.
            The default is False.
        grn_graph_directed : bool, optional
            If True, the grn graph is directed.
            The default is False.
        grn_graph_weighted : bool, optional 
            If True, the grn graph is weighted.
            The default is True.
        receptor_grn_bipartite_graph_directed : bool, optional
            If True, the receptor_grn_bipartite graph is directed.
            The default is True.
        receptor_grn_bipartite_graph_weighted : bool, optional
            If True, the receptor_grn_bipartite graph is weighted.
            The default is False.
        copy_graphs : bool, optional
            If True, the graphs are copied.
            The default is True.
        seeds : Union[list, dict], optional
            The seeds to explore the multilayer graph.
            If a list, the seeds are the nodes to explore.
            If a dictionary, the seeds are the nodes with their weights.
            The default is [].

        Returns
        -------
        Celltype
        """

        self.seeds = seeds
        self.celltype_name = celltype_name

        # Copy networks, True if no problem of memory
        if copy_graphs:
            if isinstance(receptor_graph, pd.DataFrame):
                receptor_graph = receptor_graph.copy()
            if isinstance(grn_graph, pd.DataFrame):
                grn_graph = grn_graph.copy()
            if isinstance(receptor_grn_bipartite, pd.DataFrame):
                receptor_grn_bipartite = receptor_grn_bipartite.copy()

        # Type of graphs
        receptor_graph_type = \
            str(int(receptor_graph_directed)) + \
            str(int(receptor_graph_weighted))
        grn_graph_type = \
            str(int(grn_graph_directed)) + \
            str(int(grn_graph_weighted))
        receptor_grn_bipartite_graph_type = \
            str(int(receptor_grn_bipartite_graph_directed)) + \
            str(int(receptor_grn_bipartite_graph_weighted))

        # In most cases, there is no receptor_graph.
        # The receptor_graph is then created with a fake receptor
        # that is not used in the analysis but essential since
        # all layers should have links
        if receptor_graph is None:
            warnings.warn(
                """
                No receptor_graph provided,
                an empty receptor graph will be created.
                """)
            if isinstance(receptor_grn_bipartite, str):
                warnings.warn(
                    """
                    The receptor_grn_bipartite is provided as a string,
                    it will be read now as a csv file to define
                    the fake receptor graph.
                    """)

                receptor_grn_bipartite = pd.read_csv(
                    receptor_grn_bipartite, sep=None, engine='python')
            
            receptor_graph = pd.DataFrame({
                "source": receptor_grn_bipartite.loc[:, 'source'].unique(),
                "target": ["fake" for r in range(
                    len(receptor_grn_bipartite.loc[:, 'source'].unique()))]
                }
            )
        
        # Add "_receptor" suffix to receptor nodes
        receptor_graph['source'] = receptor_graph['source'] + "_receptor"
        receptor_graph['target'] = receptor_graph['target'] + "_receptor"
        receptor_grn_bipartite['source'] = \
            receptor_grn_bipartite['source'] + "_receptor"
        
        # Format multiplex dictionary
        self.multiplexes = {
            celltype_name + "_receptor": {
                "names": ["receptor"],
                "graph_type": [receptor_graph_type],
                "layers": [receptor_graph]
            },
            celltype_name + "_grn": {
                "names": ["gene"],
                "graph_type": [grn_graph_type],
                "layers": [grn_graph]
            }
        }

        # Format receptor_grn_bipartite dictionary
        receptor_grn_bipartite = receptor_grn_bipartite.rename({
            'source': 'col2',
            'target': 'col1',
            'score': 'weight'}, axis=1)

        self.bipartites = {
            celltype_name+"_grn-" + celltype_name + "_receptor": {
                "source": celltype_name + "_receptor",
                "target": celltype_name + "_grn",
                "graph_type": receptor_grn_bipartite_graph_type,
                "edge_list_df": receptor_grn_bipartite
            }
        }

        # Prepare lamb
        if lamb is None:
            lamb = pd.DataFrame(
                np.zeros((len(self.multiplexes), len(self.multiplexes))),
                index=list(self.multiplexes.keys()),
                columns=list(self.multiplexes.keys())
            )
            lamb.loc[
                self.celltype_name + "_grn",
                self.celltype_name + "_grn"] = 1
            lamb.loc[
                self.celltype_name + "_receptor",
                self.celltype_name + "_grn"] = 1
        self.lamb = lamb

        # Fake eta
        if eta is None:
            self.eta = pd.Series(np.ones(2)/2, index=self.multiplexes.keys())
        else:
            self.eta = eta

    def Multixrank(
        self,
        self_loops=True,
        restart_proba=0.7,
        verbose=True
    ):
        """Build the HuMMuS ``Multixrank`` object for this cell type.

        Parameters
        ----------
        self_loops : bool, default=True
            Whether to add self-loops in each layer.
        restart_proba : float, default=0.7
            Random-walk restart probability.
        verbose : bool, default=True
            Passed to HuMMuS when constructing the multilayer object.

        Returns
        -------
        hummuspy.create_multilayer.Multixrank
            Configured multilayer object. Call ``random_walk_rank()`` on it to
            compute node scores, or use :meth:`explore` to do both steps.
        """

        if self.seeds == []:
            warnings.warn("""
            No seeds provided to explore the multilayer, all the scoree will be null (np.nan).
            You can pass it as an argument earlier or set up the .seeds attribute.
            """)
        
        # Swap source/target for multixrank compatibility (temporary copy)
        # This keeps internal representation biological but passes multixrank the format it expects
        multiplexes_for_multixrank = {}
        for key, value in self.multiplexes.items():
            # Create copy with swapped columns
            layer_copy = value["layers"][0].rename(
                columns={"source": "target", "target": "source"}, inplace=False)
            # Add network_key column (needed by hummuspy)
            layer_copy["network_key"] = key
            
            multiplexes_for_multixrank[key] = {
                "names": value["names"],
                "graph_type": value["graph_type"],
                "layers": [layer_copy]
            }
        
        # if seeds are provided as a dictionary
        if isinstance(self.seeds, dict):
            # Every value should be a numeric value (int or float)
            print("Seeds are provided as a dictionary with weights per seed.")
            if not all(
                [isinstance(v, int) or isinstance(v, float)
                    for v in self.seeds.values()]):
                raise ValueError("Seeds should be a dictionary with numeric values (weight per seed).")

            print("Creating a multixrank object with seeds as a dictionary.")
            # Create the multixrank object
            multilayer = hummuspy.create_multilayer.Multixrank(
                multiplex=multiplexes_for_multixrank,
                bipartite=self.bipartites,
                lamb=self.lamb.values.T,
                seeds=[],
                self_loops=self_loops,
                restart_proba=restart_proba,
                eta=self.eta.tolist(),
                verbose=verbose)

            print("Identifying produced ligands in response to the perturbation.")
            # Create a probability vector from the seeds' names and weights
            node_list = [node for node_list
                         in multilayer.multiplexall_node_list2d
                         for node in node_list]
            prox_vector = np.zeros(len(node_list))

            # Position is matching node order in the multilayer.
            node_arr = np.array(node_list)
            for seed, value in self.seeds.items():
                idx = np.where(node_arr == seed)[0]
                if idx.size:
                    prox_vector[idx[0]] = value
            # Values are normalized.
            multilayer.pr = prox_vector/prox_vector.sum()
#            multilayer.seed = hummuspy.create_multilayer.Seed(
#                "",
#                seeds=new_seeds,
#                multiplexall = multilayer.multiplexall_obj)

        else:
            #  should be a list or an array
            if not isinstance(self.seeds, (list, np.ndarray)):
                raise ValueError("Seeds should be a list or an array.")

            # Create the multixrank object
            multilayer = hummuspy.create_multilayer.Multixrank(
                multiplex=multiplexes_for_multixrank,
                bipartite=self.bipartites,
                lamb=self.lamb.values.T,
                seeds=self.seeds,
                self_loops=self_loops,
                restart_proba=restart_proba,
                eta=self.eta.tolist(),
                verbose=verbose)

        return multilayer

    def explore(
        self,
        seeds: Union[List, Dict, None] = None,
        self_loops=True,
        restart_proba=0.7,
        verbose=True,
        **random_walk_rank_kwargs
    ):
        """Build the multilayer object and run random walk with restart.

        This is the convenience method for the usual workflow::

            results = celltype.explore(seeds=["MYC"], restart_proba=0.7)

        It stores the created HuMMuS object in ``self.multilayer`` and the
        resulting scores in ``self.results`` before returning the scores.

        Parameters
        ----------
        seeds : list, dict, optional
            Seeds to use for this run. If provided, ``self.seeds`` is updated
            before building the multilayer object. Dictionaries are interpreted
            as weighted seeds.
        self_loops : bool, default=True
            Whether to add self-loops in each layer.
        restart_proba : float, default=0.7
            Random-walk restart probability.
        verbose : bool, default=True
            Passed to :meth:`Multixrank`.
        **random_walk_rank_kwargs
            Extra keyword arguments forwarded to
            ``self.multilayer.random_walk_rank``.

        Returns
        -------
        pandas.DataFrame
            Random-walk ranking results returned by HuMMuS.
        """
        if seeds is not None:
            self.seeds = seeds

        self.multilayer = self.Multixrank(
            self_loops=self_loops,
            restart_proba=restart_proba,
            verbose=verbose
        )
        self.results = self.multilayer.random_walk_rank(
            **random_walk_rank_kwargs
        )
        return self.results

    def rename_celltype(
        self,
        celltype_name
    ):

        self.multiplexes = {k.replace(self.celltype_name, celltype_name, 1): v
                            for k, v in self.multiplexes.items()
                            }

        self.bipartites = {
            celltype_name+"_grn-" + celltype_name + "_receptor": self.bipartites[
                self.celltype_name+"_grn-" + self.celltype_name + "_receptor"]\
        }
        for bipartite in self.bipartites:
            self.bipartites[bipartite]["source"] = \
                self.bipartites[bipartite]["source"].replace(
                self.celltype_name, celltype_name)
            self.bipartites[bipartite]["target"] = \
                self.bipartites[bipartite]["target"].replace(
                self.celltype_name, celltype_name)

        lamb = self.lamb
        lamb.index = lamb.index.str.replace(
            self.celltype_name+"_", celltype_name+"_")
        lamb.columns = lamb.columns.str.replace(
            self.celltype_name+"_", celltype_name+"_")
        self.lamb = lamb

        eta = self.eta
        eta.index = eta.index.str.replace(
            self.celltype_name+"_", celltype_name+"_")
        self.eta = eta

        self.celltype_name = celltype_name

            
class Multicell(Celltype):
    """Multicellular multilayer network.

    ``Multicell`` combines several :class:`Celltype` objects with a
    cell-cell communication layer. It supports the same :meth:`Multixrank` and
    :meth:`explore` workflow as ``Celltype``; ``explore`` creates
    ``self.multilayer``, runs ``random_walk_rank()``, stores ``self.results``,
    and returns the result table.
    """

    def __init__(
        self,
        celltypes: Union[
            List[Celltype],
            List[Dict],
            List[Dict[str, Celltype]]
        ],
        cell_communication_graph: pd.DataFrame,
        lamb=None,
        eta=None,
        cell_communication_graph_directed: bool = False,
        cell_communication_graph_weighted: bool = True,
        # bipartite parameters can be -1, 0, 1 here
        bipartite_grn_cell_communication_directed: Union[bool, int] = False,
        bipartite_grn_cell_communication_weighted: bool = False,
        bipartite_cell_communication_receptor_directed: bool = False,
        bipartite_cell_communication_receptor_weighted: bool = False,
        copy_graphs=True,
        seeds=[],
        verbose=True
    ):

        if not isinstance(bipartite_grn_cell_communication_directed,
                          (bool, int)):
            raise ValueError(
                "bipartite_grn_cell_communication_directed must"
                "be True, False, 0, 1, or -1.")
        if not isinstance(bipartite_cell_communication_receptor_directed,
                          (bool, int)):
            raise ValueError(
                "bipartite_cell_communication_receptor_directed must"
                " be True, False, 0, 1, or -1.")
        if not isinstance(celltypes, (list, dict)):
            raise ValueError("celltypes must be a list or a dictionary.")
        if not isinstance(bipartite_grn_cell_communication_weighted,
                          (bool, int)):
            raise ValueError(
                "bipartite_grn_cell_communication_weighted must"
                "be True, False, 0, or 1.")
        if not isinstance(bipartite_cell_communication_receptor_weighted,
                          (bool, int)):
            raise ValueError(
                "bipartite_cell_communication_receptor_weighted must"
                "be True, False, 0, or 1.")
        if not isinstance(cell_communication_graph, pd.DataFrame):
            raise ValueError(
                "cell_communication_graph must be a pandas DataFrame.")
        if not isinstance(cell_communication_graph_directed, (bool, int)):
            raise ValueError(
                "cell_communication_graph_directed must be a boolean.")
        if cell_communication_graph_weighted not in [True, False, 1, 0]:
            raise ValueError(
                "cell_communication_graph_weighted must be a boolean.")

        # Check if celltypes is a dictionary and convert it to a list of Celltype objects
        # If celltypes is a dictionary, the keys will be the new names of the celltypes
        self.celltypes_names = None
        if isinstance(celltypes, dict):
            self.celltypes_names = list(celltypes.keys())
            celltypes = list(celltypes.values())
            if verbose:
                warnings.warn(
                    "The celltypes dictionary was converted to" +
                    " a list of Celltype objects.\n" +
                    "The keys of the dictionary will be the celltype names.")
        else:
            self.celltypes_names = [
                celltype.celltype_name if isinstance(celltype, Celltype)
                else celltype["celltype_name"] for celltype in celltypes
                ]
        # Loop over celltypes and create Celltype objects
        for i in range(len(celltypes)):
            if isinstance(celltypes[i], dict):
                celltype_dict = celltypes[i]
                if "lamb" not in celltypes[i].keys():
                    celltype_dict["lamb"] = None
                celltypes[i] = Celltype(**celltype_dict)
            elif isinstance(celltypes[i], Celltype):
                if self.celltypes_names!=celltypes[i].celltype_name:
                    celltypes[i].rename_celltype(self.celltypes_names[i])
            elif not isinstance(celltypes[i], Celltype):
                raise ValueError("celltypes must be a list of Celltype objects or dictionaries.")

        # Rename celltypes if celltype_new_names is not None    
        if self.celltypes_names is not None:
            for i, celltype in enumerate(celltypes):
                celltype.celltype_name = self.celltypes_names[i]

        if copy_graphs: #test if when false not all copied already
            celltypes = list(celltypes)
            cell_communication_graph = cell_communication_graph.copy()

        # Store celltypes in a dictionary
        celltypes = {celltype.celltype_name: celltype for celltype in celltypes}

        # Save general parameters
        self.celltypes_names = list(celltypes.keys())
        self.seeds = seeds

        # Create bipartites between cell communication and celltype layer
        self.bipartites = {}
        for celltype in self.celltypes_names:
            mask_source_cell = cell_communication_graph["celltype_source"] == celltype
            mask_target_cell = cell_communication_graph["celltype_target"] == celltype

            # Create bipartites for each celltype
            # Link genes of each celltype to ligands of this celltype in cell communication
            bipartite_ligand = pd.DataFrame({
                "col1": cell_communication_graph[mask_source_cell]["source"].unique(),
                "col2": cell_communication_graph[mask_source_cell]["source"].unique() + '-' + celltype,
            })
            self.bipartites[celltype + "_to_ligands"] = {
                "source": "cell_communication",
                "target": celltype + "_grn",
                "graph_type": str(int(bool(bipartite_grn_cell_communication_directed)))\
                    + str(int(bipartite_grn_cell_communication_weighted)),
                "edge_list_df": bipartite_ligand
            }
            # Link receptors of each celltype to genes of this celltype in cell communication
            bipartite_receptor = pd.DataFrame({
                "col1": cell_communication_graph[mask_target_cell][
                    "target"].unique()+"_receptor",
                "col2": cell_communication_graph[mask_target_cell][
                    "target"].unique() + '-' + celltype
            })
            self.bipartites[celltype + "_to_receptors"] = {
                "source": "cell_communication",
                "target": celltype + "_receptor",
                "graph_type": str(int(bool(bipartite_cell_communication_receptor_directed)))\
                    + str(int(bipartite_cell_communication_receptor_weighted)),
                "edge_list_df": bipartite_receptor
            }

        # Create cell communication layer
        cell_communication_graph_type = str(int(cell_communication_graph_directed)) + str(int(cell_communication_graph_weighted))
        cell_communication_graph["source"] = cell_communication_graph["source"] + '-' + cell_communication_graph["celltype_source"]
        cell_communication_graph["target"] = cell_communication_graph["target"] + '-' + cell_communication_graph["celltype_target"]
        cell_communication_graph.rename(columns={"lr_means": "weight"}, inplace=True)
        cell_communication_graph["network_key"] = "cell_communication"

        self.multiplexes = {
            "cell_communication": {
                "names":["cell_communication"],
                "graph_type":[cell_communication_graph_type],
                "layers":[cell_communication_graph]
            }
        }

        # Update with layers and bipartites of each celltype
        for celltype in self.celltypes_names:
            self.multiplexes.update(celltypes[celltype].multiplexes)
            self.bipartites.update(celltypes[celltype].bipartites)
        
        # Add network_key column to all multiplex layers for sankey plots
        for key in self.multiplexes:
            if "network_key" not in self.multiplexes[key]["layers"][0].columns:
                self.multiplexes[key]["layers"][0]["network_key"] = key
        
        # Update nodes of each celtype to add a celltype specific suffix
        for celltype in self.celltypes_names:
            for multiplex in self.multiplexes:
                if celltype in multiplex:
                    self.multiplexes[multiplex]["layers"][0]["source"] = \
                        self.multiplexes[multiplex]["layers"][0]["source"] + \
                            "::" + celltype
                    self.multiplexes[multiplex]["layers"][0]["target"] = \
                        self.multiplexes[multiplex]["layers"][0]["target"] + \
                            "::" + celltype
            for bipartite in self.bipartites:
                if celltype in self.bipartites[bipartite]["source"]:
                    self.bipartites[bipartite]["edge_list_df"]["col2"] = \
                        self.bipartites[bipartite]["edge_list_df"]["col2"] + \
                            "::" + celltype
                if celltype in self.bipartites[bipartite]["target"]:
                    self.bipartites[bipartite]["edge_list_df"]["col1"] = \
                        self.bipartites[bipartite]["edge_list_df"]["col1"] + \
                            "::" + celltype

        # Prepare lamb matrix if not provided
        if lamb is None:
            lamb = pd.DataFrame(
                np.zeros((len(self.multiplexes), len(self.multiplexes))),
                index=list(self.multiplexes.keys()),
                columns=list(self.multiplexes.keys())
            )

            grn_layer_mask = lamb.index.str.endswith("_grn")
            receptor_layer_mask = lamb.index.str.endswith("_receptor")
            ccc_layer_mask = lamb.index.str.endswith("cell_communication")

            for celltype in self.celltypes_names:
                cell_type_mask = lamb.index.str.contains(celltype)

                lamb.loc[
                    grn_layer_mask*cell_type_mask,
                    receptor_layer_mask*cell_type_mask
                    + grn_layer_mask*cell_type_mask] = 1
                lamb.loc[
                    ccc_layer_mask,
                    ccc_layer_mask
                    + grn_layer_mask*cell_type_mask] += 1
                lamb.loc[
                    receptor_layer_mask*cell_type_mask,
                    ccc_layer_mask] = 1
            lamb = lamb.transpose().div(lamb.transpose().sum(1), 0)
        self.lamb = lamb

        # Fake eta
        if eta is None:
            self.eta = pd.Series(
                np.ones(len(self.multiplexes))/len(self.multiplexes),
                index=self.multiplexes.keys())
        else:
            self.eta = eta


    def illustrate_multicell(
        self,
        figsize=(12, 10),
        azim = 90,
        elev = 20,
        display_layer_axis=True,
        display_self_proba=True,
        display_layer_names=False,
        alpha_layers = 0.5,
        cell_communication_layer_name="cell_communication"
        ):

        lamb = self.lamb

        recon.plot.illustrate_multicell(
            lamb=lamb,
            figsize=figsize,
            azim=azim,
            elev=elev,
            display_layer_axis=display_layer_axis,
            display_self_proba=display_self_proba,
            display_layer_names=display_layer_names,
            alpha_layers=alpha_layers,
            cell_communication_layer_name=cell_communication_layer_name
        )

    def explore(
        self,
        seeds: Union[List, Dict, None] = None,
        self_loops=True,
        restart_proba=0.7,
        verbose=True,
        **random_walk_rank_kwargs
    ):
        """Build the multicellular multilayer object and run RWR.

        Parameters are the same as :meth:`Celltype.explore`. Multicellular seed
        names usually include the cell-type suffix, for example
        ``"IL6::Macrophage"`` for gene-layer seeds or ``"IL6-Macrophage"`` for
        cell-communication-layer seeds.
        """
        return super().explore(
            seeds=seeds,
            self_loops=self_loops,
            restart_proba=restart_proba,
            verbose=verbose,
            **random_walk_rank_kwargs
        )


def export_multicell_to_sif(multicell, filename="multicell_network.sif"):
    """Flatten a ``Multicell`` object into a single .sif file.

    Concatenates every multiplex layer and bipartite bridge of ``multicell``
    into one edge list, then cleans it into a static, non-redundant network:
    it drops the non-biological ``_receptor``/``_to_ligands``/``_to_receptors``
    layers (kept only to satisfy MultiXrank's requirement),
    collapses node aliases into one ``gene::celltype`` identity per layer,
    sets a placeholder ``+1`` sign on every interaction (undirected for now, need to 
    be updated on the future to expand for directed interactions), and writes the result to ``filename``.

    Parameters
    ----------
    multicell : Multicell
        A built ``Multicell`` object.
    filename : str, default="multicell_network.sif"
        Path the resulting tab-separated SIF file is written to.

    Returns
    -------
    pandas.DataFrame
        The exported edge list with ``source``, ``interaction``, and
        ``target`` columns, matching the file written to ``filename``.
    """
    edges = []  # accumulator: one small DataFrame per layer/bridge, concatenated at the end

    # Export multiplex layers: for N cell types there are 2N+1 layers
    # (N x grn, N x receptor, 1 shared cell_communication)
    for key, value in multicell.multiplexes.items():
        # value = {"names":..., "graph_type":..., "layers":[df]}; "layers" is a list but only
        # ever has 1 entry here, so [0] pulls out the actual edge table for this layer
        df = value["layers"][0]
        # keep only source/target — already in correct biological order for multiplex layers
        tmp = df[["source", "target"]].copy()
        # tag every row with the layer name (e.g. "B.cells_grn") so we know which layer it came from
        tmp["interaction"] = key
        edges.append(tmp)

    # Export bipartite layers: for N cell types there are 3N bridges connecting pairs of multiplexes
    # grn-receptor (NicheNet PKN): col2=receptor (source), col1=gene (target) — swap needed
    # to_ligands / to_receptors: col1=source, col2=target — no swap
    for key, value in multicell.bipartites.items():
        # value = {"source":..., "target":..., "graph_type":..., "edge_list_df": df}
        # edge_list_df uses generic col1/col2 (not source/target) since bipartites don't
        # share one fixed naming convention
        df = value["edge_list_df"]
        if "grn" in key and "receptor" in key:
            # only true for "{celltype}_grn-{celltype}_receptor" keys: here col1/col2 are
            # positionally reversed (col2 holds the receptor, col1 holds the gene), so swap
            # to get [source, target] = [receptor, gene]
            tmp = df[["col2", "col1"]].copy()  # receptor -> gene
        else:
            # matches "{celltype}_to_ligands" and "{celltype}_to_receptors": these already
            # have col1=source, col2=target, so no swap needed
            tmp = df[["col1", "col2"]].copy()  # gene -> CC-node
        # relabel whichever two columns were picked to the standard source/target names
        tmp.columns = ["source", "target"]
        # tag every row with the bridge name (e.g. "B.cells_to_ligands")
        tmp["interaction"] = key
        edges.append(tmp)

    # stack all multiplex tables + bipartite tables into one flat edge list
    edges = pd.concat(edges, ignore_index=True)
    # reorder to standard SIF column order: source, interaction, target
    edges = edges[["source", "interaction", "target"]]

    ### Network processing

    # --- Step 1: remove layers with no real biology ---
    # to_ligands / to_receptors are pure identity bridges 
    # grn-receptor is excluded from removal — it is the one bipartite with real
    # receptor->gene biology (e.g. NicheNet PKN) and must survive.
    remove_mask = (
        edges["interaction"].str.contains("to_ligands|to_receptors", regex=True) |
        (edges["interaction"].str.endswith("_receptor") & ~edges["interaction"].str.contains("grn-"))
    )
    edges = edges[~remove_mask].copy()

    # --- Step 2: unify all remaining node aliases into one gene::celltype identity format---
    # CC nodes:            APP-Macrophages.M1 -> APP::Macrophages.M1
    # grn-receptor source: A1BG_receptor::B.cells -> A1BG::B.cells
    # grn source:          E2F1_TF::B.cells  -> E2F1::B.cells
    def cc_node_to_grn(node):
        gene, celltype = node.rsplit("-", 1)
        return f"{gene}::{celltype}"

    cc_mask = edges["interaction"] == "cell_communication"
    edges.loc[cc_mask, "source"] = edges.loc[cc_mask, "source"].apply(cc_node_to_grn)
    edges.loc[cc_mask, "target"] = edges.loc[cc_mask, "target"].apply(cc_node_to_grn)

    grn_rec_mask = edges["interaction"].str.contains("grn-", regex=False)
    edges.loc[grn_rec_mask, "source"] = edges.loc[grn_rec_mask, "source"].str.replace("_receptor::", "::", regex=False)

    edges[["source", "target"]] = edges[["source", "target"]].apply(
        lambda x: x.str.replace("_TF::", "::", regex=False)
    )

    # --- Step 3: set all interactions to +1 (sign not available yet, to be worked on in the future) ---
    edges["interaction"] = 1

    # --- Step 4: replace remaining dashes in node names (e.g. HLA-DRA::B.cells -> HLA_DRA::B.cells) to avoid parsing issues ---
    edges[["source", "target"]] = edges[["source", "target"]].apply(
        lambda x: x.str.replace("-", "_", regex=False)
    )

    edges = edges[["source", "interaction", "target"]]
    edges.to_csv(filename, sep="\t", header=True, index=False)
    return edges


def format_multicell_results(
    multicell_multixrank_results,
    celltypes: List[str],
    keep_layers: Union[str, List[str]]="gene",
    split: str="::"
):
    """Format multicellular random-walk results as cell-type profiles.

    Parameters
    ----------
    multicell_multixrank_results : pandas.DataFrame
        Result table returned by ``Multicell.explore()`` or by
        ``Multicell.Multixrank(...).random_walk_rank()``. It must contain at
        least ``node``, ``layer``, and ``score`` columns.
    celltypes : list of str
        Cell types included in the multicellular analysis.
    keep_layers : str or list of str, default="gene"
        Layer name pattern(s) to retain. The default keeps gene-layer scores.
    split : str, default="::"
        Delimiter separating node names from cell-type names.

    Returns
    -------
    pandas.DataFrame
        Matrix with genes as rows, cell types as columns, and random-walk
        scores as values.
    """

    cell_type_profiles = multicell_multixrank_results[
        multicell_multixrank_results['layer'].str.contains(keep_layers)
    ].copy()
    cell_type_profiles[['gene', 'celltype']] = cell_type_profiles['node'].str.split(split, expand=True)
    cell_type_profiles = cell_type_profiles[cell_type_profiles['celltype'].isin(celltypes)]
    cell_type_profiles = cell_type_profiles[['celltype', 'gene', 'score']]

    # pivot to have genes as index and cell types as columns
    cell_type_profiles = cell_type_profiles.pivot(index='gene', columns='celltype', values='score')

    return cell_type_profiles


def set_lambda(
    multicell=None,
    celltypes=None,
    cell_communication_layer_name="cell_communication",
    direction: Union["upstream", "downstream"]="downstream",
    strategy: Union["intracell", "intercell"]="intercell",
    celltype_to_ccc_proba = None,
    ccc_to_celltype_proba = None
):

    if direction not in {"upstream", "downstream"}:
        raise ValueError("direction must be either 'upstream' or 'downstream'.")
    if strategy not in {"intracell", "intercell"}:
        raise ValueError("strategy must be either 'intracell' or 'intercell'.")

    if multicell is None and celltypes is None:
        raise ValueError("Either multicell or celltype should not be None.")
    elif multicell is not None:
        if celltypes is not None:
            warnings.warn(
                "Both multicell and celltypes are provided."
                "multicell will be used."
            )
        multiplexes = list(multicell.multiplexes.keys())
        celltypes = multicell.celltypes_names
    else:
        for celltype in celltypes:
            multiplexes.append(f"{celltype}_receptor")
            multiplexes.append(f"{celltype}_grn")
        multiplexes.append(cell_communication_layer_name)

    lamb = pd.DataFrame(
        np.zeros((len(multiplexes), len(multiplexes))),
        index=multiplexes,
        columns=multiplexes
    )
    
    # define masks
    grn_layer_mask = lamb.index.str.endswith("_grn")
    receptor_layer_mask = lamb.index.str.endswith("_receptor")
    ccc_layer_mask = lamb.index.str.endswith(cell_communication_layer_name)

    for celltype in celltypes:
        cell_type_mask = lamb.index.str.contains(celltype)

        # Transition FROM receptor and grn TO grn
        lamb.loc[
            receptor_layer_mask*cell_type_mask
            + grn_layer_mask*cell_type_mask,
            grn_layer_mask*cell_type_mask
            ] = 1

        # Transition FROM receptor TO receptor
        lamb.loc[
            receptor_layer_mask*cell_type_mask,
            receptor_layer_mask*cell_type_mask
            ] = 1
  

    # Transition FROM cell communication TO cell communication and receptor
    lamb.loc[
        ccc_layer_mask,
        ccc_layer_mask + receptor_layer_mask,
        ] = 1

    if strategy=="intercell":
        # Transition FROM grn TO cell communication
        lamb.loc[
            grn_layer_mask,
            ccc_layer_mask
            ] = 1

    # transpose if direction is upstream
    if direction=="upstream":
        lamb = lamb.transpose()

    lamb = lamb.div(lamb.sum(1), 0)

    return lamb



# To find targets : 

# Run intracellular predictions
## --> Direct targets

# Run intercellular predictions from prodcuted ligands
# (Corresponding to intracellular genes in the CCC)
## --> Indirect targets

def multicell_targets(
    seeds,
    celltypes,
    ccc,
    grn,
    receptor_grn: Union[str, pd.DataFrame],
    receptor_layer=None,
    grn_graph_directed=False,
    grn_graph_weighted=True,
    receptor_grn_graph_directed=False,
    receptor_grn_graph_weighted=True,
    receptor_graph_directed=False,
    receptor_graph_weighted=False,
    cell_communication_graph_directed=False,
    cell_communication_graph_weighted=True,
    bipartite_grn_cell_communication_directed=False,
    bipartite_grn_cell_communication_weighted=False,
    bipartite_cell_communication_receptor_directed=False,
    bipartite_cell_communication_receptor_weighted=False,
    restart_proba=0.6,
    extend_seeds=False,
    lamb = None,
    celltype_to_ccc_proba = None,
    ccc_to_celltype_proba = None,
    ccc_proba = 0.5,
    grn_proba = 0.5,
    njobs=-1,
    verbose=True,
    keep_layers="gene",
):
    """Predict direct and indirect perturbation effects across cell types.

    The function builds a multicellular ReCoN network, runs intracellular
    random walk with restart to estimate direct effects, then propagates the
    resulting ligand signal through cell-cell communication to estimate
    indirect effects.

    Parameters
    ----------
    seeds : list or dict
        Perturbed genes or weighted seed genes.
    celltypes : list or dict
        Cell-type names, ``Celltype`` objects, dictionaries accepted by
        ``Celltype``, or a mapping from names to those objects.
    ccc : pandas.DataFrame
        Cell-cell communication table with ligand ``source``, receptor
        ``target``, ``celltype_source``, ``celltype_target``, and weight
        columns.
    grn : pandas.DataFrame or str
        Gene regulatory network shared by the cell types, unless individual
        ``Celltype`` objects are provided.
    receptor_grn : pandas.DataFrame or str
        Receptor-to-gene bipartite network, or a loadable resource/path.
    restart_proba : float, default=0.6
        Random-walk restart probability.
    extend_seeds : bool, default=False
        If True, expand each seed to all cell types in the communication layer.
    njobs : int, default=-1
        Number of parallel jobs for indirect-effect runs.
    verbose : bool or int, default=True
        Controls progress messages and HuMMuS verbosity.
    keep_layers : str or list of str, default="gene"
        Layers retained for indirect-effect computations.

    Returns
    -------
    tuple[pandas.DataFrame, pandas.DataFrame]
        ``(direct_effect, indirect_effect)``. Direct effects are indexed by
        gene and have target cell types as columns. Indirect effects are indexed
        by gene and have a ``(celltype_target, celltype_source)`` column index.
    """
    if extend_seeds:
        if type(seeds) is not list and type(seeds) is not dict:
            seeds = list(seeds)
        starting_nodes = [f"{seed}-{celltype}" for seed in seeds for celltype in celltypes]
    else:
        starting_nodes = seeds

    # loading receptor-gene links
    if isinstance(receptor_grn, str):
        from recon.data.load_data import load_receptor_genes, receptor_gene_resources
        if receptor_grn in receptor_gene_resources:
            receptor_grn = load_receptor_genes(receptor_grn)
        else:
            try:
                receptor_grn = pd.read_csv(receptor_grn, sep=None, engine='python')
            except Exception as e:
                raise ValueError("receptor_grn should be a valid resource name or a path to a csv file.") from e
    if celltypes is None or len(celltypes) == 0:
        raise ValueError("celltypes should be a non-empty list of celltype names.")
    if type(celltypes) is list:
        for i in range(len(celltypes)):
            print(f"Processing celltype {i+1}/{len(celltypes)}: {celltypes[i]}")
            if not isinstance(celltypes[i], Celltype):
                celltypes[i] = Celltype(
                    celltype_name=celltypes[i],
                    grn_graph=grn,
                    receptor_grn_bipartite=receptor_grn,
                    receptor_graph=receptor_layer,
                    grn_graph_directed=grn_graph_directed,
                    grn_graph_weighted=grn_graph_weighted,
                    receptor_grn_bipartite_graph_directed=receptor_grn_graph_directed,
                    receptor_grn_bipartite_graph_weighted=receptor_grn_graph_weighted,
                    receptor_graph_directed=receptor_graph_directed,
                    receptor_graph_weighted=receptor_graph_weighted,
                )
        celltypes = {celltype.celltype_name: celltype for celltype in celltypes}

    elif type(celltypes) is dict:
        counter=0
        for celltype in celltypes.keys():
            print(f"Processing celltype {counter+1}/{len(celltypes)}: {celltype}")
            if not isinstance(celltypes[celltype], Celltype):
                celltypes[celltype] = Celltype(
                    **celltypes[celltype]
                )

    # Create a generic multicell object to precompute the lambdas
    generic_multicell = Multicell(
        celltypes=celltypes,
        cell_communication_graph=ccc,
        cell_communication_graph_directed=cell_communication_graph_directed,
        cell_communication_graph_weighted=cell_communication_graph_weighted,
        bipartite_grn_cell_communication_directed=bipartite_grn_cell_communication_directed,
        bipartite_grn_cell_communication_weighted=bipartite_grn_cell_communication_weighted,
        bipartite_cell_communication_receptor_directed=bipartite_cell_communication_receptor_directed,
        bipartite_cell_communication_receptor_weighted=bipartite_cell_communication_receptor_weighted,
        seeds=starting_nodes,
        verbose=verbose
    )
    celltype_names = list(celltypes.keys())

    print("Computing intracellular contributions and direct effect...")
    # Intracellular direct targets
    generic_multicell.lamb = set_lambda(
        generic_multicell,
        direction="downstream",
        strategy="intracell",
    )
    if ccc_to_celltype_proba is not None:
        ccc_to_celltype_proba = ccc_to_celltype_proba[celltype_names]
        for celltype in celltype_names:
            generic_multicell.lamb.loc["cell_communication", f"{celltype}_receptor"] = \
                (1-ccc_proba)*ccc_to_celltype_proba[celltype]/ccc_to_celltype_proba.sum()
    else:
        for celltype in celltype_names:
            generic_multicell.lamb.loc["cell_communication", f"{celltype}_receptor"] = \
                (1-ccc_proba)/len(celltype_names)

    generic_multicell.lamb.loc["cell_communication", "cell_communication"] = ccc_proba
    for celltype in celltype_names:
        generic_multicell.lamb.loc[f"{celltype}_grn", f"{celltype}_grn"] = grn_proba
        generic_multicell.lamb.loc[f"{celltype}_grn", "cell_communication"] = 1 - grn_proba
    
    multilayer = generic_multicell.Multixrank(restart_proba=restart_proba, verbose=verbose if verbose>=2 else False)
    intracell = multilayer.random_walk_rank()
    intracell = intracell.sort_values(ascending=False, by="node")
    intracell = intracell[intracell["layer"] == "gene"].set_index("node") 
    
    # Extracellular indirect regulations
    extracell = {}
    generic_multicell.lamb = set_lambda(
        generic_multicell,
        direction="downstream",
        strategy="intercell",
    )
    if ccc_to_celltype_proba is not None:
        ccc_to_celltype_proba = ccc_to_celltype_proba[celltype_names]
        for celltype in celltype_names:
            generic_multicell.lamb.loc["cell_communication", f"{celltype}_receptor"] = \
                (1-ccc_proba)*ccc_to_celltype_proba[celltype]/ccc_to_celltype_proba.sum()
    else:
        for celltype in celltype_names:
            generic_multicell.lamb.loc["cell_communication", f"{celltype}_receptor"] = \
                (1-ccc_proba)/len(celltype_names)
    generic_multicell.lamb.loc["cell_communication", "cell_communication"] = ccc_proba
    
    for celltype in celltype_names:
        generic_multicell.lamb.loc[f"{celltype}_grn", f"{celltype}_grn"] = grn_proba
        generic_multicell.lamb.loc[f"{celltype}_grn", "cell_communication"] = 1 - grn_proba

    print("Computing intercellular contributions and indirect effect...")
    # Precompute once (small/serializable)
    _targets = generic_multicell.multiplexes["cell_communication"]["layers"][0]["target"].unique()

    def _compute_for_celltype(celltype, keep_layers=["gene"]):
        # Avoid shared mutable state: fresh instance per task
        gm = copy.deepcopy(generic_multicell)   # or create_multicell() if deepcopy is heavy

        cell_seeds = intracell[intracell["multiplex"] == f"{celltype}_grn"]["score"].copy()
        cell_seeds.index = cell_seeds.index.str.replace("::", "-", regex=False)
        cell_seeds = cell_seeds[cell_seeds.index.isin(_targets)]

        gm.seeds = cell_seeds.to_dict()
        multilayer = gm.Multixrank(restart_proba=restart_proba, verbose=verbose if verbose>=2 else False)

        if keep_layers is not None:
            df = (
                multilayer.random_walk_rank()
                .sort_values("node", ascending=False)
                .query("layer == 'gene'")
                .set_index("node")[["score", "multiplex"]]
            )
        df["score"] = df["score"] * cell_seeds.sum()
        return celltype, df  # tuple is easily pickled

    # n_jobs: tune to your machine; -1 = all cores
    njobs = min(njobs, len(celltype_names)) if njobs > 0 else -1
    pairs = Parallel(n_jobs=njobs, verbose=10)(
        delayed(_compute_for_celltype)(ct, keep_layers) for ct in tqdm.tqdm(celltype_names)
    )
    extracell = dict(pairs)

    # Store cell type contributions
    cell_contributions = {}

    for celltype_target in celltype_names:
        cell_contributions[celltype_target] = intracell[
            intracell["multiplex"] == f"{celltype_target}_grn"
        ]["score"].to_frame()# * celltype_to_ccc_proba[celltype_target]
        cell_contributions[celltype_target].columns = [celltype_target+"_direct"]

        for celltype_source in celltype_names:
            if celltype_to_ccc_proba is not None:
                cell_contributions[celltype_target][celltype_source] = extracell[celltype_source][
                    extracell[celltype_source]["multiplex"] == f"{celltype_target}_grn"
                ]["score"] * celltype_to_ccc_proba[celltype_target]
            else:
                cell_contributions[celltype_target][celltype_source] = extracell[celltype_source][
                    extracell[celltype_source]["multiplex"] == f"{celltype_target}_grn"
                ]["score"]

    for name, df in cell_contributions.items():
        df.index = df.index.str.split("::").str[0]  # strip "::" + key name

    cell_contributions = pd.concat(cell_contributions, axis=1)

    return summarize_indirect_effects(cell_contributions)


def summarize_indirect_effects(predictions):
    """Split combined contribution tables into direct and indirect effects."""

    celltypes = predictions.columns.levels[0].tolist()

    indirect_effect = {}
    direct_effect = {}
    for celltype in celltypes:
        indirect_effect[celltype] = predictions.loc[:, pd.IndexSlice[celltype, celltypes]]
        direct_effect[celltype] = predictions.loc[:, pd.IndexSlice[celltype, celltype+'_direct']]

    direct_effect = pd.concat(direct_effect, axis=1)
    indirect_effect = pd.concat(indirect_effect, axis=1).droplevel(0, axis=1)
    direct_effect.index.name = 'gene'
    direct_effect.columns.name = 'celltype_target'
    indirect_effect.index.name = 'gene'
    indirect_effect.columns._names = ('celltype_target', 'celltype_source')

    return direct_effect, indirect_effect


def combine_effects(direct_effect, indirect_effect, alpha=0.8, cell_comm_matrix=None):
    """Combine direct and indirect perturbation effects.

    Parameters
    ----------
    direct_effect : pandas.DataFrame
        Direct-effect scores with genes as rows and target cell types as
        columns.
    indirect_effect : pandas.DataFrame
        Indirect-effect scores with genes as rows and a two-level column index
        ``(celltype_target, celltype_source)``.
    alpha : float, default=0.8
        Weight assigned to indirect effects. ``1 - alpha`` is assigned to
        direct effects.
    cell_comm_matrix : pandas.DataFrame, optional
        Optional cell-type communication weighting matrix. Rows should be
        target cell types and columns source cell types.

    Returns
    -------
    pandas.DataFrame
        Combined effect matrix with genes as rows and cell types as columns.
    """

    indirect_effect_summed = {}

    if cell_comm_matrix is None:
        for celltype in direct_effect:
            direct_effect[celltype] = direct_effect[celltype]/direct_effect[celltype].sum()
            indirect_effect_summed[celltype] = indirect_effect.loc[:, pd.IndexSlice[celltype, :]].sum(1) \
                / indirect_effect.loc[:, pd.IndexSlice[celltype, :]].sum(1).sum() \
                * direct_effect.loc[:, celltype].sum()

    else:
        cell_comm_matrix = cell_comm_matrix.div(cell_comm_matrix.sum(1), axis=0).fillna(0)
        for celltype in direct_effect:
            if celltype not in cell_comm_matrix.index:
                raise ValueError(f"Cell type {celltype} not found in cell communication matrix index.")
            direct_effect[celltype] = direct_effect[celltype]/direct_effect[celltype].sum()
            weighted_indirect = indirect_effect.loc[:, pd.IndexSlice[celltype, :]].droplevel(0, axis=1)\
                .dot(cell_comm_matrix.loc[celltype])
            indirect_effect_summed[celltype] = weighted_indirect\
                / weighted_indirect.sum()\
                * direct_effect[celltype].sum()

    return (1 - alpha) * direct_effect + alpha * pd.DataFrame(indirect_effect_summed)
