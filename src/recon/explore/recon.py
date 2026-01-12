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
                    "a list of Celltype objects.\n" +
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


def format_multicell_results(
    multicell_multixrank_results,
    celltypes: List[str],
    keep_layers: Union[str, List[str]]="gene",
    split: str="::"
):
    """Format multicellular multixrank results as gene profiles per cell type.

    Args:
        multicell_multixrank_results: Results from multicellular multixrank random walk.
        celltypes: List of cell types included in the multicellular analysis.
        keep_layers: Layers to keep in the output profiles. Can be "gene", "grn", or a list of specific layers.
        split: Delimiter used to separate gene names and cell type names in the node identifiers.

    Returns:
        DataFrame with cell types as columns,  genes as rows, and random walk with restart scores as values.
    """

    cell_type_profiles = multicell_multixrank_results[multicell_multixrank_results['layer'].str.contains(keep_layers)]
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

    if multicell is None and celltypes is None:
        raise ValueError("Either multicell or celltype should not be None.")
    elif multicell is not None:
        if celltypes is not None:
            raise warnings.warn("Both multicell and celltypes are provided."+
                "multicell will be used.")
        else:
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
    if extend_seeds:
        if type(seeds) is not list and type(seeds) is not dict:
            seeds = list(seeds)
        starting_nodes = [f"{seed}-{celltype}" for seed in seeds for celltype in celltypes]
    else:
        starting_nodes = seeds

    # loading receptor-gene links
    if isinstance(receptor_grn, str):
        from recon.loader import load_receptor_target, receptor_target_resources
        if receptor_grn in receptor_target_resources:
            receptor_grn = load_receptor_target(receptor_grn)
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

    print("Computing intracellular contributions and direct effect...")
    # Intracellular direct targets
    generic_multicell.lamb = set_lambda(
        generic_multicell,
        direction="downstream",
        strategy="intracell",
    )
    if ccc_to_celltype_proba is not None:
        ccc_to_celltype_proba = ccc_to_celltype_proba[celltypes]
        for celltype in celltypes:
            generic_multicell.lamb.loc["cell_communication", f"{celltype}_receptor"] = \
                (1-ccc_proba)*ccc_to_celltype_proba[celltype]/ccc_to_celltype_proba.sum()
    else:
        for celltype in celltypes:
            generic_multicell.lamb.loc["cell_communication", f"{celltype}_receptor"] = \
                (1-ccc_proba)/len(celltypes)

    generic_multicell.lamb.loc["cell_communication", "cell_communication"] = ccc_proba
    for celltype in celltypes:
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
        ccc_to_celltype_proba = ccc_to_celltype_proba[celltypes]
        for celltype in celltypes:
            generic_multicell.lamb.loc["cell_communication", f"{celltype}_receptor"] = \
                (1-ccc_proba)*ccc_to_celltype_proba[celltype]/ccc_to_celltype_proba.sum()
    else:
        for celltype in celltypes:
            generic_multicell.lamb.loc["cell_communication", f"{celltype}_receptor"] = \
                (1-ccc_proba)/len(celltypes)
    generic_multicell.lamb.loc["cell_communication", "cell_communication"] = ccc_proba
    
    for celltype in celltypes:
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
    njobs = min(njobs, len(celltypes)) if njobs > 0 else -1
    pairs = Parallel(n_jobs=njobs, verbose=10)(
        delayed(_compute_for_celltype)(ct, keep_layers) for ct in tqdm.tqdm(celltypes)
    )
    extracell = dict(pairs)

    # Store cell type contributions
    cell_contributions = {}

    for celltype_target in celltypes:
        cell_contributions[celltype_target] = intracell[
            intracell["multiplex"] == f"{celltype_target}_grn"
        ]["score"].to_frame()# * celltype_to_ccc_proba[celltype_target]
        cell_contributions[celltype_target].columns = [celltype_target+"_direct"]

        for celltype_source in celltypes:
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
            indirect_effect[celltype] = indirect_effect[celltype].dot(cell_comm_matrix.loc[celltype])
            indirect_effect[celltype] = indirect_effect[celltype]\
                / indirect_effect[celltype].sum()\
                * direct_effect[celltype].sum()

    return (1 - alpha) * direct_effect + alpha * pd.DataFrame(indirect_effect_summed)