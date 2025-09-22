import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D, proj3d
from matplotlib.patches import FancyArrowPatch
from matplotlib import rcParams
from ..utils import split_layer_name

rcParams['font.family'] = 'sans-serif'
azim = np.random.randint(-180, 180)

class Arrow3D(FancyArrowPatch):
    def __init__(self, xs, ys, zs, *args, **kwargs):
        super().__init__((0, 0), (0, 0), *args, **kwargs)
        self._verts3d = xs, ys, zs

    def do_3d_projection(self, renderer=None):
        xs3d, ys3d, zs3d = self._verts3d
        xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, self.axes.M)
        self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))

        return np.min(zs)


def illustrate_multicell(
    lamb,
    figsize=(12, 10),
    azim: float = 90,
    elev: float = 20,
    display_layer_axis=True,
    display_self_proba=True,
    display_layer_names=False,
    alpha_layers: float = 0.5,
    cell_communication_layer_name: str = "cell_communication"
):

    """
    Visualize the transition probabilities between layers in a multicellular
    multilayer model.

    Parameters
    ----------
    lamb : pd.DataFrame
        A DataFrame representing the transition probabilities between layers.
    figsize : tuple, optional
        Size of the figure (width, height). Default is (12, 10).
    azim : int, optional
        Azimuth angle for the 3D plot view. Default is 90.
    elev : int, optional
        Elevation angle for the 3D plot view. Default is 20.
    display_layer_axis : bool, optional
        Whether to display the layer axis. Default is True.
    display_self_proba : bool, optional
        Whether to display self-transition probabilities. Default is True.
    display_layer_names : bool, optional
        Whether to display layer names on the plot. Default is False.
    alpha_layers : float, optional
        Transparency level for the layer planes. Default is 0.5.
    cell_communication_layer_name : str, optional
        Name of the cell communication layer. Default is "cell_communication".

    Returns
    -------
    None
    """

    # Adjust the viewing angle to be more horizontal
    fig = plt.figure(figsize=(12, 10), dpi=200)
    ax = fig.add_subplot(111, projection='3d')

    # Set a lower elevation angle to create a more horizontal view
    ax.view_init(elev=elev, azim=azim)  # Set elevation to 20 for a more side view

    layers = lamb.index.values
    lamb = lamb.loc[layers, layers]
    cell_groups = set([split_layer_name(layer)[0] for layer in layers
                       if layer != cell_communication_layer_name])
    cell_groups.add(cell_communication_layer_name)

    # Define the x-coordinate for each cell type's column
    celltypes = set([split_layer_name(layer)[0] for layer in layers
                     if layer != cell_communication_layer_name])
    x_positions = {celltype: float(i) for i, celltype in enumerate(celltypes)}
    x_max, x_min = max(x_positions.values()), min(x_positions.values())
    x_positions[cell_communication_layer_name] = (x_min + x_max) / 2

    layer_types = list(set([split_layer_name(layer)[-1] for layer in layers
                            if layer != cell_communication_layer_name]))
    layer_heights = {
        layer_type: i for i, layer_type in enumerate(layer_types)
    }
    layer_heights[cell_communication_layer_name] = len(layer_heights)
    layer_types.append(cell_communication_layer_name)

    z_layers = {
        layer: layer_heights[split_layer_name(layer)[-1]]
        for layer in layers if layer != cell_communication_layer_name}
    z_layers[cell_communication_layer_name] = layer_heights[
        cell_communication_layer_name]

    if len(cell_groups) < 10:
        celltype_colors = {
            celltype: color for celltype, color in zip(
                cell_groups, plt.cm.tab10.colors)}
    elif len(cell_groups) < 20:
        celltype_colors = {
            celltype: color for celltype, color in zip(
                cell_groups, plt.cm.tab20.colors)}
    else:
        celltype_colors = {
            celltype: plt.cm.tab20.colors[i] for i, celltype in enumerate(
                cell_groups)}

    # Draw the "cell_communication" layer as a large plane spanning across cols
    Z = z_layers[cell_communication_layer_name]

    u = np.linspace(x_min - 0.4, x_max + 0.4, 2)
    v = np.linspace(-0.3, 0.3, 2)
    U, V = np.meshgrid(u, v)
    W = Z * np.ones_like(U)
    ax.plot_surface(
        U, V, W,
        color=celltype_colors[cell_communication_layer_name],
        label=cell_communication_layer_name,
        alpha=alpha_layers)

    if display_layer_names:
        ax.text(
            (x_max+x_min)/2+0.3, 0.3, Z,  # Position of the layer name
            cell_communication_layer_name,  # Layer name
            "x", color="black", ha="center", va="top", fontsize=10)

    # Draw planes for each layer in their respective columns
    for layer in layers:
        if layer == cell_communication_layer_name:
            continue  # Skip since cell_communication is already drawn

        celltype = split_layer_name(layer)[0]
        x = x_positions.get(celltype, 0)
        z = z_layers[layer]
        color = celltype_colors.get(celltype, "gray")

        # Draw each layer's plane at specified x and z positions
        u = np.linspace(x - 0.4, x + 0.4, 2)
        v = np.linspace(-0.3, 0.3, 2)
        U, V = np.meshgrid(u, v)
        W = z * np.ones_like(U)
        ax.plot_surface(U, V, W, color=color, alpha=alpha_layers)
        if display_layer_names:
            ax.text(
                x, 0.3, z,  # Position of the layer name
                layer,  # Layer name
                "x", color="black", ha="center", va="top", fontsize=10)

    # Draw edges based on probabilities with positions and heights specified
    for i, layer_from in enumerate(lamb.index):
        for j, layer_to in enumerate(lamb.columns):
            prob = lamb.iloc[i, j]
            if i == j:
                if layer_from != cell_communication_layer_name:
                    x1, z1 = x_positions.get(
                        split_layer_name(
                            layer_from)[0], 0), z_layers[layer_from]
                else:
                    x1, z1 = x_positions.get(
                        layer_from, 0), z_layers[layer_from]
                if display_self_proba:
                    ax.text(
                        x1, 0, z1,  # Position of the probability of self-transition
                        round(prob, 2),  # Displayed probability of self-transition
                        "x", color="black", ha="center", va="center",
                        fontsize=10)
            elif prob != 0:
                x1, z1 = x_positions.get(
                    split_layer_name(layer_from)[0], 0), z_layers[layer_from]
                x2, z2 = x_positions.get(
                    split_layer_name(layer_to)[0], 0), z_layers[layer_to]
                if layer_from == cell_communication_layer_name:
                    x1 = x2
                elif layer_to == cell_communication_layer_name:
                    x2 = x1

                if x1 > x2:
                    y1 = 0.2
                    y2 = 0.2
                elif x1 <= x2:
                    y1 = -0.2
                    y2 = -0.2

                if z1 > z2:
                    y1, y2 = -0.0, -0.0,
                    x1 -= 0.2
                    x2 -= 0.2
                elif z1 <= z2:
                    y1, y2 = 0.0, 0.0
                    x1 += 0.2
                    x2 += 0.2

                ax.add_artist(
                    Arrow3D(
                        [x1, x2], [y1, y2], [z1, z2],
                        mutation_scale=15, lw=prob * 5,
                        arrowstyle="-|>", color="black", alpha=0.8))

    # Set plot limits and labels
    ax.set_xlim(x_min-0.4, x_max+0.4)
    ax.set_ylim(-0.5, 0.3)
    ax.set_zlim(-0.0, len(layer_heights)-1)

    ax.set_xlabel("Cell Types")
    ax.set_ylabel("")
    ax.set_zlabel("Layer Height")

    ax.set_title("Layer Transition Probabilities in Multilayer Graph (Adjusted View)")

    # Generate x-ticks and labels for cell types
    x_dict = {
        celltype: x for celltype, x in x_positions.items()
        if celltype != cell_communication_layer_name}
    ax.set_xticks(list(x_dict.values()))
    ax.set_xticklabels(list(x_dict.keys()))

    # Generate z-ticks and labels for layer types
    z_dict = {layer_type: z for layer_type, z in layer_heights.items()}
    ax.set_zticks(list(z_dict.values()))
    ax.set_zticklabels(list(z_dict.keys()))

    # Hide y-axis
    ax.set_yticks([])  # Hide x-axis ticks
    ax.set_yticklabels([])  # Hide x-axis tick labels

    tmp_planes = ax.zaxis._PLANES
    ax.zaxis._PLANES = (tmp_planes[3], tmp_planes[2],
                        tmp_planes[5], tmp_planes[4],
                        tmp_planes[0], tmp_planes[1])

    if not display_layer_axis:
        ax.set_zticks([])  # Hide z-axis ticks
        ax.set_zticklabels([])  # Hide z-axis tick labels
        ax.set_zlabel("")
        ax.set_xticks([])  # Hide x-axis ticks
        ax.set_xticklabels([])  # Hide x-axis tick labels
        ax.set_xlabel("")

    # Hide the grid and axis lines
    ax.grid(False)
    ax.xaxis.pane.set_alpha(0.0)
    ax.yaxis.pane.set_alpha(0.0)
    ax.zaxis.pane.set_alpha(0.0)

    # Hide axis lines while keeping ticks
    ax.xaxis.line.set_color((1.0, 1.0, 1.0, 0.0))  # Hide x-axis line
    ax.yaxis.line.set_color((1.0, 1.0, 1.0, 0.0))  # Hide y-axis line
    ax.zaxis.line.set_color((1.0, 1.0, 1.0, 0.0))  # Hide z-axis line

    ax.set_box_aspect((np.sqrt(len(lamb)), 1, 1.5))
    plt.show()
