import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import mark_inset
import numpy as np

def plot_celltype_comparison(df, celltype1, celltype2, quantile=0.999):
    """
    Scatter plot of two cell types, highlighting outliers where
    celltype1 >> celltype2 (far to the right of diagonal), with a side zoom panel.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame with at least the two columns.
    celltype1 : str
        Column name for X-axis.
    celltype2 : str
        Column name for Y-axis.
    quantile : float
        Quantile threshold for defining outliers (default=0.999).
    """
    # Compute mask of top outliers (genes overexpressed in celltype1 relative to celltype2)
    diff = (df[celltype1] - df[celltype2]).sort_values(ascending=False)
    mask = diff > diff.quantile(quantile)

    fig, ax = plt.subplots(figsize=(6, 6))

    # Main scatter
    ax.scatter(
        df.loc[~mask, celltype1],
        df.loc[~mask, celltype2],
        color="gray", alpha=0.7, label="Normal", s=5
    )
    ax.scatter(
        df.loc[mask, celltype1],
        df.loc[mask, celltype2],
        color="blue", alpha=0.8, label=f"{celltype1} >> {celltype2}", s=20
    )

    # Add diagonal reference line
    lims = [
        min(df[celltype1].min(), df[celltype2].min()),
        max(df[celltype1].max(), df[celltype2].max()),
    ]
    ax.plot(lims, lims, "k--", alpha=0.5)

    ax.set_xlabel(f"{celltype1} effect")
    ax.set_ylabel(f"{celltype2} effect")
    ax.legend()
    ax.set_title(f"{celltype2} vs {celltype1} effect")

    # ---------------------------
    # Side zoom panel
    # ---------------------------
    axins = fig.add_axes([1.05, 0.15, 0.6, 0.6])  # panel on right side

    # Scatter again but with bigger points
    axins.scatter(
        df.loc[mask, celltype1],
        df.loc[mask, celltype2],
        color="blue", alpha=0.8, s=60
    )

    # Annotate outliers with index (gene names)
    for idx, row in df.loc[mask].iterrows():
        axins.annotate(
            str(idx),
            (row[celltype1], row[celltype2]),
            xytext=(0, 8),
            textcoords="offset points",
            fontsize=9, color="blue",
            ha="center", va="bottom"
        )

    # Set zoomed view limits (pad around outliers)
    xvals = df.loc[mask, celltype1]
    yvals = df.loc[mask, celltype2]
    axins.set_xlim(0, xvals.max() + 0.1 * np.abs(xvals.max()))
    axins.set_ylim(0, yvals.max() + 0.1 * np.abs(yvals.max()))
    axins.set_xlabel(f"{celltype1} effect", fontsize=10)
    axins.set_ylabel(f"{celltype2} effect", fontsize=10)
    axins.tick_params(axis="both", which="major", labelsize=8)
    axins.set_title(f"Top perturbed genes in {celltype1}", fontsize=12)

    # Try connector box (if inset is within bounds)
    try:
        mark_inset(ax, axins, loc1=2, loc2=4, fc="none", ec="0.5")
    except Exception:
        pass

    plt.show()
