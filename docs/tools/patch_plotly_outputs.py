#!/usr/bin/env python3
"""
Convert Plotly outputs stored as
application/vnd.plotly.v1+json
into text/html so Sphinx/MyST-NB can render them.

No notebook execution.
"""

import json
import nbformat
import plotly.io as pio
from pathlib import Path

SOURCE_DIR = Path("source")

def patch_notebook(nb_path: Path) -> bool:
    nb = nbformat.read(nb_path, as_version=4)
    changed = False

    for cell in nb.cells:
        if cell.get("cell_type") != "code":
            continue

        for out in cell.get("outputs", []):
            data = out.get("data", {})
            if "application/vnd.plotly.v1+json" not in data:
                continue

            fig_spec = data["application/vnd.plotly.v1+json"]

            # Convert Plotly JSON → HTML
            fig = pio.from_json(
                json.dumps(fig_spec)
                if isinstance(fig_spec, dict)
                else fig_spec
            )

            html = fig.to_html(
                full_html=False,
                include_plotlyjs="cdn"
            )

            data["text/html"] = html
            changed = True

    if changed:
        nbformat.write(nb, nb_path)

    return changed


def main():
    notebooks = list(SOURCE_DIR.rglob("*.ipynb"))
    patched = 0

    for nb in notebooks:
        if patch_notebook(nb):
            patched += 1
            print(f"[plotly] patched {nb}")

    if patched:
        print(f"[plotly] updated {patched} notebook(s)")
    else:
        print("[plotly] nothing to patch")


if __name__ == "__main__":
    main()
