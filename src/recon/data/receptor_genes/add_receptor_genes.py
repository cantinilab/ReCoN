import pandas as pd

# For mouse and human, we used the receptor to gene weights from the Nichenet prior knowledge network
# MOUSE
mouse_receptor_grn = pd.read_csv("old_multicell/nichenet/nnls_bipartite_receptor_gene_mouse.tsv", sep = "\t", header=None)
mouse_receptor_grn.columns = ["receptor", "grn", "score"]
mouse_receptor_grn = mouse_receptor_grn[mouse_receptor_grn["score"] > 0.005]

mouse_receptor_grn = mouse_receptor_grn.rename(columns={"grn": "target", "receptor": "source", "score": "weight"})[["source", "target", "weight"]]
mouse_receptor_grn.source = mouse_receptor_grn.source.str.partition("_receptor")[0].str.capitalize()
mouse_receptor_grn.to_parquet("ReCoN/data/receptor_genes/mouse_receptor_gene_from_NichenetPKN.parquet", index=False)

# HUMAN
human_receptor_grn = pd.read_csv("old_multicell/nichenet/nnls_bipartite_receptor_gene_human.tsv", sep = "\t", header=None)
human_receptor_grn.columns = ["receptor", "grn", "score"]
human_receptor_grn = human_receptor_grn[human_receptor_grn["score"] > 0.005]
human_receptor_grn = human_receptor_grn.rename(columns={"grn": "target", "receptor": "source", "score": "weight"})[["source", "target", "weight"]]
human_receptor_grn.source = human_receptor_grn.source.str.partition("_receptor")[0].str.upper()

human_receptor_grn.to_parquet("ReCoN/data/receptor_genes/human_receptor_gene_from_NichenetPKN.parquet", index=False)
