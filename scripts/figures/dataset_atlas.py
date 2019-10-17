import networkx as nx
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.pyplot import cm

from typing import Tuple, List, Optional
from itertools import islice

import sys
import os
sys.path.append(os.path.dirname(os.path.realpath(__file__)) + "/..")
from graph_io import read_metis_graph

try:
    import pydot
    from networkx.drawing.nx_pydot import graphviz_layout
except ImportError:
    raise ImportError("The dataset atlas needs pydot for the layout.")

import warnings
warnings.filterwarnings("ignore")


def read(*args, **kwargs) -> Tuple[nx.Graph, np.ndarray]:
    return read_metis_graph(*args, **kwargs)


def plot_atlas(paths: List[Path], output: Path, *, max_size: int = None, number_of_graphs: Optional[int] = None):
    graphs = (G for (G, S) in map(read, paths))

    if max_size is not None:
        graphs = (G for G in graphs if G.number_of_nodes() < max_size)

    UU = nx.disjoint_union_all(list(islice(graphs, number_of_graphs)))

    fig, ax = plt.subplots(figsize=(10, 10))

    ax.axis("off")
    ax.set_aspect("equal")

    # ‘dot’, ‘twopi’, ‘fdp’, ‘sfdp’, ‘circo’
    pos = graphviz_layout(UU, prog="neato")

    edge_color = np.array([c for (u, v, c) in UU.edges.data('cost')])
    edge_color = np.log1p(edge_color)

    nx.draw_networkx_nodes(UU, pos, ax=ax, node_size=30)
    nx.draw_networkx_edges(UU, pos, ax=ax, edge_color=edge_color, edge_cmap=cm.Blues)

    plt.savefig(output)


def main():
    ROOT = Path("../..")
    OUTPUT_DIR = ROOT / "scripts" / "figures"

    BIO_GRAPHS = list((ROOT / "data" / "bio").glob("*size-??.graph"))
    plot_atlas(BIO_GRAPHS, OUTPUT_DIR / "bio-atlas-15-50.pdf", max_size=15, number_of_graphs=50)
    plot_atlas(BIO_GRAPHS, OUTPUT_DIR / "bio-atlas-15-30.pdf", max_size=15, number_of_graphs=30)
    plot_atlas(BIO_GRAPHS, OUTPUT_DIR / "bio-atlas-15-10.pdf", max_size=15, number_of_graphs=10)

    BA_GRAPHS = list((ROOT / "data" / "barabasi-albert").glob("*.graph"))
    plot_atlas(BA_GRAPHS, OUTPUT_DIR / "barabasi-albert-atlas.pdf")

    DD_GRAPHS = list((ROOT / "data" / "duplication-divergence").glob("*.graph"))
    plot_atlas(DD_GRAPHS, OUTPUT_DIR / "duplication-divergence-atlas-15-*.pdf", max_size=15)
    plot_atlas(DD_GRAPHS, OUTPUT_DIR / "duplication-divergence-atlas-30-50.pdf", max_size=30, number_of_graphs=50)
    plot_atlas(DD_GRAPHS, OUTPUT_DIR / "duplication-divergence-atlas-30-20.pdf", max_size=30, number_of_graphs=20)

    PC_GRAPHS = list((ROOT / "data" / "powerlaw-cluster").glob("*.graph"))
    plot_atlas(PC_GRAPHS, OUTPUT_DIR / "powerlaw-cluster-atlas.pdf")


if __name__ == "__main__":
    main()
