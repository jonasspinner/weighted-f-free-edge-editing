import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import string

import pydot
from networkx.drawing.nx_pydot import graphviz_layout


def plot_small_p3_editing_example():
    np.random.seed(1)

    G_edited = nx.disjoint_union_all([nx.complete_graph(n) for n in [3, 4]])
    G_edited = nx.relabel_nodes(G_edited, dict(zip(G_edited, string.ascii_lowercase)))
    edits = [("b", "g"), ("d", "e")]

    deletions = set(edits) - set(G_edited.edges())
    inserts = set(edits) - deletions

    G = G_edited.copy()
    G.add_edges_from(deletions)
    G.remove_edges_from(inserts)

    graphs = [
        ("P3-editing-example-G", G, set(G.edges()) - set(edits), inserts, deletions, [("b", "e")],
         [("b", "e"), ("b", "g"), ("e", "g"), ("d", "e")]),
        ("P3-editing-example-G-edited", G_edited, G_edited.edges(), (), (), (), ())
    ]

    for name, G, unedited, inserted, deleted, non_edges, fat in graphs:
        fig, ax = plt.subplots(figsize=(4, 4))
        ax.margins(0.15)
        ax.set_aspect("equal")
        ax.axis(False)

        def width(edges):
            return [3 if e in fat else 1 for e in edges]

        pos = graphviz_layout(G, prog="neato")
        nx.draw_networkx_nodes(G, pos, node_color="lightgrey", node_size=400, ax=ax)
        nx.draw_networkx_labels(G, pos, labels={v: f"${v}$" for v in G}, ax=ax)

        edge_lists = [
            (unedited,  "black", "solid"),
            (inserted,  "green", "solid"),
            (deleted,   "red",   "solid"),
            (non_edges, "grey", "dashed")
        ]

        for edgelist, edge_color, style in edge_lists:
            nx.draw_networkx_edges(G, pos, edgelist=edgelist, edge_color=edge_color, style=style, ax=ax, width=width(edgelist))

        fig.tight_layout()
        plt.savefig(f"{name}.pdf", bbox_inches="tight", pad_inches=0)


def plot_small_graphs():
    P5 = nx.path_graph(5)
    C5 = nx.cycle_graph(5)
    _5 = [(-0.5,-0.688 ), (-0.809,0.2628), (0,0.85065), (0.809,0.2628), (0.5,-0.688)]

    necktie = nx.from_edgelist([(0, 1), (1, 2), (2, 3), (2, 4), (3, 4)])
    bowtie = nx.from_edgelist([(0, 1), (0, 2), (1, 2), (2, 3), (2, 4), (3, 4)])

    graphs = [
        ("P3", nx.path_graph(3), [(0, 0), (1, 0), (0.5, 0.5**(1/3))]),
        ("P4", nx.path_graph(4), [(0, 0), (1, 0), (1, 1), (0, 1)]),
        ("C4", nx.cycle_graph(4), [(0, 0), (1, 0), (1, 1), (0, 1)]),
        ("2K2", nx.from_edgelist([(0, 1), (2, 3)]), [(0, 0), (0, 1), (1, 1), (1, 0)]),
        ("P5", P5, _5),
        ("C5", C5, _5),
        ("necktie", necktie, [(0, -1.5), (0, -0.75), (0, 0), (-0.5, 0.5**(1/3)), (0.5, 0.5**(1/3))]),
        ("bowtie", bowtie, [(-0.5**(1/3), 0.5), (-0.5**(1/3), -0.5), (0, 0), (0.5**(1/3), -0.5), (0.5**(1/3), 0.5)])
    ]

    for name, G, pos in graphs:
        fig, ax = plt.subplots(figsize=(2, 2))
        ax.set_aspect("equal")
        ax.axis(False)
        ax.margins(0.15)

        nx.draw_networkx(G, pos, with_labels=False, node_color="lightgrey", node_size=300, width=1, ax=ax)

        fig.tight_layout()
        plt.savefig(f"{name}.pdf", bbox_inches="tight", pad_inches=0)


def main():
    plot_small_p3_editing_example()
    plot_small_graphs()


if __name__ == '__main__':
    main()
