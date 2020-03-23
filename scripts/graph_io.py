import networkx as nx
import numpy as np
from pathlib import Path

from typing import Union, Tuple, Optional, List


def read_metis_graph(path: Union[str, Path]) -> Tuple[nx.Graph, np.ndarray]:
    """
    Read a metis file from `path`. The instance must be a fully connected undirected weighted graph.
    Returns a graph with edges where c_{uv} > 0 and a matrix with the similarity scores.

    References
    ----------
        [1]: https://people.sc.fsu.edu/~jburkardt/data/metis_graph/metis_graph.html
    """
    path = Path(path)
    with path.open("r") as f:
        lines = f.read().split("\n")

    n, m, fmt, *_ = list(map(int, lines[0].split(" "))) + [0]
    graph = nx.Graph()
    graph.name = path.name
    graph.add_nodes_from(range(n))

    S = np.zeros((n, n))

    if fmt == 0 or fmt is None:
        # unweighted, default

        S[:] = 1

        for u, line in enumerate(lines[1:]):
            for v in [int(v) - 1 for v in line.split()]:
                graph.add_edge(u, v)

    elif fmt == 1:
        # weighted
        assert m == n * (n - 1) // 2

        for u, line in enumerate(lines[1:]):
            it = iter(line.split())
            for v, c in [(int(v) - 1, float(c)) for (v, c) in zip(it, it)]:
                if c >= 0:
                    graph.add_edge(u, v, cost=c)
                S[u, v] = S[v, u] = c

    else:
        raise Exception(f"fmt \"{fmt}\" not supported")

    return graph, S


def write_metis_graph(path: Union[str, Path], graph: nx.Graph, costs: Optional[np.ndarray] = None, *,
                      comments: Optional[List[str]] = None) -> None:
    """
    Writes graph to file. Either writes a unweighted (fmt = 0) or a weighted instance (fmt = 1) when `costs` is given.

    For weighted instances the output is a fully connected graph. Only the upper triangular adjacency matrix is being
    written, i.e. all vertex pairs where u < v. The edge weight is negative if no edge is present in the graph.

    For unweighted instances the output is the graph directly.

    Parameters
    ----------
    path : Path
    graph : Graph
    costs : symmetric costs matrix, optional
    comments : list of comments, optional

    References
    ----------
        [1]: https://people.sc.fsu.edu/~jburkardt/data/metis_graph/metis_graph.html
    """
    with Path(path).open('w') as file:
        if comments is not None:
            for comment in comments:
                file.write(f"% {comment}\n")

        if costs is None:
            n = graph.number_of_nodes()
            m = graph.number_of_edges()
            fmt = 0

            file.write(f"{n} {m} {fmt}\n")

            for u in range(n):
                file.write(" ".join([f"{v + 1}" for v in graph.neighbors(u)]) + "\n")

        else:
            n = graph.number_of_nodes()
            m = n * (n - 1) // 2
            fmt = 1

            file.write(f"{n} {m} {fmt}\n")

            for u in range(n):
                file.write(" ".join([f"{v + 1} {costs[u][v] if graph.has_edge(u, v) else -costs[u][v]}"
                                     for v in range(u + 1, n)]) + "\n")
