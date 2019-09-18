import networkx as nx
import sys
from pathlib import Path
import matplotlib.pyplot as plt
import numpy as np
from typing import Tuple, List
import argparse

from graph_io import write_metis_graph


DEFAULTS = dict(
    barabasi=dict(
        m=4
    ),
    powerlaw_cluster=dict(
        m=4,
        p=0.4
    ),
    duplication_divergence=dict(
        p=0.2
    ),
    mu1=10,
    mu2=-10,
    sigma1=5,
    sigma2=5
)


def generate_similarity_matrix(graph: nx.Graph,
                               mu1: float = DEFAULTS['mu1'], mu2: float = DEFAULTS['mu2'],
                               sigma1: float = DEFAULTS['sigma1'], sigma2: float = DEFAULTS['sigma2']) -> np.ndarray:
    """
    Generates a similarity matrix for the given graph. Similarity scores are drawn are either drawn from N(mu1, sigma1)
    for edges or N(mu2, sigma2) for non-edges. N(mu, sigma) is the normal distribution with mean mu and standard
    deviation sigma.

    The similarity matrix is symmetric an has a zero diagonal.

    Parameters
    ----------
    graph : nx.Graph
    mu1, mu2 : number
        The mean value for edges and non-edges similarities respectively.
    sigma1, sigma2 : number
        The standard deviation for edges and non-edges similarities respectively.

    Returns
    -------
    S : np.ndarray
    """
    n = graph.number_of_nodes()
    S = np.zeros((n, n))

    for u in graph.nodes:
        for v in graph.nodes:
            if u < v:
                if graph.has_edge(u, v):
                    S[u, v] = S[v, u] = np.random.normal(mu1, sigma1)
                else:
                    S[u, v] = S[v, u] = np.random.normal(mu2, sigma2)
    return S


def build_instance(S: np.ndarray) -> Tuple[nx.Graph, np.ndarray]:
    """
    Builds an instance from a similarity matrix.

    Parameters
    ----------
    S : np.ndarray

    Returns
    -------
    G : nx.Graph
    costs : np.ndarray
    """
    G = nx.from_numpy_array(S > 0)
    costs = np.abs(S)
    return G, costs


def generate_barabasi_albert_instance(directory: Path, nr: int, n: int,
                                      m: int = DEFAULTS['barabasi']['m'],
                                      seed=None) -> None:
    """
    Generate a graph from the barabasi albert model and build an instance from it.
    """
    graph = nx.generators.barabasi_albert_graph(n, m, seed=seed)
    similarity = generate_similarity_matrix(graph)
    G, cost = build_instance(similarity)

    path = directory / f"barabasi-albert-nr-{nr}-size-{n}-m-{m}-seed-{seed}.graph"
    print(f"Generated {path.name}")
    write_metis_graph(path, G, cost)


def generate_powerlaw_cluster_instance(directory: Path, nr: int, n: int,
                                       m: int = DEFAULTS['powerlaw_cluster']['m'],
                                       p: float = DEFAULTS['powerlaw_cluster']['p'],
                                       seed=None) -> None:
    """
    Generate a graph from the powerlaw cluster model and build an instance from it.
    """
    graph = nx.generators.powerlaw_cluster_graph(n, m, p, seed=seed)
    similarity = generate_similarity_matrix(graph)
    G, cost = build_instance(similarity)

    path = directory / f"powerlaw-cluster-nr-{nr}-size-{n}-m-{m}-p-{p}-seed-{seed}.graph"
    print(f"Generated {path.name}")
    write_metis_graph(path, G, cost)


def generate_duplication_divergence_instance(directory: Path, nr: int, n: int,
                                             p: float = DEFAULTS['duplication_divergence']['p'],
                                             seed=None) -> None:
    """
    Generate a graph from the duplication divergence model and build an instance from it.
    """
    graph = nx.generators.duplication_divergence_graph(n, p, seed=seed)
    similarity = generate_similarity_matrix(graph)
    G, cost = build_instance(similarity)

    path = directory / f"duplication-divergence-nr-{nr}-size-{n}-p-{p}-seed-{seed}.graph"
    print(f"Generated {path.name}")
    write_metis_graph(path, G, cost)


def generate_dataset(path: Path, model: str, num_nodes: List[int], seed=None) -> None:
    """
    Generate instances from the given model.
    """
    path.mkdir(exist_ok=True)

    for i, n in enumerate(num_nodes):
        nr = i + 1
        if model == "barabasi-albert":
            generate_barabasi_albert_instance(path, nr, n, seed=seed)
        elif model == "powerlaw-cluster":
            generate_powerlaw_cluster_instance(path, nr, n, seed=seed)
        elif model == "duplication-divergence":
            generate_duplication_divergence_instance(path, nr, n, seed=seed)
        else:
            raise RuntimeError("invalid model")


def plot(seed=None):
    G0 = nx.generators.barabasi_albert_graph(100, 4, seed=seed)
    G0 = nx.generators.powerlaw_cluster_graph(100, 4, 0.4, seed=seed)
    G0 = nx.generators.duplication_divergence_graph(100, 0.2, seed=seed)
    S = generate_similarity_matrix(G0)
    G, costs = build_instance(S)

    fig, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(figsize=(15, 5), ncols=3, nrows=2)
    nx.draw(G0, pos=nx.spring_layout(G0), node_size=10, ax=ax1)
    ax2.matshow(nx.adjacency_matrix(G0).todense())
    ax3.matshow(S)
    nx.draw(G, pos=nx.spring_layout(G), node_size=10, ax=ax4)
    ax5.matshow(nx.adjacency_matrix(G).todense())
    ax6.matshow(costs)
    plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="Generate synthetic instances. Graphs are generated from the given model. "
                    "A instance is created by drawing similarities scores "
                    f"from N(mu={DEFAULTS['mu1']}, sigma={DEFAULTS['sigma1']}) for edges and "
                    f"from N(mu={DEFAULTS['mu2']}, sigma={DEFAULTS['sigma2']}) for non-edges. "
                    "The instance has an edge uv if s(uv) > 0 and cost(uv) = |sim(uv)|.",
        conflict_handler="resolve")

    parser.add_argument("--seed", type=int, default=0,
                        help="Seed for random sampling.")
    parser.add_argument("--dir", type=str, default=".",
                        help="Path for output dicectory.")

    model_group = parser.add_mutually_exclusive_group(required=True)
    model_group.add_argument("--barabasi-albert", action='store_true',
                             help=f"Generate a barabasi albert graph with parameters "
                                  f"m={DEFAULTS['barabasi']['m']}.")
    model_group.add_argument("--powerlaw-cluster", action='store_true',
                             help=f"Generate a powerlaw cluster graph with parameters "
                                  f"m={DEFAULTS['powerlaw_cluster']['m']}, p={DEFAULTS['powerlaw_cluster']['p']}.")
    model_group.add_argument("--duplication-divergence", action='store_true',
                             help=f"Generate a duplication divergence graph with parameters "
                                  f"p={DEFAULTS['duplication_divergence']['p']}.")

    parser.add_argument("vertices", type=int,
                        help="Number of vertices.")

    parser.add_argument("--step-size", type=int, default=1)
    parser.add_argument("--num-steps", type=int, default=1)

    namespace = parser.parse_args()

    model = ""
    if namespace.barabasi_albert:
        model = "barabasi-albert"
    elif namespace.powerlaw_cluster:
        model = "powerlaw-cluster"
    elif namespace.duplication_divergence:
        model = "duplication-divergence"

    num_vertices = [namespace.vertices + i * namespace.step_size for i in range(namespace.num_steps)]

    generate_dataset(Path(namespace.dir), model, num_vertices, seed=namespace.seed)


if __name__ == "__main__":
    main()
