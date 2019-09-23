import numpy as np
from pathlib import Path
import networkx as nx
from concurrent.futures import ThreadPoolExecutor
import argparse

import yaml
yaml.add_representer(np.float64, lambda dumper, data: dumper.represent_float(data))
yaml.add_representer(np.ndarray, lambda dumper, data: dumper.represent_list(data))

from typing import Dict, Any, Tuple, Sequence

from graph_io import read_metis_graph


def cost_stats(costs: np.ndarray) -> Dict[str, Any]:
    stats = dict()

    stats['mean'] = np.mean(costs)
    stats['std'] = np.std(costs)
    stats['min'] = np.min(costs)
    stats['max'] = np.max(costs)
    stats['median'] = np.median(costs)
    stats['quantiles'] = np.percentile(costs, 100 * np.linspace(0, 1, 11))

    return stats


def graph_stats(G: nx.Graph) -> Dict[str, Any]:
    stats = dict()

    n, m = G.number_of_nodes(), G.number_of_edges()

    stats['number_of_vertices'] = n
    stats['number_of_edges'] = m
    stats['complexity'] = n * m
    stats['density'] = 2 * m / (n * (n - 1))

    stats['connected_components'] = []
    for G_hat in (G.subgraph(c) for c in nx.connected_components(G)):
        component_stats = dict()

        component_stats['number_of_vertices'] = G_hat.number_of_nodes()
        component_stats['number_of_edges'] = G_hat.number_of_edges()
        component_stats['diameter'] = nx.diameter(G_hat, usebounds=True)
        component_stats['radius'] = nx.radius(G_hat, usebounds=True)
        component_stats['center_size'] = len(nx.center(G_hat, usebounds=True))
        component_stats['periphery_size'] = len(nx.periphery(G_hat, usebounds=True))

        stats['connected_components'] += [component_stats]

    stats['number_of_connected_components'] = len(stats['connected_components'])

    stats['average_clustering_coefficient'] = nx.average_clustering(G)

    return stats


def instance_stats(instance: Tuple[nx.Graph, np.ndarray]) -> Dict[str, Any]:
    G, S = instance
    stats = dict()

    stats['name'] = G.name
    stats['graph'] = graph_stats(G)
    stats['costs'] = cost_stats(S[np.triu_indices(S.shape[0], 1)])

    return stats


def compute_metadata(paths: Sequence[Path], output_path: Path) -> None:
    def handle_path(path):
        print(f"collecting stats on {path}")
        instance = read_metis_graph(path)
        return instance_stats(instance)

    with ThreadPoolExecutor() as executor, Path(output_path).open('w') as file:
        stats = list(executor.map(handle_path, paths))
        yaml.dump(stats, file, default_flow_style=False)


def main():
    parser = argparse.ArgumentParser(
            description="Generates metadata for instances in a given directory.")

    parser.add_argument("dir", type=str, default=".",
                        help="Path for input directory.")

    parser.add_argument("--pattern", type=str, default="*.graph",
                        help="Pattern for files in the directory.")

    options = parser.parse_args()

    input_dir = Path(options.dir)
    compute_metadata(list(input_dir.glob(options.pattern)), input_dir / f"{input_dir.resolve().name}.metadata.yaml")


if __name__ == '__main__':
    main()
