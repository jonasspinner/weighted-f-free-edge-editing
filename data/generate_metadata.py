import numpy as np
from pathlib import Path
import networkx as nx
from concurrent.futures import ThreadPoolExecutor

import yaml
yaml.add_representer(np.float64, lambda dumper, data: dumper.represent_float(data))
yaml.add_representer(np.ndarray, lambda dumper, data: dumper.represent_list(data))

from typing import Dict, Any, Tuple

from graph_io import read_metis_graph


def cost_stats(costs: np.ndarray) -> Dict[str, Any]:
    stats = dict()

    stats['mean'] = np.mean(costs)
    stats['std'] = np.std(costs)
    stats['min'] = np.min(costs)
    stats['max'] = np.max(costs)
    stats['median'] = np.quantile(costs, 0.5)
    stats['quantiles'] = np.quantile(costs, np.linspace(0, 1, 11))

    return stats


def graph_stats(G: nx.Graph) -> Dict[str, Any]:
    stats = dict()

    stats['number_of_vertices'] = G.number_of_nodes()
    stats['number_of_edges'] = G.number_of_edges()
    stats['complexity'] = G.number_of_nodes() * G.number_of_edges()
    stats['density'] = G.number_of_edges() / G.number_of_nodes() ** 2

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


def generate(paths, output_path):
    def handle_path(path):
        print(f"collecting stats on {path}")
        instance = read_metis_graph(path)
        return instance_stats(instance)

    with ThreadPoolExecutor() as executor, Path(output_path).open('w') as file:
        stats = list(executor.map(handle_path, paths))
        yaml.dump(stats, file, default_flow_style=False)


if __name__ == '__main__':
    generate(Path('bio').glob('*.metis'), Path('bio/metadata2.yaml'))
