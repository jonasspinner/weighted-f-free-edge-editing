import argparse
import yaml
from pathlib import Path


from graph_io import write_metis_graph, read_metis_graph


def convert_instance_to_unweighted(input_path: Path, output_path: Path):
    graph, scores = read_metis_graph(input_path)

    origin_comment = f"unweighted instance based on {graph.name}"

    write_metis_graph(output_path, graph, comments=[origin_comment])


def generate_dataset_unweighted(config_path: Path, input_dir: Path, output_dir: Path):
    with config_path.open() as file:
        config = yaml.safe_load(file)

    weighted_dataset: str = config["weighted_dataset"]

    output_dir.mkdir()
    for input_path in (input_dir / weighted_dataset).glob("*.graph"):
        output_path = output_dir / f"{input_path.stem}-unweighted.graph"
        convert_instance_to_unweighted(input_path, output_path)


def main():
    parser = argparse.ArgumentParser(
        description="Generate an unweighted version of a dataset.",
        conflict_handler="resolve")

    parser.add_argument("config", type=str,
                        help="Config file specifying the dataset.")
    parser.add_argument("--input-dir", type=str,
                        help="Path for input directory. Default is the directory of the config file.")
    parser.add_argument("--output-dir", type=str,
                        help="Path for output directory. "
                             "Default is the directory with the same name as the config file.")

    options = parser.parse_args()

    config_path = Path(options.config)

    input_dir = config_path.parent
    if options.input_dir is not None:
        input_dir = Path(options.input_dir)

    output_dir = config_path.parent / config_path.stem
    if options.output_dir is not None:
        output_dir = Path(options.output_dir)

    generate_dataset_unweighted(config_path, input_dir, output_dir)


if __name__ == "__main__":
    main()