import argparse
import yaml
from pathlib import Path
import shutil


def generate_dataset_subset(config_path: Path, input_dir: Path, output_dir: Path):
    with config_path.open() as file:
        config = yaml.safe_load(file)

    superset = config["superset"]
    instances = config["instances"]

    if output_dir.is_dir() and all((output_dir / instance).is_file() for instance in instances):
        return

    output_dir.mkdir()
    for instance in instances:
        shutil.copy(str(input_dir / superset / instance), str(output_dir))


def main():
    parser = argparse.ArgumentParser(
        description="Generate a subset of a dataset.",
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

    generate_dataset_subset(config_path, input_dir, output_dir)


if __name__ == "__main__":
    main()