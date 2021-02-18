import argparse
from subprocess import run
from tempfile import NamedTemporaryFile
from pathlib import Path
from os import environ
from typing import Tuple, Dict, List
import numpy as np

BINARY_ENV_NAME = "PEACE_WEIGHTED_CLUSTER_EDITING_BINARY"


def convert_pace_format_file(input_path: Path) -> str:
    with input_path.open() as input_file:
        lines = input_file.readlines()

    p_str, problem, n, m = lines[0].split()
    assert p_str == "p"
    assert problem == "cep"
    n, m = int(n), int(m)
    edges = [tuple(line.split()) for line in lines[1:]]
    edges = [(int(u) - 1, int(v) - 1) for u, v in edges]

    costs = -np.ones((n, n), dtype=float)
    for u, v in edges:
        costs[u, v] = 1
        costs[v, u] = 1

    instance_str = f"{n}\n"
    instance_str += "\n".join([f"{u}" for u in range(n)]) + "\n"
    for u in range(n - 1):
        instance_str += "\t".join(map(str, costs[u, u + 1:])) + "\n"
    return instance_str


def peace(input_path: Path, binary_path: Path) -> Tuple[Dict, List[List[str]]]:
    pace_format_instance = convert_pace_format_file(input_path)

    with NamedTemporaryFile() as output_file, NamedTemporaryFile() as input_file:
        input_file.write(pace_format_instance.encode("utf8"))
        input_file.flush()

        out = run([binary_path, "--mode", "2", input_file.name, output_file.name], capture_output=True)

        if out.returncode != 0:
            raise RuntimeError(
                f"command failed with error code {out.returncode} and stderr '{out.stderr.decode('utf-8')}'.")

        meta_info = out.stdout.decode("utf8")
        solution = output_file.read().decode("utf8")

    statistic_lines = meta_info[meta_info.find("STATISTICS\n") + 11:].split("\n")[:4]
    statistics = dict([tuple(line.split(": ")) for line in statistic_lines])
    solution_lines = solution.split("\n")[:-1]
    component_lines = solution_lines[:-1]
    components = [line.split(": ")[1].split(" ")[:-1] for line in component_lines]

    return statistics, components


def parse_args() -> Tuple[Path, Path]:
    parser = argparse.ArgumentParser(
        description="Execute PEACE cluster editing.")

    parser.add_argument("--binary", type=str, default=None,
                        help="Path for PEACE binary.")

    parser.add_argument("--instance", type=str, default="../../data/pace2021-exact/exact031.gr",
                        help="Path to instance in PACE format.")

    options = parser.parse_args()

    binary_path_str = options.binary if options.binary is not None else environ.get(BINARY_ENV_NAME)
    if binary_path_str is None:
        raise RuntimeError(
            f"binary must be specified as argument or environment variable {BINARY_ENV_NAME} must be specified.")
    binary_path = Path(binary_path_str)
    if not binary_path.exists():
        raise RuntimeError(f"binary {binary_path} does not exist.")

    instance_path = Path(options.instance)
    if not instance_path.exists():
        raise RuntimeError(f"instance {instance_path} does not exist.")

    return instance_path, binary_path


def main():
    instance_path, binary_path = parse_args()

    statistics, components = peace(instance_path, binary_path)

    print(f"cost = {statistics['costs']}")
    print(f"time = {statistics['time']}")
    print("components = ")
    for component in components:
        print(f"\t{' '.join(component)}")


if __name__ == '__main__':
    main()
