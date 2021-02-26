#!/usr/bin/python3
from argparse import ArgumentParser

from pathlib import Path
import numpy as np


def convert_pace_to_input_format(input_path: Path) -> str:
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

    instance_str = f"{n} {n * (n - 1) // 2} 1\n"
    for u in range(n):
        instance_str += " ".join([f"{v + 1} {costs[u, v]}" for v in range(n) if v != u]) + "\n"
    return instance_str


if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument("input", type=str)
    parser.add_argument("output", type=str)

    args = parser.parse_args()

    txt = convert_pace_to_input_format(Path(args.input))
    with Path(args.output).open("w") as output_file:
        output_file.write(txt)

