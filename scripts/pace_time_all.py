#!/usr/bin/env python3
from pathlib import Path
from subprocess import run, TimeoutExpired
from tempfile import NamedTemporaryFile
from pace_to_graph import convert_pace_to_input_format
import yaml


def main():
    timeout = 60

    instances = list(Path("../data/pace2021-exact").glob("*.gr"))

    print("Graph,k,Total Time [s],Solved")
    for instance in sorted(instances, key=lambda x: x.name):
        try:
            with NamedTemporaryFile() as input_file:
                txt = convert_pace_to_input_format(instance)
                input_file.write(txt.encode("utf8"))
                input_file.flush()
                out = run(["../cmake-build-release/fpt", "--input", input_file.name,
                           "--F", "P3", "--multiplier", "1",
                           "--lower-bound", "WeightedPackingLocalSearch", "--selector", "MostAdjacentSubgraphs",
                           "--search-strategy", "Exponential",
                           "--all", "1", "--verbosity", "0", "--timelimit", "-1"], timeout=timeout, capture_output=True)
                d = yaml.safe_load(out.stdout.decode("utf8"))

            time = d["time"] / 10 ** 9
            cost = d["solution_cost"]
            solved = True
        except TimeoutExpired:
            time = timeout
            cost = ""
            solved = False
        print(f"{instance.name},{cost},{time},{solved}")


if __name__ == '__main__':
    try:
        main()
    except KeyboardInterrupt:
        pass
