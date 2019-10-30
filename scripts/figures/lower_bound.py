import numpy as np
import matplotlib.pyplot as plt
import subprocess
import yaml
from itertools import product
from pathlib import Path

from typing import List, Dict, Any


def run(path: Path, selectors: List[str] = None, lower_bounds: List[str] = None, search_strategies: List[str] = None, timelimit=10):
    if search_strategies is None:
        search_strategies = ["Exponential"]
    if lower_bounds is None:
        lower_bounds = ["SortedGreedy"]
    if selectors is None:
        selectors = ["MostAdjacentSubgraphs"]

    for selector, lower_bound, search_strategy in product(selectors, lower_bounds, search_strategies):
        try:
            out = subprocess.run([
                "../../cmake-build-release/fpt",
                "--input", path,
                "--search-strategy", search_strategy, "--all", "1",
                "--permutation", str(0),
                "--multiplier", str(100),
                "--selector", selector,
                "--lower-bound", lower_bound,
                "--F", "C4P4",
                "--verbosity", str(0),
                "--timelimit", str(timelimit)], capture_output=True)  # timeout=4*timelimit
            doc = yaml.safe_load(out.stdout.decode())
            if doc is None:
                print(out)
                continue
            yield doc
        except subprocess.TimeoutExpired:
            print(f"timeout {selector}, {lower_bound}, {search_strategy}")
            continue
        except RuntimeError as e:
            print(e)


def plot_experiment_docs(docs: List[Dict[str, Any]], output_path: Path) -> None:
    doc = docs[0]
    name = doc["instance"]["name"].split("/")[-1]

    k_final = doc["solution_cost"]

    fig, ax = plt.subplots(figsize=(15, 5))
    ax.grid(True)
    ax.set_yscale("log")

    for doc in docs:
        k = np.array(doc["stats"]["k"])
        calls = np.array(doc["stats"]["calls"])
        time = np.array(doc["stats"]["time"]) / 10**9
        config = doc["config"]

        x = k
        y = time

        ax.plot(x, y, "o", label="{0} {1}".format(config["selector"], config["lower_bound"], config["search_strategy"]))
        if k_final != -1:
            ax.axvline(x=k_final, c="black")

    ax.set_title(name)
    ax.set_xlim((0, None))
    ax.set_ylabel("Time [s]")
    ax.set_xlabel("Editing cost $k$")
    #ax.set_ylim((1, None))
    ax.legend()
    plt.savefig(output_path)
    plt.show()


def main():
    paths = [Path("../../data/bio/bio-nr-1590-size-56.graph")]

    for path in paths:
        print(path)
        docs = []
        selectors = ["MostAdjacentSubgraphs"]  # + ["MostMarkedPairs"]
        lower_bounds = ["SortedGreedy", "LocalSearch", "Greedy"] + ["Trivial"]
        search_strategies = ["IncrementByMultiplier"]  # + ["Exponential"]  # + ["PrunedDelta"]
        for doc in run(path, selectors=selectors, lower_bounds=lower_bounds, search_strategies=search_strategies):
            docs += [doc]

        if len(docs) > 0:
            plot_experiment_docs(docs, Path(path.stem + ".pdf"))


if __name__ == '__main__':
    main()
