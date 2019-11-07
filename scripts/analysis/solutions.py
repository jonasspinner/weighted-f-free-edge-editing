import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from typing import List

plt.rcParams['axes.axisbelow'] = True


def read_data(ilp_paths, fpt_paths) -> pd.DataFrame:
    ilp_df = pd.concat(map(pd.read_pickle, ilp_paths))
    fpt_df = pd.concat(map(pd.read_pickle, fpt_paths))

    ilp_df["name"] = "Basic"
    ilp_df.loc[ilp_df["single_constraints"], "name"] = "Single"
    ilp_df.loc[ilp_df["sparse_constraints"], "name"] = "Sparse"

    fpt_df["name"] = fpt_df.apply(lambda row: f"{row['selector']} {row['lower_bound']} {row['search_strategy']}", axis=1)

    headers = list(set(ilp_df.columns) & set(fpt_df.columns))

    df = pd.concat([ilp_df[headers], fpt_df[headers]])
    return df


def plot_solved_by_time_curve(df, output_path: Path, *, names : List[str] = None, labels : List[str] = None,
                              min_number_of_solutions: int = None):
    if min_number_of_solutions is None:
        min_number_of_solutions = 0
    if names is None:
        names = list(df["name"].unique())
    if labels is None:
        labels = names

    d = dict()
    for name in names:
        g = df.loc[df["name"] == name]
        g = g.loc[g["solutions"].apply(lambda x: len(x[0]["edits"]) >= min_number_of_solutions if len(x) != 0 else True)]
        solved = g["solution_cost"] != -1
        t = pd.Series(g["total_time"])
        t[~solved] = t.max() * 1.5
        d[name] = t.values

    fig, ax = plt.subplots(figsize=(6, 4))
    ax.set_xscale("log")
    ax.grid(True)

    for name, label in zip(names, labels):
        ax.plot(np.sort(d[name]) / 10**9, range(len(d[name])), label=label)

    for y in (0, len(list(d.values())[0])):
        ax.axhline(y=y, c="darkgrey")
    ax.set_ylim((-50, None))
    ax.set_xlim((10**-3, 10**2))
    ax.set_ylabel("Number of solved instances")
    ax.set_xlabel("Total Time [s]")

    ax.legend(loc="upper left")
    # fig.legend(loc="upper left", bbox_to_anchor=(0.9, 0.9))
    plt.savefig(output_path)


def main():
    ilp_paths = list((Path.cwd() / "../../experiments/C4P4/").glob("ilp*/*.solutions.df.gzip"))
    fpt_paths = list((Path.cwd() / "../../experiments/C4P4/").glob("fpt*/*.solutions.df.gzip"))

    df = read_data(ilp_paths, fpt_paths)

    subset_df = df[df["dataset"] == "bio-C4P4-subset"]
    bio_df = df[df["dataset"] == "bio"]

    plot_solved_by_time_curve(subset_df, Path("solved-curve-ilp-vs-fpt.pdf"),
                              names=["Basic", "Single", "Sparse", "MostAdjacentSubgraphs SortedGreedy Exponential"],
                              labels=["ILP", "ILP Single", "ILP Sparse", "FPT"], min_number_of_solutions=10)
    plot_solved_by_time_curve(bio_df, Path("solved-curve-search-strategies.pdf"),
                              names=["MostAdjacentSubgraphs SortedGreedy Exponential", "MostAdjacentSubgraphs SortedGreedy PrunedDelta", "MostAdjacentSubgraphs SortedGreedy IncrementByMinCost", "MostAdjacentSubgraphs SortedGreedy IncrementByMultiplier"],
                              labels=["Exponential", "PrunedDelta", "Increment by mininum cost", "Increment by 1"], min_number_of_solutions=10)
    plot_solved_by_time_curve(subset_df, Path("solved-curve-lower-bounds.pdf"),
                              names=["MostAdjacentSubgraphs Greedy Exponential", "MostAdjacentSubgraphs LocalSearch Exponential", "MostAdjacentSubgraphs SortedGreedy Exponential", "MostAdjacentSubgraphs Trivial Exponential"],
                              labels=["Greedy", "LocalSearch", "SortedGreedy", "No lower bound"], min_number_of_solutions=10)
    plot_solved_by_time_curve(subset_df, Path("solved-curve-selectors.pdf"),
                              names=["MostAdjacentSubgraphs SortedGreedy Exponential", "FirstFound SortedGreedy Exponential", "MostMarkedPairs SortedGreedy Exponential"],
                              labels=["MostAdjacentSubgraphs", "FirstFound", "MostMarkedPairs"], min_number_of_solutions=10)

    # plot_lower_bound_quality(df)


if __name__ == "__main__":
    main()
