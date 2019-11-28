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

    fpt_df["total_calls"] = fpt_df["calls"].apply(sum)
    ilp_df["total_calls"] = np.nan

    fpt_df["last_time"] = fpt_df["time"].str[-1].astype(float) / 10**9
    fpt_df["last_k"] = fpt_df["k"].str[-1].astype(float)
    fpt_df["last_calls"] = fpt_df["calls"].str[-1].astype(float)
    ilp_df["last_time"] = np.nan
    ilp_df["last_k"] = np.nan
    ilp_df["last_calls"] = np.nan

    headers = list(set(ilp_df.columns) & set(fpt_df.columns))

    df = pd.concat([ilp_df[headers], fpt_df[headers]])

    df["total_time"] = df["total_time"] / 10**9
    df.loc[df["total_time"] < 0, "total_time"] = np.nan
    df.loc[df["total_time"] > 100, "solved"] = False
    return df


def plot_solved_by_time_curve(df, output_path: Path, *, names : List[str] = None, labels : List[str] = None,
                              min_number_of_solutions: int = None, y: str = "time"):
    if min_number_of_solutions is None:
        min_number_of_solutions = 0
    if names is None:
        names = list(df["name"].unique())
    if labels is None:
        labels = names
    y_label = dict(total_time="Total Time [s]",
                   total_calls="Total Calls",
                   last_time="Time of last search step [s]",
                   last_calls="Number of calls of last search step")[y]

    d = dict()
    for name in names:
        g = df.loc[df["name"] == name]
        g = g.loc[g["solutions"].apply(lambda x: len(x[0]["edits"]) >= min_number_of_solutions if len(x) != 0 else True)]
        solved = g["solved"]
        t = pd.Series(g[y]) # .astype(float)
        t[~solved] = np.nan  # t.max() * 1.5
        d[name] = t.values

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.set_xscale("log")
    ax.grid(True)

    for name, label in zip(names, labels):
        ax.plot(np.sort(d[name]), range(len(d[name])), label=label)

    for y_max in (0, len(list(d.values())[0])):
        ax.axhline(y=y_max, c="k")
    ax.set_ylim((-50, None))

    if "time" in y:
        ax.set_xlim((10**-3, 10**2))
    ax.set_ylabel("Number of solved instances")
    ax.set_xlabel(y_label)

    ax.legend(loc="best")
    # fig.legend(loc="upper left", bbox_to_anchor=(0.9, 0.9))

    plt.savefig(output_path)
    plt.close()


def main():
    ilp_paths = list((Path.cwd() / "../../experiments/C4P4/").glob("ilp*/*.solutions.df.gzip"))
    fpt_paths = list((Path.cwd() / "../../experiments/C4P4/").glob("fpt*/*.solutions.df.gzip"))

    df = read_data(ilp_paths, fpt_paths)

    subset_df = df[df["dataset"] == "bio-C4P4-subset"]
    bio_df = df[df["dataset"] == "bio"]

    for y in ["total_time", "total_calls", "last_time", "last_calls"]:
        plot_solved_by_time_curve(bio_df, Path(f"solved-curve-ilp-vs-fpt-bio-{y}.pdf"),
                                  names=["Sparse", "MostAdjacentSubgraphs SortedGreedy Fixed"],
                                  labels=["ILP Sparse", "FPT, known $k^*$"], min_number_of_solutions=0, y=y)
        plot_solved_by_time_curve(bio_df, Path(f"solved-curve-ilp-{y}.pdf"),
                                  names=["Basic", "Single", "Sparse"],
                                  labels=["ILP", "ILP, single constraint", "ILP, sparse constraints"], min_number_of_solutions=10, y=y)
        plot_solved_by_time_curve(bio_df, Path(f"solved-curve-fpt-{y}.pdf"),
                                  names=["MostAdjacentSubgraphs SortedGreedy Exponential", "MostAdjacentSubgraphs SortedGreedy Fixed", "FirstFound Trivial IncrementByMultiplier"],
                                  labels=["FPT, exponential search strategy", "FPT, known $k^*$", "FPT, base"], min_number_of_solutions=10, y=y)
        plot_solved_by_time_curve(bio_df, Path(f"solved-curve-ilp-vs-fpt-{y}.pdf"),
                                  names=["Basic", "Sparse", "MostAdjacentSubgraphs SortedGreedy Exponential", "MostAdjacentSubgraphs SortedGreedy Fixed", "FirstFound Trivial IncrementByMultiplier"],
                                  labels=["ILP", "ILP, sparse constraints", "FPT, exponential search strategy", "FPT, known $k^*$", "FPT, base"], min_number_of_solutions=10, y=y)
        plot_solved_by_time_curve(bio_df, Path(f"solved-curve-search-strategies-{y}.pdf"),
                                  names=["MostAdjacentSubgraphs SortedGreedy Exponential", "MostAdjacentSubgraphs SortedGreedy PrunedDelta", "MostAdjacentSubgraphs SortedGreedy IncrementByMinCost", "MostAdjacentSubgraphs SortedGreedy IncrementByMultiplier", "MostAdjacentSubgraphs SortedGreedy Fixed"],
                                  labels=["Exponential", "Prune preventention", "Increment by minimum cost", "Increment by 1", "Known $k^*$"], min_number_of_solutions=10, y=y)
        plot_solved_by_time_curve(bio_df, Path(f"solved-curve-lower-bounds-fixed-{y}.pdf"),
                                  names=["MostAdjacentSubgraphs LocalSearch Fixed", "MostAdjacentSubgraphs SortedGreedy Fixed", "MostAdjacentSubgraphs Trivial Fixed"],
                                  labels=["Local search", "Greedy lower bound", "No lower bound"], min_number_of_solutions=10, y=y)
        plot_solved_by_time_curve(subset_df, Path(f"solved-curve-lower-bounds-exponential-{y}.pdf"),
                                  names=["MostAdjacentSubgraphs Greedy Exponential", "MostAdjacentSubgraphs LocalSearch Exponential", "MostAdjacentSubgraphs SortedGreedy Exponential", "MostAdjacentSubgraphs Trivial Exponential"],
                                  labels=["Simple packing", "Local search", "Greedy lower bound", "No lower bound"], min_number_of_solutions=10, y=y)
        plot_solved_by_time_curve(bio_df, Path(f"solved-curve-selectors-fixed-{y}.pdf"),
                                  names=["MostAdjacentSubgraphs SortedGreedy Fixed", "FirstFound SortedGreedy Fixed", "MostMarkedPairs SortedGreedy Fixed"],
                                  labels=["Most adjacent subgraphs", "First subgraph found", "Most marked vertex pairs"], min_number_of_solutions=10, y=y)

    # plot_lower_bound_quality(df)


if __name__ == "__main__":
    main()
