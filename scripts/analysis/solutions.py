import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from typing import List, Optional

plt.rcParams['axes.axisbelow'] = True


def read_data(ilp_paths, fpt_paths) -> pd.DataFrame:
    ilp_df = pd.concat(map(pd.read_pickle, ilp_paths))
    fpt_df = pd.concat(map(pd.read_pickle, fpt_paths))

    ilp_df["name"] = "Basic"
    ilp_df.loc[ilp_df["single_constraints"], "name"] = "Single"
    ilp_df.loc[ilp_df["sparse_constraints"], "name"] = "Sparse"

    fpt_df["name"] = fpt_df.apply(lambda row: f"{row['selector']} {row['lower_bound']} {row['search_strategy']}",
                                  axis=1)

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
    # df.loc[df["total_time"] > df["timelimit"], "solved"] = False
    return df


def plot_solved_by_time_curve(df, output_path: Path, *, names: List[str] = None, labels: List[str] = None,
                              min_number_of_edits: int = None, y: str = "time", y_max: Optional[int] = None,
                              title: Optional[str] = None):
    if min_number_of_edits is None:
        min_number_of_edits = 0
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
        print(f"{name} {len(g)} {g['solved'].sum()}")
        g = g.loc[g["solutions"].apply(lambda x: len(x[0]["edits"]) >= min_number_of_edits if len(x) != 0 else True)]
        solved = g["solved"]
        t = pd.Series(g[y])  # .astype(float)
        t[~solved] = np.nan  # t.max() * 1.5
        d[name] = t.values
        print(f"{name} {len(g)} {len(g[solved])} {solved.sum()} {len(d[name])}")

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.set_xscale("log")
    ax.grid(True)

    print(names)
    for name, label in zip(names, labels):
        ax.plot(np.sort(d[name]), range(len(d[name])), label=label)

    if y_max is None:
        # y_max = max([c.max() for _, c in df.groupby("name").count().iterrows()])
        y_max = max([v.shape[0] for v in d.values()])
    for hline_y in (0, y_max):
        ax.axhline(y=hline_y, c="k")
    ax.set_ylim((-50, None))

    if "time" in y:
        ax.set_xlim((10**-3, 10**3))
    ax.set_ylabel("Number of solved instances")
    ax.set_xlabel(y_label)

    ax.legend(loc="best")
    # fig.legend(loc="upper left", bbox_to_anchor=(0.9, 0.9))

    if title is not None:
        ax.set_title(title)

    plt.savefig(output_path)
    plt.close()


def main():
    ilp_paths = list((Path.cwd() / "../../experiments/C4P4/").glob("ilp*/*.solutions.df.gzip"))
    fpt_paths = list((Path.cwd() / "../../experiments/C4P4/").glob("fpt*/*.solutions.df.gzip"))

    print("Ignore {}".format([x for x in fpt_paths if "all=0" in str(x)]))
    fpt_paths = [x for x in fpt_paths if "all=0" not in str(x)]

    df = read_data(ilp_paths, fpt_paths)

    # Filter datasets
    bio_df = df[df["dataset"] == "bio"]
    subset_df = df[df["dataset"] == "bio-C4P4-subset"]
    unweighted_df = df[df["dataset"] == "bio-unweighted"]

    subset_instances = set(subset_df["instance"])
    subset_df = pd.concat([subset_df, bio_df[bio_df["instance"].apply(lambda i: i in subset_instances)]])
    # subset_df = subset_df.drop_duplicates(subset=subset_df.dtypes[subset_df.dtypes != "object"].index.values)
    subset_df = subset_df.drop_duplicates(subset=["instance", "name"])

    # plot_solved_by_time_curve(bio_df, Path(f"solved-curve-ilp-vs-fpt-geq-0.pdf"),
    #                           names=["Sparse", "MostAdjacentSubgraphs SortedGreedy Exponential"],
    #                           labels=["ILP (sparse)", "FPT (greedy, most adjacent, exponential)"],
    #                           min_number_of_edits=0, y="total_time")
    # plot_solved_by_time_curve(bio_df, Path(f"solved-curve-ilp-vs-fpt-geq-10.pdf"),
    #                           names=["Sparse", "MostAdjacentSubgraphs SortedGreedy Exponential"],
    #                           labels=["ILP (sparse)", "FPT (greedy, most adjacent, exponential)"],
    #                           min_number_of_edits=10, y="total_time")
    #
    # plot_solved_by_time_curve(bio_df, Path(f"solved-curve-search-strategies-total_time.pdf"),
    #                           names=["MostAdjacentSubgraphs SortedGreedy Exponential",
    #                                  "MostAdjacentSubgraphs SortedGreedy PrunedDelta",
    #                                  "MostAdjacentSubgraphs SortedGreedy IncrementByMinCost",
    #                                  "MostAdjacentSubgraphs SortedGreedy IncrementByMultiplier"],
    #                           labels=["Exponential", "Prune preventention", "Increment by minimum cost",
    #                                   "Increment by 1"],
    #                           min_number_of_edits=10, y="total_time")
    #
    # plot_solved_by_time_curve(bio_df, Path(f"solved-curve-search-strategies-last_time.pdf"),
    #                           names=["MostAdjacentSubgraphs SortedGreedy Exponential",
    #                                  "MostAdjacentSubgraphs SortedGreedy PrunedDelta",
    #                                  "MostAdjacentSubgraphs SortedGreedy IncrementByMinCost",
    #                                  "MostAdjacentSubgraphs SortedGreedy IncrementByMultiplier"],
    #                           labels=["Exponential", "Prune preventention", "Increment by minimum cost",
    #                                   "Increment by 1"],
    #                           min_number_of_edits=10, y="last_time")
    #
    # plot_solved_by_time_curve(bio_df, Path(f"solved-curve-.pdf"),
    #                           names=["MostAdjacentSubgraphs SortedGreedy Exponential",
    #                                  "MostAdjacentSubgraphs SortedGreedy PrunedDelta",
    #                                  "MostAdjacentSubgraphs SortedGreedy IncrementByMinCost",
    #                                  "MostAdjacentSubgraphs SortedGreedy IncrementByMultiplier"],
    #                           labels=["Exponential", "Prune preventention", "Increment by minimum cost",
    #                                   "Increment by 1"],
    #                           min_number_of_edits=10, y="last_time")

    # No. 06
    n06a_names = [f"MostAdjacentSubgraphs {lb} Fixed"
                  for lb in ["LocalSearch", "SortedGreedy", "Trivial", "LPRelaxation", "NPS_MWIS_Solver",
                             "LSSWZ_MWIS_Solver"]]
    n06a_labels = ["Local search", "Greedy lower bound", "No lower bound", "Relaxation", "NPS MWIS", "LSSWZ MWIS"]

    n06b_names = [f"MostAdjacentSubgraphs {lb} Exponential"
                  for lb in ["Greedy", "LocalSearch", "SortedGreedy", "Trivial", "LPRelaxation"]]
    n06b_labels = ["Simple packing", "Local search", "Greedy lower bound", "No lower bound", "Relaxation"]

    for col in ["total_time", "total_calls"]:
        plot_solved_by_time_curve(subset_df, Path(f"06-fixed-{col}.pdf"), names=n06a_names, labels=n06a_labels,
                                  min_number_of_edits=10, y=col,
                                  title="Relaxation vs. Packing, weighted bio")
        plot_solved_by_time_curve(subset_df, Path(f"06-exponential-{col}.pdf"), names=n06b_names, labels=n06b_labels,
                                  min_number_of_edits=10, y=col)

    # No. 11
    n11_names = [f"MostAdjacentSubgraphs {lb} Fixed"
                 for lb in ["LocalSearch", "LPRelaxation", "NPS_MWIS_Solver", "LSSWZ_MWIS_Solver"]]
    n11_labels = ["Local search", "Relaxation", "NPS MWIS", "LSSWZ MWIS"]

    plot_solved_by_time_curve(unweighted_df, Path(f"11-fixed-total_time.pdf"), names=n11_names, labels=n11_labels,
                              min_number_of_edits=10, y="total_time",
                              title="Relaxation vs. Packing, unweighted bio")
    plot_solved_by_time_curve(unweighted_df, Path(f"11-fixed-total_calls.pdf"), names=n11_names, labels=n11_labels,
                              min_number_of_edits=10, y="total_calls",
                              title="Relaxation vs. Packing, unweighted bio")

    # No. 12
    # Weighted Single Edge Editing vs. Most Adjacent Subgraphs
    n12_names = [f"{sel} {lb} Fixed"
                 for sel in ["MostAdjacentSubgraphs", "SingleEdgeEditing"]
                 for lb in ["LocalSearch", "SortedGreedy", "LPRelaxation", "NPS_MWIS_Solver"]]
    n12_labels = [f"{lb}, {sel}"
                  for sel in ["Subgraph", "Single Edge"]
                  for lb in ["Local search", "Greedy", "Relaxation", "NPS MWIS"]]

    plot_solved_by_time_curve(subset_df, Path(f"12-fixed-total_time.pdf"), names=n12_names, labels=n12_labels,
                              min_number_of_edits=10, y="total_time", y_max=1058,
                              title="Single Edge Editing vs. Most Adjacent Subgraphs, weighted bio")
    plot_solved_by_time_curve(subset_df, Path(f"12-fixed-total_calls.pdf"), names=n12_names, labels=n12_labels,
                              min_number_of_edits=10, y="total_calls", y_max=1058,
                              title="Single Edge Editing vs. Most Adjacent Subgraphs, weighted bio")

    for y in []:  # ["total_time", "total_calls", "last_time", "last_calls"]:
        plot_solved_by_time_curve(bio_df, Path(f"solved-curve-ilp-vs-fpt-bio-{y}.pdf"),
                                  names=["Sparse", "MostAdjacentSubgraphs SortedGreedy Fixed"],
                                  labels=["ILP Sparse", "FPT, known $k^*$"], min_number_of_edits=0, y=y)
        plot_solved_by_time_curve(bio_df, Path(f"solved-curve-ilp-{y}.pdf"),
                                  names=["Basic", "Single", "Sparse"],
                                  labels=["ILP", "ILP, single constraint", "ILP, sparse constraints"],
                                  min_number_of_edits=10, y=y)
        plot_solved_by_time_curve(bio_df, Path(f"solved-curve-fpt-{y}.pdf"),
                                  names=["MostAdjacentSubgraphs SortedGreedy Exponential",
                                         "MostAdjacentSubgraphs SortedGreedy Fixed",
                                         "FirstFound Trivial IncrementByMultiplier"],
                                  labels=["FPT, exponential search strategy", "FPT, known $k^*$", "FPT, base"],
                                  min_number_of_edits=10, y=y)
        plot_solved_by_time_curve(bio_df, Path(f"solved-curve-ilp-vs-fpt-{y}.pdf"),
                                  names=["Basic", "Sparse",
                                         "MostAdjacentSubgraphs SortedGreedy Exponential",
                                         "MostAdjacentSubgraphs SortedGreedy Fixed",
                                         "FirstFound Trivial IncrementByMultiplier"],
                                  labels=["ILP", "ILP, sparse constraints", "FPT, exponential search strategy",
                                          "FPT, known $k^*$", "FPT, base"],
                                  min_number_of_edits=10, y=y)
        plot_solved_by_time_curve(bio_df, Path(f"solved-curve-search-strategies-{y}.pdf"),
                                  names=["MostAdjacentSubgraphs SortedGreedy Exponential",
                                         "MostAdjacentSubgraphs SortedGreedy PrunedDelta",
                                         "MostAdjacentSubgraphs SortedGreedy IncrementByMinCost",
                                         "MostAdjacentSubgraphs SortedGreedy IncrementByMultiplier",
                                         "MostAdjacentSubgraphs SortedGreedy Fixed"],
                                  labels=["Exponential", "Prune preventention", "Increment by minimum cost",
                                          "Increment by 1", "Known $k^*$"],
                                  min_number_of_edits=10, y=y)
        plot_solved_by_time_curve(bio_df, Path(f"solved-curve-lower-bounds-fixed-{y}.pdf"),
                                  names=["MostAdjacentSubgraphs LocalSearch Fixed",
                                         "MostAdjacentSubgraphs SortedGreedy Fixed",
                                         "MostAdjacentSubgraphs Trivial Fixed"],
                                  labels=["Local search", "Greedy lower bound", "No lower bound"],
                                  min_number_of_edits=10, y=y)
        plot_solved_by_time_curve(subset_df, Path(f"solved-curve-lower-bounds-exponential-{y}.pdf"),
                                  names=["MostAdjacentSubgraphs Greedy Exponential",
                                         "MostAdjacentSubgraphs LocalSearch Exponential",
                                         "MostAdjacentSubgraphs SortedGreedy Exponential",
                                         "MostAdjacentSubgraphs Trivial Exponential"],
                                  labels=["Simple packing", "Local search", "Greedy lower bound", "No lower bound"],
                                  min_number_of_edits=10, y=y)
        plot_solved_by_time_curve(bio_df, Path(f"solved-curve-selectors-fixed-{y}.pdf"),
                                  names=["MostAdjacentSubgraphs SortedGreedy Fixed", "FirstFound SortedGreedy Fixed",
                                         "MostMarkedPairs SortedGreedy Fixed"],
                                  labels=["Most adjacent subgraphs", "First subgraph found",
                                          "Most marked vertex pairs"],
                                  min_number_of_edits=10, y=y)

    # plot_lower_bound_quality(df)


if __name__ == "__main__":
    main()
