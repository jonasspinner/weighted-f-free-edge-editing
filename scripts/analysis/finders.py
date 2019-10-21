import yaml
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import argparse
from hashlib import sha1

from typing import Sequence

plt.rcParams['axes.axisbelow'] = True


def read_data(benchmark_paths: Sequence[Path], meta_path: Path) -> pd.DataFrame:
    dfs = []
    for path in benchmark_paths:
        with path.open() as file:
            dfs += [pd.DataFrame(yaml.safe_load_all(file))]
    df = pd.concat(dfs, ignore_index=True)

    df = df.astype({k: "category" for k in ["commit_hash", "finder", "forbidden_subgraphs", "type"]})
    df["time_mean"] = df["time_mean"].astype(float)
    df = pd.concat([df.drop(["instance"], axis=1), df["instance"].apply(pd.Series)], axis=1)

    df["name"] = df["name"].str.replace(".*/", "")
    df.rename(columns={"time": "raw_time"}, inplace=True)

    df["time_mean"] = df["raw_time"].apply(lambda x: np.nan if len(x) == 0 else np.mean(x[1:])) / 10**9
    df["time_std"] = df["raw_time"].apply(lambda x: np.nan if len(x) == 0 else np.std(x[1:])) / 10**9

    with meta_path.open() as file:
        meta_df = pd.DataFrame(yaml.safe_load(file))
    meta_df = pd.concat([meta_df.drop(["costs", "graph"], axis=1), meta_df["graph"].apply(pd.Series)], axis=1)
    meta_df = meta_df.drop(["connected_components"], axis=1)

    df = df.join(meta_df.set_index("name"), on="name")

    return df


def plot_time_to_find_all(df, *, forbidden_subgraphs="C4P4", plot_timeout=False):
    fsg_df = df[(df["forbidden_subgraphs"] == forbidden_subgraphs) & (df["type"] == "find_all_subgraphs")].copy()
    fig, (ax1, ax2) = plt.subplots(ncols=2, sharey="all", figsize=(10, 5))

    ax1.grid(True)
    ax2.grid(True)

    if plot_timeout:
        not_finished = fsg_df["time_mean"].max() * 1.1
        fsg_df.loc[fsg_df["time_mean"].isnull(), "time_mean"] = not_finished
        fsg_df.loc[fsg_df["time_std"].isnull(), "time_std"] = 0

        d = df.groupby("name")["count"].max()
        fsg_df.loc[fsg_df["count"] == -1, "count"] = fsg_df.loc[fsg_df["count"] == -1, "name"].apply(lambda n: d[n])

        ax1.axhline(y=not_finished, c="black", zorder=-1)
        ax2.axhline(y=not_finished, c="black", zorder=-1)

    for finder, group_df in fsg_df.groupby("finder"):
        if len(group_df) == 0: continue

        x1 = group_df["number_of_vertices"]
        x2 = group_df["count"]
        y = group_df["time_mean"]
        yerr = group_df["time_std"]

        # ax1.errorbar(x1, y, yerr=yerr, fmt="o", label=finder)
        ax1.scatter(x1, y, marker="o", label=finder, zorder=1)
        ax2.scatter(x2, y, marker="o", label=finder, zorder=1)

    ax1.set_xlabel("Number of vertices")
    ax2.set_xlabel("Number of forbidden subgraphs")
    ax1.set_ylabel("Time [s]")
    ax2.legend(loc="best", fancybox=False)  # frameon=False
    fig.tight_layout()
    plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="Generate a subset of a dataset.",
        conflict_handler="resolve")

    parser.add_argument("--benchmarks-dir", type=str,
                        help="Path for input directory.")
    parser.add_argument("--dataset", type=str, default="bio",
                        help="Name of dataset.")
    parser.add_argument("--data-dir", type=str, help="Path to datasets.")
    parser.add_argument("--save-df", type=bool, default=True)

    options = parser.parse_args()

    input_dir = Path.home() / "experiments/experiments/"
    if options.benchmarks_dir is not None:
        input_dir = Path(options.benchmarks_dir)

    meta_path = Path(f"../../data/{options.dataset}/{options.dataset}.metadata.yaml")

    benchmarks_paths = input_dir.glob(f"finder*/{options.dataset}.benchmarks.yaml")
    benchmarks_paths = list(benchmarks_paths)
    h = sha1("#".join(str(path.absolute()) for path in sorted(benchmarks_paths)).encode("utf8")).hexdigest()[:10]

    df_path = Path(f"{h}.df")
    if df_path.exists():
        df = pd.read_pickle(str(df_path))
    else:
        df = read_data(benchmarks_paths, meta_path)
        if options.save_df:
            df.to_pickle(str(df_path))

    plot_time_to_find_all(df, forbidden_subgraphs="C4P4")
    plot_time_to_find_all(df, forbidden_subgraphs="P3")


if __name__ == '__main__':
    main()
