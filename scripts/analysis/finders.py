import yaml
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import argparse

from typing import Sequence

plt.rcParams['axes.axisbelow'] = True


def read_data(benchmark_paths: Sequence[Path], meta_path: Path) -> pd.DataFrame:
    df = pd.concat(map(pd.read_pickle, benchmark_paths), ignore_index=True)

    # TODO: remove
    df.loc[df["forbidden_subgraphs"].isnull(), "forbidden_subgraphs"] = "C5P5"

    with meta_path.open() as file:
        meta_df = pd.DataFrame(yaml.safe_load(file))
    meta_df = pd.concat([meta_df.drop(["costs", "graph"], axis=1), meta_df["graph"].apply(pd.Series)], axis=1)
    meta_df = meta_df.drop(["connected_components"], axis=1)

    df = df.join(meta_df.rename(columns={"name": "instance"}).set_index("instance"), on="instance")

    return df


def plot_time_to_find_all(df: pd.DataFrame, output_dir: Path, *, forbidden_subgraphs: str = "C4P4",
                          plot_timeout: bool = False, max_number_of_vertices: int = None):
    if max_number_of_vertices is None:
        max_number_of_vertices = np.inf
    fsg_df = df[
        (df["forbidden_subgraphs"] == forbidden_subgraphs) &
        (df["finder_benchmark_type"] == "find_all_subgraphs") &
        (df["number_of_vertices"] <= max_number_of_vertices)].copy()
    fsg_df["timedout"] = fsg_df["time_mean"].isnull()

    fig, (ax1, ax2) = plt.subplots(ncols=2, sharey="all", figsize=(8, 4))

    ax1.grid(True)
    ax2.grid(True)
    ax2.ticklabel_format(axis="x", scilimits=(3, 5))

    if plot_timeout:
        not_finished = fsg_df["time_mean"].max() * 1.1
        fsg_df.loc[fsg_df["time_mean"].isnull(), "time_mean"] = not_finished
        fsg_df.loc[fsg_df["time_std"].isnull(), "time_std"] = 0

        d = df.groupby("instance")["count"].max()
        fsg_df.loc[fsg_df["count"] == -1, "count"] = fsg_df.loc[fsg_df["count"] == -1, "instance"].apply(lambda n: d[n])

        ax1.axhline(y=not_finished, c="black", zorder=-1)
        ax2.axhline(y=not_finished, c="black", zorder=-1)

    s = 10
    colors = [f"C{i}" for i in range(10)]

    for color, (finder, group_df) in zip(colors, [x for x in fsg_df.groupby("finder") if len(x[1]) != 0]):

        x1 = group_df["number_of_vertices"]
        x2 = group_df["count"]
        y = group_df["time_mean"]

        ax1.scatter(x1, y, label=finder, s=s, c=color)
        ax2.scatter(x2, y, label=finder, s=s, c=color)

        for ax, x, d in zip([ax1, ax2], [x1, x2], [(int(forbidden_subgraphs[-1]),), (1,)]):
            p, *_ = np.linalg.lstsq(np.vstack([x[~group_df["timedout"]]**i for i in d]).T, y[~group_df["timedout"]], rcond=None)

            x_pred = np.linspace(x.min(), 2 * x.max(), 20)
            y_pred = np.vstack([x_pred**i for i in d]).T @ p

            ax.plot(x_pred, y_pred, "k--", alpha=0.25, c=color)

    ax1.set_xlabel("Number of vertices")
    ax2.set_xlabel("Number of forbidden subgraphs")
    ax1.set_ylabel("Time [s]")

    for ax in (ax1, ax2):
        y = fsg_df["time_mean"]
        eps = y.max() / 20
        ax.set_ylim((0, y.max() + eps))
    for ax, col in zip((ax1, ax2), ("number_of_vertices", "count")):
        x = fsg_df[col]
        eps = (x.max() - x.min()) / 20
        ax.set_xlim((x.min() - eps, x.max() + eps))

    ax2.legend(loc="best", fancybox=False)  # frameon=False
    fig.tight_layout()

    plt.savefig(output_dir / f"finder-benchmark-{forbidden_subgraphs}-n_max={max_number_of_vertices}.pdf")
    # plt.show()


def main():
    parser = argparse.ArgumentParser(
        description="Generate a subset of a dataset.",
        conflict_handler="resolve")

    parser.add_argument("--benchmarks-dir", type=str,
                        help="Path for input directory.")
    parser.add_argument("--dataset", type=str, default="bio",
                        help="Name of dataset.")
    parser.add_argument("--data-dir", type=str, help="Path to datasets.")

    options = parser.parse_args()

    input_dir = Path.cwd() / "../../experiments/"
    if options.benchmarks_dir is not None:
        input_dir = Path(options.benchmarks_dir)

    data_dir = Path.cwd() / ".." / ".." / "data"
    if options.data_dir is not None:
        data_dir = Path(options.data_dir)

    output_dir = Path.cwd()

    meta_path = data_dir / f"{options.dataset}/{options.dataset}.metadata.yaml"

    benchmarks_paths = list(input_dir.glob(f"finder*/{options.dataset}.benchmarks.df.gzip"))

    df = read_data(benchmarks_paths, meta_path)

    for fsg in df["forbidden_subgraphs"].unique():
        for n in [100, 800]:
            if (df["forbidden_subgraphs"] == fsg).any():
                plot_time_to_find_all(df, output_dir, forbidden_subgraphs=fsg, plot_timeout=False,
                                      max_number_of_vertices=n)


if __name__ == '__main__':
    main()
