import yaml
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
import argparse

from typing import Sequence, Tuple, Optional

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


BBox = Tuple[Tuple[float, float], Tuple[float, float]]


def plot_time_to_find_all(df: pd.DataFrame, output_dir: Path, *, forbidden_subgraphs: str = "C4P4",
                          plot_timeout: bool = False, bboxes: Tuple[Tuple[BBox, BBox], Tuple[BBox, BBox]],
                          finders: Sequence[str] = [], labels: Sequence[str] = []):
    fsg_df = df[
        (df["forbidden_subgraphs"] == forbidden_subgraphs) &
        (df["finder_benchmark_type"] == "find_all_subgraphs")].copy()
    fsg_df["timedout"] = fsg_df["time_mean"].isnull()

    fig, axes = plt.subplots(nrows=2, ncols=2, figsize=(10, 6))

    for (ax1, ax2), (bbox_n, bbox_count) in zip(axes, bboxes):

        for ax in (ax1, ax2):
            ax.grid(True)
            #ax.ticklabel_format(axis="y", scilimits=(-1, 1))
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

        for color, finder, label in zip(colors, finders, labels):
            group_df = fsg_df[fsg_df["finder"] == finder]

            x1 = group_df["number_of_vertices"]
            x2 = group_df["count"]
            y = group_df["time_mean"]

            ax1.scatter(x1, y, label=label, s=s, c=color)
            ax2.scatter(x2, y, label=label, s=s, c=color)

            for ax, x, d, bbox in zip([ax1, ax2], [x1, x2], [(int(forbidden_subgraphs[-1]),), (1,)], [bbox_n, bbox_count]):
                x_train, y_train = x[~group_df["timedout"] & (x < 1.5 * bbox[0][1])], y[~group_df["timedout"] & (x < 1.5 * bbox[0][1])]
                p, *_ = np.linalg.lstsq(np.vstack([x_train**i for i in d]).T, y_train, rcond=None)

                x_pred = np.linspace(0, 1.5 * bbox[0][1], 40)
                y_pred = np.vstack([x_pred**i for i in d]).T @ p

                ax.plot(x_pred, y_pred, "k--", alpha=0.25, c=color)

            print(f"{finder} {np.isnan(y).sum()} / {y.shape[0]}")

        ax1.set_xlabel("Number of vertices")
        ax2.set_xlabel("Number of forbidden subgraphs")
        ax1.set_ylabel("Time [s]")

        def expand(lim: Tuple[float, float], epsilon: float = 1/20) -> Optional[Tuple[float, float]]:
            if lim is None: return None
            if None in lim: return lim
            d = (lim[1] - lim[0]) * epsilon
            return lim[0] - d, lim[1] + d

        for ax, bbox in [(ax1, bbox_n), (ax2, bbox_count)]:
            if bbox is not None:
                xlim, ylim = bbox
                ax.set_xlim(expand(xlim))
                ax.set_ylim(expand(ylim))

    axes[0][1].legend(loc="best", fancybox=False)  # frameon=False
    fig.tight_layout()

    plt.savefig(output_dir / f"finder-benchmark-{forbidden_subgraphs}.pdf")
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

    for fsg, bboxes, finders, labels in [
        ("P3",   ((((0, 600), (0, 0.25)),  ((0, 0.5 * 10**7), (0, 0.25))),
                  (((0, 200), (0, 0.025)), ((0, 0.4 * 10**6), (0, 0.025)))), ["NaiveP3", "OuterP3", "CenterP3"], ["Naive", "Fill from outer vertices", "Edge expansion"]),
        ("C4P4", ((((0, 600), (0, 10)),    ((0, 0.4 * 10**8), (0, 10))),
                  (((0, 200), (0, 0.15)),  ((0, 0.6 * 10**7), (0, 0.5)))),   ["NaiveC4P4", "CenterC4P4", "EndpointRecC4P4"], ["Naive", "Midpoint", "Endpoint"]),
        ("C5P5", ((((0, 600), (0, 10)),    ((0, 2 * 10**8),   (0, 10))),
                  (((0, 200), (0, 0.5)),   ((0, 0.2 * 10**6), (0, 0.5)))),   ["NaiveRecC5P5", "CenterRecC5P5", "EndpointRecC5P5"], ["Naive", "Midpoint", "Endpoint"])
    ]:
        plot_time_to_find_all(df, output_dir, forbidden_subgraphs=fsg, plot_timeout=False,
                              bboxes=bboxes, finders=finders, labels=labels)


if __name__ == '__main__':
    main()
