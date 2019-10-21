import yaml
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from hashlib import sha1

plt.rcParams['axes.axisbelow'] = True


def read_data(paths) -> pd.DataFrame:

    dfs = []
    for path in paths:
        with open(path) as file:
            docs = list(yaml.safe_load_all(file))
        dfs += [pd.DataFrame(docs)]

    df = pd.concat(dfs, ignore_index=True, sort=True)

    df = pd.concat([df, df["instance"].apply(pd.Series), df["config"].apply(pd.Series)], axis=1)
    df.drop(["instance", "config"], axis=1, inplace=True)

    df = df.astype({k: "category" for k in ["commit_hash", "selector", "lower_bound", "search_strategy", "forbidden_subgraphs"]})

    df.loc[df["time"] != -1, "time"] = df.loc[df["time"] != -1, "time"] / 10**9  # convert to seconds

    return df


def plot_solved_by_time_curve(df, *, selectors=None, lower_bounds=None, search_strategies=None, min_number_of_solutions=10):
    if min_number_of_solutions is None:
        min_number_of_solutions = 0

    d = dict()
    for (k, g) in df.groupby(["selector", "lower_bound", "search_strategy"]):
        t = pd.Series(g.loc[g["solutions"].apply(lambda x: len(x[0]["edits"]) >= min_number_of_solutions if len(x) != 0 else True), "time"])
        t[t == -1] = t.max() * 1.5
        d[k] = t.values

    fig, ax = plt.subplots(figsize=(15, 10))
    ax.set_xscale("log")
    ax.grid(True)

    for k in d:
        if all(b is None or a in b for a, b in zip(k, [selectors, lower_bounds, search_strategies])):
            ax.plot(np.sort(d[k]), range(len(d[k])), label="{0} {1} {2}".format(*k))

    ax.axhline(y=len(list(d.values())[0]), c="black")
    ax.set_ylim((-50, None))
    ax.set_xlim((None, 100))
    ax.set_ylabel("Number of solved instances")
    ax.set_xlabel("Total Time [s]")

    fig.legend(loc="lower right")
    plt.show()


def plot_lower_bound_quality(df):
    fig, ax = plt.subplots()
    ax.grid(True)

    df = df[~df["stats"].isnull() & (df["search_strategy"] != "Fixed")]

    fmts = iter(list(".1234+x_") + [4, 5])

    for lb, lb_group in df.groupby("lower_bound"):
        lb_group = lb_group.sort_values(by="solution_cost")
        init_k = lb_group["stats"].apply(lambda s: s["k"][0] if len(s["k"]) != 0 else 0)
        init_time = lb_group["stats"].apply(lambda s: s["time"][0] if "time" in s and len(s["time"]) != 0 else 0)
        num_vertices = lb_group["name"].str.replace(".*/", "").apply(lambda x: int(x[:-6].split("-")[-1]))

        ax.scatter(init_time, 1 - init_k / lb_group["solution_cost"], label=lb, marker=next(fmts), s=100)

    fig.legend()
    fig.tight_layout()
    plt.show()


def main():
    solution_paths = (Path.home() / "experiments/experiments/C4P4/").glob("fpt*/all.solutions.yaml")
    solution_paths = list(solution_paths)
    h = sha1("#".join(str(path.absolute()) for path in sorted(solution_paths)).encode("utf8")).hexdigest()[:10]

    df_path = Path(f"{h}.df")
    if df_path.exists():
        df = pd.read_pickle(str(df_path))
    else:
        df = read_data(solution_paths)
        df.to_pickle(str(df_path))

    plot_solved_by_time_curve(df, search_strategies=["IncrementByMultiplier"], min_number_of_solutions=0)

    plot_lower_bound_quality(df)


if __name__ == "__main__":
    main()
