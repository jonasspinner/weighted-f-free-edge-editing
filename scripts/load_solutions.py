import pandas as pd
import numpy as np
import yaml
from pathlib import Path
from pandas.api.types import CategoricalDtype
from collections import defaultdict

from typing import List, Dict, Any, Callable, Union

DEFAULT_CATEGORY: Callable[[], Union[str, CategoricalDtype]] = lambda: "category"
CATEGORIES = defaultdict(DEFAULT_CATEGORY,
                         forbidden_subgraphs=CategoricalDtype([
                             "P3", "P4", "P5", "P6", "C4P4", "C5P5", "C6P6", ", C4_C5_2K2", "C4_C5_P5_Bowtie_Necktie"]),
                         selector=CategoricalDtype([
                             "FirstFound", "LeastWeight", "MostMarkedPairs", "MostAdjacentSubgraphs",
                             "SingleEdgeEditing"]),
                         lower_bound=CategoricalDtype([
                             "Trivial", "Greedy", "SortedGreedy", "LocalSearch", "LPRelaxation", "NPS_MWIS_Solver",
                             "LSSWZ_MWIS_Solver"]),
                         search_strategy=CategoricalDtype([
                             "IncrementByMultiplier", "PrunedDelta", "IncrementByMinCost", "Exponential", "Fixed"]),
                         dataset=CategoricalDtype([
                             "barabasi-albert", "bio", "bio-C4P4-subset", "bio-subset-A", "duplication-divergence",
                             "misc", "powerlaw-cluster", "bio-subset-B", "bio-unweighted"]),
                         finder=CategoricalDtype(
                             ["OuterP3", "CenterP3", "CenterRecP3", "EndpointRecP3", "NaiveRecP3", "NaiveP3"] +
                             ["CenterC4P4", "CenterRecC4P4", "EndpointRecC4P4", "NaiveRecC4P4", "NaiveC4P4"] +
                             ["CenterRecC5P5", "EndpointRecC5P5", "NaiveRecC5P5"] +
                             ["CenterRecC6P6", "EndpointRecC6P6", "NaiveRecC6P6"]),
                         finder_benchmark_type=CategoricalDtype([
                             "find_all_subgraphs", "find_one_subgraph"])
                         )


def load_data(path: Path) -> List[Dict[str, Any]]:
    with path.open() as file:
        docs = list(yaml.safe_load_all(file))
    return docs


def build_fpt_dataframe(docs: List[Dict[str, Any]]) -> pd.DataFrame:
    df = pd.DataFrame(docs)

    df["solved"] = df["solutions"].apply(lambda x: len(x) > 0)
    df.rename(columns={"time": "total_time"}, inplace=True)

    df = pd.concat([df.drop(["config", "instance", "stats"], axis=1)] + [df[col].apply(pd.Series) for col in
                                                                         ["config", "instance", "stats"]], axis=1)

    df[["dataset", "instance"]] = df["name"].str.split("/", expand=True)[[1, 2]]

    df = df.astype({k: CATEGORIES[k] for k in
                    ["commit_hash", "forbidden_subgraphs", "selector", "lower_bound", "search_strategy", "dataset"]})

    df["n"] = df["instance"].str[:-6].str.split("-").str[4].astype(int)

    return df


def build_ilp_dataframe(docs: List[Dict[str, Any]]) -> pd.DataFrame:
    df = pd.DataFrame(docs)

    df["solved"] = df["solutions"].apply(lambda x: len(x) > 0)
    df.rename(columns={"time": "total_time"}, inplace=True)

    df = pd.concat([df.drop(["config", "instance"], axis=1)] +
                   [df[col].apply(pd.Series) for col in ["config", "instance"]], axis=1)

    df[["dataset", "instance"]] = df["name"].str.split("/", expand=True)[[1, 2]]

    df = df.astype({k: CATEGORIES[k] for k in ["commit_hash", "forbidden_subgraphs", "dataset"]})

    df["n"] = df["instance"].str[:-6].str.split("-").str[-1].astype(int)

    return df


def build_finder_dataframe(docs: List[Dict[str, Any]]) -> pd.DataFrame:
    df = pd.DataFrame(docs)

    df.rename(columns={"time": "raw_time", "type": "finder_benchmark_type"}, inplace=True)

    df = pd.concat([df.drop(["instance"], axis=1), df["instance"].apply(pd.Series)], axis=1)

    df[["dataset", "instance"]] = df["name"].str.split("/", expand=True)[[1, 2]]

    df = df.astype({k: CATEGORIES[k]
                    for k in ["commit_hash", "finder", "forbidden_subgraphs", "finder_benchmark_type"]})
    df["time_mean"] = df["time_mean"].astype(float)

    df["time_mean"] = df["raw_time"].apply(lambda x: np.nan if len(x) == 0 else np.mean(x[1:]) / 10**9)
    df["time_std"] = df["raw_time"].apply(lambda x: np.nan if len(x) == 0 else np.std(x[1:]) / 10**9)

    df["n"] = df["instance"].str[:-6].str.split("-").str[-1].astype(int)

    return df


def save_dataframe(path: Path, df: pd.DataFrame) -> None:
    df.to_pickle(str(path))


def convert_docs_to_dataframes(input_experiments_path: Path, output_experiments_path: Path, pattern: str,
                               build_dataframe: Callable[[List[Dict[str, Any]]], pd.DataFrame]) -> None:
    paths = [p.relative_to(input_experiments_path) for p in input_experiments_path.glob(pattern)]

    for path in paths:
        input_path = input_experiments_path / path
        output_path = output_experiments_path / path.with_name(path.name.replace(".yaml", ".df.gzip"))

        if output_path.exists():
            print(f"file already exists {path}")
            continue

        output_path.parent.mkdir(parents=True, exist_ok=True)

        docs = load_data(input_path)
        df = build_dataframe(docs)
        save_dataframe(output_path, df)

        print(f"converted {path}")


def convert_fpt_solutions_to_dataframes(input_experiments_path: Path, output_experiments_path: Path,
                                        pattern: str = "*/fpt*/*.solutions.yaml"):
    convert_docs_to_dataframes(input_experiments_path, output_experiments_path,
                               pattern, build_fpt_dataframe)


def convert_ilp_solutions_to_dataframes(input_experiments_path: Path, output_experiments_path: Path,
                                        pattern: str = "*/ilp*/*.solutions.yaml"):
    convert_docs_to_dataframes(input_experiments_path, output_experiments_path,
                               pattern, build_ilp_dataframe)


def convert_finder_benchmarks_to_dataframes(input_experiments_path: Path, output_experiments_path: Path,
                                            pattern: str = "finder*/*.benchmarks.yaml"):
    convert_docs_to_dataframes(input_experiments_path, output_experiments_path,
                               pattern, build_finder_dataframe)


def main() -> None:
    input_experiments_path = Path.home() / "ba" / "experiments"
    output_experiments_path = Path.cwd() / ".." / "experiments"

    convert_fpt_solutions_to_dataframes(input_experiments_path, output_experiments_path)
    convert_fpt_solutions_to_dataframes(input_experiments_path / "calls-experiment",
                                        output_experiments_path / "calls-experiment",
                                        "*/fpt*/*.results.yaml")
    convert_ilp_solutions_to_dataframes(input_experiments_path, output_experiments_path)
    convert_finder_benchmarks_to_dataframes(input_experiments_path, output_experiments_path)


if __name__ == '__main__':
    main()
