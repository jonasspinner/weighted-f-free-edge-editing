from glob import glob
from pathlib import Path
import subprocess
import yaml


BIO_GRAPHS = [path.stem for path in Path("data/bio").glob("*.graph")]
BIO_GRAPHS_C4P4_SUBSET = [path.stem for path in Path("data/bio-C4P4-subset").glob("*.graph")]
BIO_GRAPHS_SUBSET_A = [path.stem for path in Path("data/bio-subset-A").glob("*.graph")]


def get_dataset_files(name):
    if name == "bio":
        return BIO_GRAPHS
    elif name == "bio-C4P4-subset":
        return BIO_GRAPHS_C4P4_SUBSET
    elif name == "bio-subset-A":
        return BIO_GRAPHS_SUBSET_A
    else:
        return []

MULTIPLIER = [100]
PERMUTATION = [0]
TIMELIMITS = [100]
FSG = ["P3", "C4P4"]
ILP_CONSTRAINTS = ["basic", "sparse", "single"]
ILP_NUM_THREADS = [1]

FPT_SELECTOR = ["FirstFound", "MostMarkedPairs", "MostAdjacentSubgraphs"] # + ["LeastWeight"]
FPT_LOWER_BOUND = ["Trivial", "LocalSearch", "SortedGreedy"] + ["Greedy"] # + ["LPRelaxation"]
FPT_PRE_MARK = [0]
FPT_SEARCH_STRATEGY = ["IncrementByMultiplier"] + ["PrunedDelta", "IncrementByMinCost"] + ["Exponential"]

FINDERS = ["CenterRecC4P4", "CenterRecP3", "EndpointRecC4P4", "EndpointRecP3", "CenterC4P4", "CenterP3"] + ["NaiveC4P4", "NaiveP3"] + ["CenterRecC5P5", "EndpointRecC5P5"]


PRELIM_TIMELIMITS = [100]
PRELIM_ILP_CONSTRAINTS = ILP_CONSTRAINTS
PRELIM_FPT_LOWER_BOUND = ["Trivial", "LocalSearch", "SortedGreedy"]
PRELIM_FPT_FPT_SELECTOR = ["FirstFound", "MostMarkedPairs", "MostAdjacentSubgraphs"]
PRELIM_FPT_SEARCH_STRATEGY = ["IncrementByMultiplier"]


rule all:
        input:
                expand("experiments/finder-benchmark.finder={finder}/bio.benchmarks.yaml", finder=FINDERS),

                expand("experiments/C4P4/ilp.timelimit={timelimit}.threads={threads}.constraints={constraints}/bio-C4P4-subset.solutions.yaml",
                       timelimit=TIMELIMITS, threads=ILP_NUM_THREADS, constraints=ILP_CONSTRAINTS),
                expand("experiments/C4P4/fpt.timelimit={timelimit}.selector={selector}.lower-bound={lower_bound}.all=1.pre-mark={pre_mark}.search-strategy={search_strategy}/bio-C4P4-subset.solutions.yaml",
                       timelimit=TIMELIMITS, selector=FPT_SELECTOR, lower_bound=FPT_LOWER_BOUND, pre_mark=FPT_PRE_MARK, search_strategy=["Exponential"]),

                # expand("experiments/{fsg}/fpt.timelimit={timelimit}.selector={selector}.lower-bound={lower_bound}.all=1.pre-mark={pre_mark}.search-strategy={search_strategy}/bio.solutions.yaml",
                #        fsg=FSG, timelimit=[1000], selector=["MostAdjacentSubgraphs"], lower_bound=["SortedGreedy"], pre_mark=FPT_PRE_MARK, search_strategy=["PrunedDelta"]),
                # expand("experiments/{fsg}/ilp.timelimit={timelimit}.threads={threads}.constraints={constraints}/bio.solutions.yaml",
                #        fsg=FSG, timelimit=[1000], threads=ILP_NUM_THREADS, constraints=["sparse"]),

                expand("experiments/calls-experiment/C4P4/fpt.timelimit={timelimit}.selector={selector}.lower-bound={lower_bound}.all=1.pre-mark={pre_mark}.search-strategy={search_strategy}/bio-subset-A.yaml",
                       timelimit=[10], selector=["MostAdjacentSubgraphs"], lower_bound=FPT_LOWER_BOUND, pre_mark=[0], search_strategy=["IncrementByMultiplier"]),

                expand("data/{dataset}/{dataset}.metadata.yaml", dataset=["bio", "bio-C4P4-subset", "bio-subset-A"])


rule preliminary:
        input:
                expand("experiments/{fsg}/ilp.timelimit={timelimit}.threads={threads}.constraints={constraints}/bio.solutions.yaml",
                       fsg=FSG, timelimit=PRELIM_TIMELIMITS, threads=[1], constraints=PRELIM_ILP_CONSTRAINTS),
                expand("experiments/{fsg}/fpt.timelimit={timelimit}.selector={selector}.lower-bound={lower_bound}.all=0.pre-mark={pre_mark}.search-strategy=Fixed/bio.solutions.yaml",
                       fsg=FSG, timelimit=PRELIM_TIMELIMITS, selector=PRELIM_FPT_FPT_SELECTOR, lower_bound=PRELIM_FPT_LOWER_BOUND, pre_mark=FPT_PRE_MARK),
                expand("experiments/{fsg}/fpt.timelimit={timelimit}.selector={selector}.lower-bound={lower_bound}.all=1.pre-mark={pre_mark}.search-strategy={search_strategy}/bio.solutions.yaml",
                       fsg=FSG, timelimit=PRELIM_TIMELIMITS, selector=PRELIM_FPT_FPT_SELECTOR, lower_bound=PRELIM_FPT_LOWER_BOUND, pre_mark=FPT_PRE_MARK, search_strategy=PRELIM_FPT_SEARCH_STRATEGY)
        output:
                "experiments/preliminary_rule"
        shell: "touch {output}"


rule ilp:
        input:
                "data/{dataset}/{graph}.graph"
        output:
                "experiments/{fsg}/ilp.timelimit={timelimit}.threads={threads}.constraints={constraints}/{dataset}/{graph}.{multiplier}.{permutation}.solution.yaml"
        params:
                hard_timeout = lambda wildcards, output: int(1.1 * int(wildcards.timelimit))
        run:
                constraint_args = dict(basic=[], sparse=["--sparse-constraints", "1"], single=["--single-constraints", "1"])
                try:
                    subprocess.run(f"cmake-build-release/ilp "
                                   f"--input {input} --permutation {wildcards.permutation} --multiplier {wildcards.multiplier} --F {wildcards.fsg} --output {output} "
                                   f"--timelimit {wildcards.timelimit} "
                                   f"--num-threads {wildcards.threads}".split(" ") + constraint_args[wildcards.constraints], timeout=params.hard_timeout)
                except subprocess.TimeoutExpired:
                    pass

rule collect_ilp:
        input:
                lambda wildcards: expand("experiments/{{fsg}}/"
                                          "ilp.timelimit={{timelimit}}.threads={{threads}}.constraints={{constraints}}/"
                                          "{{dataset}}/{graph}.{multiplier}.{permutation}.solution.yaml",
                       graph=get_dataset_files(wildcards.dataset), multiplier=MULTIPLIER, permutation=PERMUTATION)
        output:
                "experiments/{fsg}/ilp.timelimit={timelimit}.threads={threads}.constraints={constraints}/{dataset}.solutions.yaml"
        run:
                with open(output[0], "w") as out_file:
                    for path in input:
                        with open(path) as in_file:
                            out_file.write(in_file.read())


rule copy_instance_solution:
        input:
                f"experiments/{{fsg}}/ilp.timelimit={PRELIM_TIMELIMITS[0]}.threads={ILP_NUM_THREADS[0]}.constraints={ILP_CONSTRAINTS[1]}/{{dataset}}/{{graph}}.{{multiplier}}.{PERMUTATION[0]}.solution.yaml"
        output:
                "experiments/{fsg}/solutions/{dataset}/{graph}.{multiplier}.solution.yaml"
        shell:
                "cp {input} {output}"


rule fpt_fixed_k_from_solution:
        input:
                instance = "data/{dataset}/{graph}.graph",
                solution = "experiments/{fsg}/solutions/{dataset}/{graph}.{multiplier}.solution.yaml"
        output:
                "experiments/{fsg}/fpt.timelimit={timelimit}.selector={selector}.lower-bound={lower_bound}.all=0.pre-mark={pre_mark}.search-strategy=Fixed/{dataset}/{graph}.{multiplier}.{permutation}.solution.yaml"
        params:
                hard_timeout = lambda wildcards, output: int(1.1 * int(wildcards.timelimit))
        run:
                k = yaml.safe_load(open(input.solution))["solution_cost"]
                try:
                    subprocess.run(f"cmake-build-release/fpt "
                                   f"--input {input} --permutation {wildcards.permutation} --multiplier {wildcards.multiplier} --F {wildcards.fsg} --output {output} "
                                   f"--selector {wildcards.selector} --lower-bound {wildcards.lower_bound} "
                                   f"--all 0 --pre-mark {wildcards.pre_mark} "
                                   f"--search-strategy Fixed --k {k}".split(" "), timeout=params.hard_timeout)
                except subprocess.TimeoutExpired:
                    pass

rule fpt:
        input:
                "data/{dataset}/{graph}.graph"
        output:
                "experiments/{fsg}/fpt.timelimit={timelimit}.selector={selector}.lower-bound={lower_bound}.all=1.pre-mark={pre_mark}.search-strategy={search_strategy}/{dataset}/{graph}.{multiplier}.{permutation}.solution.yaml"
        params:
                hard_timeout = lambda wildcards, output: int(4 * int(wildcards.timelimit))
        run:
                try:
                    subprocess.run(f"cmake-build-release/fpt "
                                   f"--input {input} --permutation {wildcards.permutation} --multiplier {wildcards.multiplier} --F {wildcards.fsg} --output {output} "
                                   f"--selector {wildcards.selector} --lower-bound {wildcards.lower_bound} "
                                   f"--all 1 --pre-mark {wildcards.pre_mark} "
                                   f"--search-strategy {wildcards.search_strategy} "
                                   f"--timelimit {wildcards.timelimit}".split(" "), timeout=params.hard_timeout)
                except subprocess.TimeoutExpired:
                    pass


rule collect_fpt_fixed:
        input:
                lambda wildcards: expand("experiments/{{fsg}}/"
                                         "fpt.timelimit={{timelimit}}.selector={{selector}}.lower-bound={{lower_bound}}.all=0.pre-mark={{pre_mark}}.search-strategy=Fixed/"
                                         "{{dataset}}/{graph}.{multiplier}.{permutation}.solution.yaml",
                       graph=get_dataset_files(wildcards.dataset), multiplier=MULTIPLIER, permutation=PERMUTATION)
        output:
                "experiments/{fsg}/fpt.timelimit={timelimit}.selector={selector}.lower-bound={lower_bound}.all=0.pre-mark={pre_mark}.search-strategy=Fixed/{dataset}.solutions.yaml"
        run:
                with open(output[0], "w") as out_file:
                    for path in input:
                        with open(path) as in_file:
                            out_file.write(in_file.read())

rule collect_fpt:
        input:
                lambda wildcards: expand("experiments/{{fsg}}/"
                                         "fpt.timelimit={{timelimit}}.selector={{selector}}.lower-bound={{lower_bound}}.all=1.pre-mark={{pre_mark}}.search-strategy={{search_strategy}}/"
                                         "{{dataset}}/{graph}.{multiplier}.{permutation}.solution.yaml",
                       graph=get_dataset_files(wildcards.dataset), multiplier=MULTIPLIER, permutation=PERMUTATION)
        output:
                "experiments/{fsg}/fpt.timelimit={timelimit}.selector={selector}.lower-bound={lower_bound}.all=1.pre-mark={pre_mark}.search-strategy={search_strategy}/{dataset}.solutions.yaml"
        run:
                with open(output[0], "w") as out_file:
                    for path in input:
                        with open(path) as in_file:
                            out_file.write(in_file.read())

rule finder_experiment:
        input:
                instance = "data/{dataset}/{graph}.graph"
        output:
                "experiments/finder-benchmark.finder={finder}/{dataset}/{graph}.{permutation}.benchmark.yaml"
        params:
                iterations = 10,
                hard_timeout = 100
        run:
                try:
                    subprocess.run(f"cmake-build-release/finder_benchmark --input {input.instance} --output {output} --finder {wildcards.finder} --iterations {params.iterations}".split(" "), timeout=params.hard_timeout)
                except subprocess.TimeoutExpired:
                    pass

rule collect_finder_experiment:
        input:
                lambda wildcards: expand("experiments/"
                                         "finder-benchmark.finder={{finder}}/"
                                         "{{dataset}}/{graph}.{permutation}.benchmark.yaml", graph=get_dataset_files(wildcards.dataset), permutation=PERMUTATION)
        output:
                "experiments/finder-benchmark.finder={finder}/{dataset}.benchmarks.yaml"
        run:
                with open(output[0], "w") as out_file:
                    for path in input:
                        with open(path) as in_file:
                            out_file.write(in_file.read())


rule metadata:
        input:
                "data/{dataset}"
        output:
                "data/{dataset}/{dataset}.metadata.yaml"
        shell:
                "python3 scripts/compute_metadata.py {input}"


rule calls_experiment:
        input:
                "data/{dataset}/{graph}.graph"
        output:
                "experiments/calls-experiment/{fsg}/fpt.timelimit={timelimit}.selector={selector}.lower-bound={lower_bound}.all=1.pre-mark={pre_mark}.search-strategy={search_strategy}/{dataset}/{graph}.{multiplier}.{permutation}.result.yaml"
        run:
                subprocess.run(f"cmake-build-release/fpt "
                               f"--input {input} --permutation {wildcards.permutation} --multiplier {wildcards.multiplier} --F {wildcards.fsg} --output {output} "
                               f"--selector {wildcards.selector} --lower-bound {wildcards.lower_bound} "
                               f"--all 1 --pre-mark {wildcards.pre_mark} "
                               f"--search-strategy {wildcards.search_strategy} "
                               f"--timelimit {wildcards.timelimit}".split(" "))

rule collect_calls_experiment:
        input:
                lambda wildcards: expand("experiments/calls-experiment/{{fsg}}/fpt.timelimit={{timelimit}}.selector={{selector}}.lower-bound={{lower_bound}}.all=1.pre-mark={{pre_mark}}.search-strategy={{search_strategy}}/"
                                          "{{dataset}}/{graph}.{multiplier}.{permutation}.result.yaml",
                                          graph=get_dataset_files(wildcards.dataset), multiplier=MULTIPLIER, permutation=PERMUTATION)
        output:
                "experiments/calls-experiment/{fsg}/fpt.timelimit={timelimit}.selector={selector}.lower-bound={lower_bound}.all=1.pre-mark={pre_mark}.search-strategy={search_strategy}/{dataset}.results.yaml"
        run:
                with open(output[0], "w") as out_file:
                    for path in input:
                        with open(path) as in_file:
                            out_file.write(in_file.read())