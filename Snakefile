from glob import glob
import subprocess
import yaml

BIO_GRAPHS = [path.split("/")[-1][:-6] for path in glob("data/bio/*.graph")]
# BIO_GRAPHS = BIO_GRAPHS[:10]
MULTIPLIER = [100]
PERMUTATION = [0]
TIMELIMITS = [100]
FSG = ["C4P4"]
ILP_CONSTRAINTS = ["basic", "sparse", "single"]
ILP_NUM_THREADS = [1]

FPT_SELECTOR = ["FirstFound", "MostMarkedPairs", "MostAdjacentSubgraphs"] # + ["LeastWeight"]
FPT_LOWER_BOUND = ["Trivial", "LocalSearch", "SortedGreedy"] # + ["Greedy", "LPRelaxation"]
FPT_PRE_MARK = [0]
FPT_SEARCH_STRATEGY = ["IncrementByMultiplier"] # + ["PrunedDelta", "Exponential", "IncrementByMinCost"]

FINDERS = ["CenterRecC4P4", "CenterRecP3", "EndpointRecC4P4", "EndpointRecP3", "CenterC4P4", "CenterP3"] + ["NaiveC4P4", "NaiveP3"] # + ["CenterRecC5P5", "EndpointRecC5P5"]


rule all:
        input:
                expand("experiments/{fsg}/ilp.timelimit={timelimit}.threads={threads}.constraints={constraints}/bio/{graph}.{multiplier}.{permutation}.solution.yaml",
                       fsg=FSG, timelimit=TIMELIMITS, threads=ILP_NUM_THREADS, constraints=ILP_CONSTRAINTS, graph=BIO_GRAPHS, multiplier=MULTIPLIER, permutation=PERMUTATION),
                expand("experiments/finder-benchmark.finder={finder}/all.benchmark.yaml", finder=FINDERS),
                expand("experiments/{fsg}/fpt.timelimit={timelimit}.selector={selector}.lower-bound={lower_bound}.all=0.pre-mark={pre_mark}.search-strategy=Fixed/all.solutions.yaml",
                       fsg=FSG, timelimit=TIMELIMITS, selector=FPT_SELECTOR, lower_bound=FPT_LOWER_BOUND, pre_mark=FPT_PRE_MARK),
                expand("experiments/{fsg}/fpt.timelimit={timelimit}.selector={selector}.lower-bound={lower_bound}.all=1.pre-mark={pre_mark}.search-strategy={search_strategy}/all.solutions.yaml",
                       fsg=FSG, timelimit=TIMELIMITS, selector=FPT_SELECTOR, lower_bound=FPT_LOWER_BOUND, pre_mark=FPT_PRE_MARK, search_strategy=FPT_SEARCH_STRATEGY),
                "data/bio/bio.metadata.yaml"

rule ilp:
        input:
                "data/{dataset}/{graph}.graph"
        output:
                "experiments/{fsg}/ilp.timelimit={timelimit}.threads={threads}.constraints={constraints}/{dataset}/{graph}.{multiplier}.{permutation}.solution.yaml"
        params:
                hard_timeout=lambda wildcards, output: int(1.1 * int(wildcards.timelimit))
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
                expand("experiments/{{fsg}}/ilp.timelimit={{timelimit}}.threads={{threads}}.constraints={{constraints}}/bio/{graph}.{multiplier}.{permutation}.solution.yaml",
                       graph=BIO_GRAPHS, multiplier=MULTIPLIER, permutation=PERMUTATION)
        output:
                "experiments/{fsg}/ilp.timelimit={timelimit}.threads={threads}.constraints={constraints}/all.solutions.yaml"
        run:
                with open(output[0], "w") as out_file:
                    for path in input:
                        with open(path) as in_file:
                            out_file.write(in_file.read())


rule copy_instance_solution:
        input:
                f"experiments/{{fsg}}/ilp.timelimit={TIMELIMITS[0]}.threads={ILP_NUM_THREADS[0]}.constraints={ILP_CONSTRAINTS[0]}/{{dataset}}/{{graph}}.{{multiplier}}.{PERMUTATION[0]}.solution.yaml"
        output:
                "experiments/solutions/{fsg}/{dataset}/{graph}.{multiplier}.solution.yaml"
        shell:
                "cp {input} {output}"


rule fpt_fixed_k_from_solution:
        input:
                instance = "data/{dataset}/{graph}.graph",
                solution = "experiments/solutions/{fsg}/{dataset}/{graph}.{multiplier}.solution.yaml"
        output:
                "experiments/{fsg}/fpt.timelimit={timelimit}.selector={selector}.lower-bound={lower_bound}.all=0.pre-mark={pre_mark}.search-strategy=Fixed/{dataset}/{graph}.{multiplier}.{permutation}.solution.yaml"
        params:
                hard_timeout=lambda wildcards, output: int(1.1 * int(wildcards.timelimit))
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
                hard_timeout=lambda wildcards, output: int(1.1 * int(wildcards.timelimit))
        run:
                try:
                    subprocess.run(f"cmake-build-release/fpt "
                                   f"--input {input} --permutation {wildcards.permutation} --multiplier {wildcards.multiplier} --F {wildcards.fsg} --output {output} "
                                   f"--selector {wildcards.selector} --lower-bound {wildcards.lower_bound} "
                                   f"--all 1 --pre-mark {wildcards.pre_mark} "
                                   f"--search-strategy {wildcards.search_strategy}".split(" "), timeout=params.hard_timeout)
                except subprocess.TimeoutExpired:
                    pass


rule collect_fpt_fixed:
        input:
                expand("experiments/{{fsg}}/fpt.timelimit={{timelimit}}.selector={{selector}}.lower-bound={{lower_bound}}.all=0.pre-mark={{pre_mark}}.search-strategy=Fixed/bio/{graph}.{multiplier}.{permutation}.solution.yaml",
                       graph=BIO_GRAPHS, multiplier=MULTIPLIER, permutation=PERMUTATION)
        output:
                "experiments/{fsg}/fpt.timelimit={timelimit}.selector={selector}.lower-bound={lower_bound}.all=0.pre-mark={pre_mark}.search-strategy=Fixed/all.solutions.yaml"
        run:
                with open(output[0], "w") as out_file:
                    for path in input:
                        with open(path) as in_file:
                            out_file.write(in_file.read())

rule collect_fpt:
        input:
                expand("experiments/{{fsg}}/fpt.timelimit={{timelimit}}.selector={{selector}}.lower-bound={{lower_bound}}.all=1.pre-mark={{pre_mark}}.search-strategy={{search_strategy}}/bio/{graph}.{multiplier}.{permutation}.solution.yaml",
                       graph=BIO_GRAPHS, multiplier=MULTIPLIER, permutation=PERMUTATION)
        output:
                "experiments/{fsg}/fpt.timelimit={timelimit}.selector={selector}.lower-bound={lower_bound}.all=1.pre-mark={pre_mark}.search-strategy={search_strategy}/all.solutions.yaml"
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
                iterations=10
        run:
                try:
                    subprocess.run(f"cmake-build-release/finder_benchmark --input {input.instance} --output {output} --finder {wildcards.finder} --iterations {params.iterations}".split(" "), timeout=2 * params.iterations + 2)
                except subprocess.TimeoutExpired:
                    pass

rule collect_finder_experiment:
        input:
                expand("experiments/finder-benchmark.finder={{finder}}/bio/{graph}.{permutation}.benchmark.yaml", graph=BIO_GRAPHS, permutation=PERMUTATION)
        output:
                "experiments/finder-benchmark.finder={finder}/all.benchmark.yaml"
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
