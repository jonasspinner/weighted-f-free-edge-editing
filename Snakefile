from glob import glob
import subprocess
import yaml

GRAPHS = [path.split("/")[-1][:-6] for path in glob("data/bio/*.graph")]
# GRAPHS = GRAPHS[:10]
MULTIPLIERS = [100]
PERMUTATIONS = [0]
TIMELIMITS = [100]
FSG = ["P4C4"]
ILP_CONSTRAINTS = ["basic", "sparse", "single"]
ILP_SINGLE_CONSTRAINTS = [0]
ILP_NUM_THREADS = [1]


# FINDERS = ["CenterRecC5P5", "CenterRecC4P4", "CenterRecP3", "EndpointRecC5P5", "EndpointRecC4P4", "EndpointRecP3", "CenterC4P4", "CenterP3", "NaiveC4P4", "NaiveP3"]
FINDERS = ["CenterRecC4P4", "CenterRecP3", "EndpointRecC4P4", "EndpointRecP3", "CenterC4P4", "CenterP3"]


rule all:
        input:
                expand("experiments/{fsg}/ilp.timelimit={timelimit}.threads={threads}.constraints={constraints}/bio/{graph}.{multiplier}.{permutation}.solution.yaml", fsg=FSG, timelimit=TIMELIMITS, threads=ILP_NUM_THREADS, constraints=ILP_CONSTRAINTS, graph=GRAPHS, multiplier=MULTIPLIERS, permutation=PERMUTATIONS),
                # expand("experiments/{fsg}/fpt/bio/{graph}.{multiplier}.{permutation}.solution.yaml", fsg=FSG, timelimit=TIMELIMITS, threads=ILP_NUM_THREADS, sparse=ILP_SPARSE_CONSTRAINTS, graph=GRAPHS, multiplier=MULTIPLIERS, permutation=PERMUTATIONS),
                expand("experiments/finder-benchmark.finder={finder}/all.benchmark.yaml", finder=FINDERS)
                "data/bio/bio.metadata.yaml"

rule ilp:
        input:
                "data/{dataset}/{graph}.graph"
        output:
                "experiments/{fsg}/ilp.timelimit={timelimit}.threads={threads}.constraints={constraints}/{dataset}/{graph}.{multiplier}.{permutation}.solution.yaml"
        params:
                soft_timeout=lambda wildcards, output: int(1.1 * int(wildcards.timelimit))
        run:
                constraint_args = dict(basic=[], sparse=["--sparse-constraints", "1"], single=["--single-constraints", "1"])
                try:
                    subprocess.run(f"cmake-build-release/ilp --input {input}"
                                    "--multiplier {wildcards.multiplier}"
                                    "--permutation {wildcards.permutation}"
                                    "--output {output}"
                                    "--timelimit {wildcards.timelimit}"
                                    "--num-threads {wildcards.threads}".split(" ") + constraint_args[wildcards.constraints], timeout=params.soft_timeout)
                except subprocess.TimeoutExpired:
                    pass

rule copy_instance_solution:
        input:
                f"experiments/{{fsg}}/ilp.timelimit={TIMELIMITS[0]}.threads={ILP_NUM_THREADS[0]}.constraints={ILP_CONSTRAINTS[0]}/{{dataset}}/{{graph}}.{{multiplier}}.{PERMUTATIONS[0]}.solution.yaml"
        output:
                "experiments/solutions/{fsg}/{dataset}/{graph}.{multiplier}.solution.yaml"
        shell:
                "cp {input} {output}"


rule fpt_fixed_k_from_solution:
        input:
                instance = "data/{dataset}/{graph}.graph",
                solution = "experiments/solutions/{fsg}/{dataset}/{graph}.{multiplier}.solution.yaml"
        output:
                "experiments/{fsg}/fpt/{dataset}/{graph}.{multiplier}.{permutation}.solution.yaml"
        run:
                k = yaml.safe_load(open(input.solution))["solution_cost"]
                shell(f"cmake-build-release/fpt --input {input.instance} --multiplier {wildcards.multiplier} --permutation {wildcards.permutation} --k {k} --output {output}")


rule finder_experiment:
        input:
                instance = "data/{dataset}/{graph}.graph"
        output:
                "experiments/finder-benchmark.finder={finder}/{dataset}/{graph}.{permutation}.benchmark.yaml"
        run:
                try:
                    subprocess.run(f"cmake-build-release/finder_benchmark --input {input.instance} --output {output} --finder {wildcards.finder} --iterations 10".split(" "), timeout=10)
                except subprocess.TimeoutExpired:
                    pass

rule collect_finder_experiment:
        input:
                expand("experiments/finder-benchmark.finder={finder}/bio/{graph}.{permutation}.benchmark.yaml", fsg=FSG, finder=FINDERS, graph=GRAPHS, permutation=PERMUTATIONS)
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
