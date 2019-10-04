from glob import glob
GRAPHS = [path.split("/")[-1][:-6] for path in glob("data/bio/*.graph")]
GRAPHS = GRAPHS[:100]
MULTIPLIERS = [10, 100]
PERMUTATIONS = [0, 1]
TIMELIMITS = [100]
FSG = ["P4C4"]
SPARSE_CONSTRAINTS = [0]
EXTENDED_CONSTRAINTS = [0]
NUM_THREADS = [1]


rule all:
        input:
                expand("experiments/{fsg}/ilp.timelimit={timelimit}.threads={threads}.sparse={sparse}.extended={extended}/bio/{graph}.{multiplier}.{permutation}.solution.yaml", fsg=FSG, timelimit=TIMELIMITS, threads=NUM_THREADS, sparse=SPARSE_CONSTRAINTS, extended=EXTENDED_CONSTRAINTS, graph=GRAPHS, multiplier=MULTIPLIERS, permutation=PERMUTATIONS),
                expand("experiments/{fsg}/fpt/bio/{graph}.{multiplier}.{permutation}.solution.yaml", fsg=FSG, timelimit=TIMELIMITS, threads=NUM_THREADS, sparse=SPARSE_CONSTRAINTS, extended=EXTENDED_CONSTRAINTS, graph=GRAPHS, multiplier=MULTIPLIERS, permutation=PERMUTATIONS),
                "data/bio/bio.metadata.yaml"

rule ilp:
        input:
                "data/{dataset}/{graph}.graph"
        output:
                "experiments/{fsg}/ilp.timelimit={timelimit}.threads={threads}.sparse={sparse}.extended={extended}/{dataset}/{graph}.{multiplier}.{permutation}.solution.yaml"
        shell:
                "build-release/ilp --input {input} --multiplier {wildcards.multiplier} --permutation {wildcards.permutation} --output {output} --timelimit {wildcards.timelimit} --num-threads {wildcards.threads} --sparse-constraints {wildcards.sparse} --extended-constraints {wildcards.extended}"

rule fpt:
        input:
                "data/{dataset}/{graph}.graph"
        output:
                "experiments/{fsg}/fpt/{dataset}/{graph}.{multiplier}.{permutation}.solution.yaml"
        shell:
                "build-release/fpt --input {input} --multiplier {wildcards.multiplier} --permutation {wildcards.permutation} --output {output}"


rule metadata:
        input:
                "data/{dataset}"
        output:
                "data/{dataset}/{dataset}.metadata.yaml"
        shell:
                "python3 scripts/compute_metadata.py {input}"