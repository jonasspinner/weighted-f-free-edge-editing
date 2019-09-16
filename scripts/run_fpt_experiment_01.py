import subprocess
import yaml


def main(solutions_path):
    results = []
    try:
        with open(solutions_path, "r") as file:
            for doc in yaml.safe_load_all(file):
                path = doc["instance"]
                k = doc["solution_cost"]

                if k < 0 or doc["time"][0] > 1000000:
                    continue

                print(doc)

                args = ["../cmake-build-release/fpt_experiment"]
                options = {"input": f"../{path}", "k": k, "warmup-steps": 10, "steps": min([k, 10])}

                for k, v in options.items():
                    args += [f"--{k}", str(v)]

                result = subprocess.run(args, capture_output=True)
                result_doc = yaml.safe_load(result.stdout)
                results.append(result_doc)
                print(result_doc)
    except KeyboardInterrupt:
        pass
    finally:
        with open("output.yaml", "w") as file:
            yaml.safe_dump_all(results, file)


if __name__ == "__main__":
    main("../output/bio_solutions.yaml")
