import urllib.request
import yaml
from tempfile import NamedTemporaryFile, TemporaryDirectory
from pathlib import Path
import zipfile
import re
from shutil import rmtree


def output(text: str):
    print(text)


def transform(input_path: Path, output_path: Path):
    with input_path.open('r') as input_file:
        lines = input_file.read().split('\n')
    n = int(lines[0])
    m = n * (n - 1) // 2
    fmt = 1

    M = [[0.0] + [float(value) for value in line.split('\t')] for line in lines[n + 1:-2]] + [[0.0]]
    S = [[M[min([i, j])][max([i, j]) - min([i, j])] for j in range(n)] for i in range(n)]

    output(f"\twriting {output_path.name}")

    with output_path.open('w') as output_file:
        output_file.write(f"{n} {m} {fmt}\n")
        for i in range(n):
            output_file.write(" ".join([f"{j+1} {S[i][j]}" for j in range(n) if i != j]))
            output_file.write("\n")


def download_bio_dataset():
    folder_path = Path("bio")
    if folder_path.exists():
        output(f"{folder_path.name} already exists")
        pass

    output("loading config")

    with Path("bio.yaml").open('r') as file:
        config = yaml.safe_load(file)

    output("downloading data")

    response = urllib.request.urlopen(config['download_url'])
    CHUNK = 16 * 1024
    with NamedTemporaryFile() as file:
        for chunk in iter(lambda: response.read(CHUNK), b''):
            file.write(chunk)
        with zipfile.ZipFile(file.name, 'r') as zip_ref:
            zip_ref.extractall(folder_path)

    output("transforming files")

    for file_path in folder_path.glob("**/*.cm"):
        match = re.match(r"cost_matrix_component_nr_(\d+)_size_(\d+)_cutoff_10.0.cm", file_path.name)
        number, size = match.group(1, 2)
        if int(size) >= 1000:
            continue

        output_path = folder_path / f"bio-nr-{number}-size-{size}.metis"
        try:
            transform(file_path, output_path)
        except MemoryError:
            output(f"\t{output_path.name} encountered memory error")

    output("deleting original data")

    rmtree(folder_path / "biological")


def main():
    download_bio_dataset()


if __name__ == '__main__':
    main()
