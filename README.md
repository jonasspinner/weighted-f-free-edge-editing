# Weighted F-free Edge Editing

This is the source code for a bachelor thesis regarding the Weighted F-free Edge Editing Problem. We implement a FPT and an ILP algorithm to solve the problem exactly. The FPT algorithm and ideas are based on the following paper for (unweighted) F-free Edge Editing:

```
Lars Gottesbüren, Michael Hamann, Philipp Schoch, Ben Strasser, Dorothea Wagner, and Sven Zühlsdorf.
Engineering Exact F-free  Edge Editing. Unpublished, 2019.
```

The license of their work can be found under `external/GHSSWZ_LICENSE`.

A detailed introduction to the problem and a description of the algorithms are given in the thesis.

## How to run

```bash
mkdir build && cd build
cmake .. -DCMAKE_BUILD_TYPE=Release
make
```


## Dependencies

### [BOOST](https://www.boost.org/)

### [GUROBI](https://www.gurobi.com/) *(optional)*

### [YAML CPP](https://github.com/jbeder/yaml-cpp/)

```bash
sudo apt install libyaml-cpp-dev
```


## FPT Algorithm

```
Title:	FPT Algorithm for Weighted F-free Edge Editing
Input:
	A graph G = (V, E),
	a set of forbidden subgraphs F,
	a cost function c: (V choose 2) -> R_{>= 0},
	a maximum editing cost k in R,
	the set of marked vertex pairs M subset of (V choose 2) and
	the set of currently edited vertex pairs L.
	
	// F and c are fixed throughout the algorithm.
Output:
	Solutions to the F-free Editing Problem} with total editing costs of at most k.

if k < LowerBound(G, k, M):
	return

S := FindSubgraph(G, M)
if S is none:
	output L
	return

m := {}						// locally marked vertex pairs
for all e in (S choose 2):
	if e not in M:
		M += {e}; m += {e}; L += {e}
		G.toggle(e)
		Edit(G, k - c(e), M, L) 	// recursive call
		G.toggle(e)
		L -= {e}

M -= m						// unmark vertex pairs

```

## How to cite

```
Jonas Spinner.
Weighted F-free Edge Editing. Bachelor's thesis, Karlsruhe, 2019.
```


```bibtex
@thesis{Spinner2019,
	title = {{Weighted} {F}-free {Edge} {Editing}},
	language = {en},
	school = {Karlsruhe Institute of Technology},
	author = {Spinner, Jonas},
	year = {2019},
	type={Bachelor's thesis}
}
```
