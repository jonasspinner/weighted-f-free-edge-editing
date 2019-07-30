cmake_minimum_required(VERSION 3.14)
project(Weighted_F-free_Edge_Editing)

set(CMAKE_CXX_STANDARD 17)

add_executable(graphedit src/main.cpp src/graph/Graph.h src/graph/VertexPairMap.h src/graph/VertexMap.h src/finder/NaiveP3.h src/graph/GraphIO.h src/finder/CenterC4P4.h src/lower_bound/ExplicitMWISSolver.h src/graph/Subgraph.h src/Editor.h src/interfaces/LowerBoundI.h src/interfaces/FinderI.h src/interfaces/ConsumerI.h src/lower_bound/NoLowerBound.h src/selector/FirstEditable.h src/selector/LeastWeight.h src/interfaces/SelectorI.h src/Configuration.h src/Definitions.h src/Instance.h src/Statistics.h src/finder/CenterP3.h src/lower_bound/IteratedLocalSearch.h src/Edits.h src/Bound.h src/tests/FinderTests.h src/lower_bound/GreedyLowerBound.h src/tests/Tests.h src/finder/NaiveC4P4.h)
