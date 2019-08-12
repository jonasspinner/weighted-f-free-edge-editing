cmake_minimum_required(VERSION 3.13)
project(Weighted_F-free_Edge_Editing)

set(CMAKE_CXX_STANDARD 17)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wold-style-cast -Wcast-align -Wunused -Woverloaded-virtual -Wpedantic -Wsign-conversion -Wmisleading-indentation -Wduplicated-cond -Wduplicated-branches -Wlogical-op -Wnull-dereference -Wuseless-cast -Wdouble-promotion -Wformat=2")
# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wold-style-cast -Wcast-align -Wunused -Woverloaded-virtual -Wpedantic -Wconversion -Wsign-conversion -Wmisleading-indentation -Wduplicated-cond -Wduplicated-branches -Wlogical-op -Wnull-dereference -Wuseless-cast -Wdouble-promotion -Wformat=2")
# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic")

find_package(Boost REQUIRED)

if(Boost_FOUND)
    message("Boost Found")
else()
    error("Boost Not Found")
endif()


add_executable(graphedit src/main.cpp src/graph/Graph.h src/graph/VertexPairMap.h src/graph/VertexMap.h src/finder/NaiveP3.h src/graph/GraphIO.h src/finder/CenterC4P4.h src/lower_bound/ExplicitMWISSolver.h src/graph/Subgraph.h src/Editor.h src/interfaces/LowerBoundI.h src/interfaces/FinderI.h src/interfaces/ConsumerI.h src/lower_bound/NoLowerBound.h src/selector/FirstEditable.h src/selector/LeastWeight.h src/interfaces/SelectorI.h src/Configuration.h src/definitions.h src/Instance.h src/Statistics.h src/finder/CenterP3.h src/lower_bound/IteratedLocalSearch.h src/Edits.h src/Bound.h src/tests/FinderTests.h src/lower_bound/GreedyLowerBound.h src/tests/Tests.h src/finder/NaiveC4P4.h src/finder/Center.h src/consumer/SubgraphStats.h src/tests/GraphTests.h src/tests/SubgraphTests.h src/tests/EditorTests.h)

target_link_libraries(graphedit PRIVATE Boost::boost)