cmake_minimum_required(VERSION 3.13)


project(Weighted_F-free_Edge_Editing)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")
set(CMAKE_CXX_STANDARD 17)

# Compiler flags
add_compile_options(-Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wold-style-cast -Wcast-align -Wunused
        -Woverloaded-virtual -Wpedantic -Wsign-conversion -Wmisleading-indentation -Wduplicated-cond
        -Wduplicated-branches -Wlogical-op -Wnull-dereference -Wuseless-cast -Wdouble-promotion -Wformat=2)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -Wno-comment")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG -Wno-unused-parameter")


# TMP_DIR
if (NOT TMP_DIR)
    if (NOT DEFINED ENV{TMP_DIR})
        set(TMP_DIR /tmp CACHE PATH "Path to tmp directory")
        message(STATUS "TMP_DIR is not defined. Using ${TMP_DIR} instead.")
    else ()
        set(TMP_DIR $ENV{TMP_DIR} CACHE PATH "Path to tmp directory")
    endif ()
endif ()

# Generate version.h
add_custom_target(
        version
        ${CMAKE_COMMAND}
        -D TMP_DIR=${TMP_DIR}
        -D SRC=${CMAKE_CURRENT_LIST_DIR}/src/version.h.in -D DST=${CMAKE_CURRENT_LIST_DIR}/src/version.h
        -P "${CMAKE_CURRENT_LIST_DIR}/cmake/FillVersionH.cmake"
)


# Include directories

# Boost
find_package(Boost COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIR})


# YAML
find_package(YAML-CPP REQUIRED)
include_directories(${YAML_CPP_INCLUDE_DIR})

# Gurobi
find_package(GUROBI)
if (GUROBI_FOUND)
    include_directories(SYSTEM ${GUROBI_INCLUDE_DIRS})
    add_definitions("-DGUROBI_FOUND")
endif (GUROBI_FOUND)


# NPS MWIS
find_package(NPS_MWIS)
if (NPS_MWIS_FOUND)
    add_definitions(-DNPS_MWIS_FOUND)
    set(NPS_MWIS_SRC_FILES extern/nps_mwis/src/ArgPack.cpp extern/nps_mwis/src/bossa_timer.cpp extern/nps_mwis/src/Graph.cpp extern/nps_mwis/src/Solution.cpp)
    set_source_files_properties(${NPS_MWIS_SRC_FILES} PROPERTIES COMPILE_FLAGS "-w")
else (NPS_MWIS_FOUND)
    message(STATUS "Not found NPS_MWIS: Download with extern/nps_mwis/get.sh")
endif (NPS_MWIS_FOUND)

# LSSWZ MWIS
find_package(LSSWZ_MWIS)
if (LSSWZ_MWIS_FOUND)
    add_definitions(-DLSSWZ_MWIS_SCRIPT=\"${LSSWZ_MWIS_SCRIPT}\")
    add_definitions(-DLSSWZ_MWIS_FOUND)
    # add_subdirectory(extern/lsswz_mwis/code)
endif ()


# Source files

set(TEST_SRC_FILES src/tests/EditorTests.cpp src/tests/GraphTests.cpp src/tests/test_utils.cpp)
set(FINDER_SRC_FILES src/finder/finder_utils.cpp src/finder/CenterC4P4.cpp src/finder/Center.cpp src/finder/SplitCluster.cpp src/finder/SplitGraph.cpp src/finder/CenterP3.cpp src/finder/NaiveP3.cpp src/finder/Endpoint.cpp src/finder/Naive.cpp src/finder/OuterP3.cpp)
set(SELECTOR_SRC_FILES src/selector/selector_utils.cpp)
set(SEARCH_STRATEGY_SRC_FILES src/search_strategy/search_strategy_utils.cpp)

set(LOWER_BOUND_SRC_FILES src/lower_bound/lower_bound_utils.cpp src/lower_bound/Trivial.cpp src/lower_bound/SortedGreedy.cpp src/lower_bound/Greedy.cpp src/lower_bound/LocalSearch.cpp)
if (NPS_MWIS_FOUND)
    set(LOWER_BOUND_SRC_FILES ${LOWER_BOUND_SRC_FILES} src/lower_bound/NPS_MWIS_Solver.cpp)
endif (NPS_MWIS_FOUND)
if (LSSWZ_MWIS_FOUND)
    set(LOWER_BOUND_SRC_FILES ${LOWER_BOUND_SRC_FILES} src/lower_bound/LSSWZ_MWIS_Solver.cpp)
endif (LSSWZ_MWIS_FOUND)

set(SRC_FILES ${FINDER_SRC_FILES} ${SELECTOR_SRC_FILES} ${LOWER_BOUND_SRC_FILES} ${SEARCH_STRATEGY_SRC_FILES} src/Configuration.cpp src/graph/Subgraph.cpp src/options.cpp ${NPS_MWIS_SRC_FILES})


add_executable(fpt apps/run_fptsolver.cpp ${SRC_FILES})
add_executable(tests apps/run_tests.cpp ${SRC_FILES} ${TEST_SRC_FILES})
add_executable(finder_benchmark apps/run_finder_benchmark.cpp ${SRC_FILES})
add_executable(lower_bound_benchmark apps/run_lower_bound_benchmark.cpp ${SRC_FILES})
add_executable(fpt_experiment apps/run_fpt_experiment.cpp ${SRC_FILES})


add_dependencies(fpt version)
add_dependencies(tests version)
add_dependencies(finder_benchmark version)
add_dependencies(lower_bound_benchmark version)
add_dependencies(fpt_experiment version)


if (GUROBI_FOUND)
    add_executable(ilp apps/run_ilpsolver.cpp ${SRC_FILES})
    add_dependencies(ilp version)
    target_link_libraries(ilp LINK_PUBLIC ${Boost_LIBRARIES} ${YAML_CPP_LIBRARIES})
endif (GUROBI_FOUND)

# Boost, YAML
target_link_libraries(fpt LINK_PUBLIC ${Boost_LIBRARIES} ${YAML_CPP_LIBRARIES})
target_link_libraries(tests LINK_PUBLIC ${Boost_LIBRARIES} ${YAML_CPP_LIBRARIES})
target_link_libraries(finder_benchmark LINK_PUBLIC ${Boost_LIBRARIES} ${YAML_CPP_LIBRARIES})
target_link_libraries(lower_bound_benchmark LINK_PUBLIC ${Boost_LIBRARIES} ${YAML_CPP_LIBRARIES})
target_link_libraries(fpt_experiment LINK_PUBLIC ${Boost_LIBRARIES} ${YAML_CPP_LIBRARIES})

# Gurobi
# TODO: Only link if necessary.
if (GUROBI_FOUND)
    target_link_libraries(fpt LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARIES} ${GUROBI_CXX_LIBRARY})
    target_link_libraries(ilp LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARIES} ${GUROBI_CXX_LIBRARY})
    target_link_libraries(tests LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARIES} ${GUROBI_CXX_LIBRARY})
    target_link_libraries(finder_benchmark LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARIES} ${GUROBI_CXX_LIBRARY})
    target_link_libraries(lower_bound_benchmark LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARIES} ${GUROBI_CXX_LIBRARY})
    target_link_libraries(fpt_experiment LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARIES} ${GUROBI_CXX_LIBRARY})
endif (GUROBI_FOUND)
