cmake_minimum_required(VERSION 3.13)


project(Weighted_F-free_Edge_Editing)

include_directories(SYSTEM "${PROJECT_SOURCE_DIR}/extern/robin-hood-hashing/src/include")


list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake/modules/")
set(CMAKE_CXX_STANDARD 17)

# Compiler flags
add_compile_options(-Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wold-style-cast -Wcast-align -Wunused
        -Woverloaded-virtual -Wpedantic -Wsign-conversion -Wmisleading-indentation -Wduplicated-cond
        -Wduplicated-branches -Wlogical-op -Wnull-dereference -Wuseless-cast -Wdouble-promotion -Wformat=2 -Wno-comment)

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O3 -march=native -DNDEBUG -Wno-unused-parameter")


# TMP_DIR
if (NOT TMP_DIR)
    if (NOT DEFINED ENV{TMP_DIR})
        if (NOT DEFINED ENV{TMPDIR})
            set(TMP_DIR /tmp CACHE PATH "Path to tmp directory")
            message(WARNING "TMP_DIR is not defined.")
        else ()
            set(TMP_DIR $ENV{TMPDIR} CACHE PATH "Path to tmp directory")
        endif()
    else ()
        set(TMP_DIR $ENV{TMP_DIR} CACHE PATH "Path to tmp directory")
    endif ()
    message(STATUS "Using ${TMP_DIR} as TMP_DIR.")
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
    file(GLOB NPS_MWIS_SRC_FILES extern/nps_mwis/src/*.cpp)
    list(FILTER NPS_MWIS_SRC_FILES EXCLUDE REGEX ".*/main.cpp")
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
file(GLOB TEST_SRC_FILES src/tests/*.cpp)
file(GLOB LEGACY_FINDER_SRC_FILES src/legacy/finder/*.cpp)
file(GLOB SELECTOR_SRC_FILES src/selector/*.cpp)
file(GLOB SEARCH_STRATEGY_SRC_FILES src/search_strategy/*.cpp)

file(GLOB LOWER_BOUND_SRC_FILES src/lower_bound/*.cpp)
#list(APPEND LOWER_BOUND_SRC_FILES src/lower_bound/weighted_packing/WeightedPackingLocalSearch2.cpp)
list(FILTER LOWER_BOUND_SRC_FILES EXCLUDE REGEX "(.*/NPS_MWIS_Solver.cpp)|(.*/LSSWZ_MWIS_Solver.cpp)")

if (NPS_MWIS_FOUND)
    list(APPEND LOWER_BOUND_SRC_FILES src/lower_bound/NPS_MWIS_Solver.cpp)
endif (NPS_MWIS_FOUND)
if (LSSWZ_MWIS_FOUND)
    list(APPEND LOWER_BOUND_SRC_FILES src/lower_bound/LSSWZ_MWIS_Solver.cpp)
endif (LSSWZ_MWIS_FOUND)

set(SRC_FILES ${SELECTOR_SRC_FILES} ${LOWER_BOUND_SRC_FILES} ${SEARCH_STRATEGY_SRC_FILES} src/Configuration.cpp src/options.cpp ${NPS_MWIS_SRC_FILES})


add_executable(fpt apps/run_fptsolver.cpp ${SRC_FILES})
add_executable(tests apps/tests/run_tests.cpp ${SRC_FILES} ${TEST_SRC_FILES})
add_executable(legacy_finder_benchmark apps/legacy/run_legacy_finder_benchmark.cpp ${SRC_FILES} ${LEGACY_FINDER_SRC_FILES} )
add_executable(lower_bound_benchmark apps/benchmarks/run_lower_bound_benchmark.cpp ${SRC_FILES})
add_executable(fpt_experiment apps/run_fpt_experiment.cpp ${SRC_FILES})
add_executable(lower_bound_inspection apps/run_lower_bound_inspection.cpp ${SRC_FILES})
add_executable(instance_converter apps/run_instance_converter.cpp)
add_executable(instance_packing_problem apps/run_instance_packing_problem.cpp src/options.cpp)
add_executable(templated_vs_legacy_subgraph_test apps/run_templated_vs_legacy_subgraph_test.cpp src/Configuration.cpp src/options.cpp src/legacy/finder/CenterC4P4.cpp)


add_dependencies(fpt version)
add_dependencies(tests version)
add_dependencies(legacy_finder_benchmark version)
add_dependencies(lower_bound_benchmark version)
add_dependencies(fpt_experiment version)


if (GUROBI_FOUND)
    add_executable(ilp apps/run_ilpsolver.cpp src/Configuration.cpp src/options.cpp)
    add_dependencies(ilp version)
    target_link_libraries(ilp LINK_PUBLIC ${Boost_LIBRARIES} ${YAML_CPP_LIBRARIES})
endif (GUROBI_FOUND)

# Boost, YAML
target_link_libraries(fpt LINK_PUBLIC ${Boost_LIBRARIES} ${YAML_CPP_LIBRARIES})
target_link_libraries(tests LINK_PUBLIC ${Boost_LIBRARIES} ${YAML_CPP_LIBRARIES})
target_link_libraries(legacy_finder_benchmark LINK_PUBLIC ${Boost_LIBRARIES} ${YAML_CPP_LIBRARIES})
target_link_libraries(lower_bound_benchmark LINK_PUBLIC ${Boost_LIBRARIES} ${YAML_CPP_LIBRARIES})
target_link_libraries(fpt_experiment LINK_PUBLIC ${Boost_LIBRARIES} ${YAML_CPP_LIBRARIES})
target_link_libraries(lower_bound_inspection LINK_PUBLIC ${Boost_LIBRARIES} ${YAML_CPP_LIBRARIES})
target_link_libraries(instance_converter LINK_PUBLIC ${Boost_LIBRARIES} ${YAML_CPP_LIBRARIES})
target_link_libraries(instance_packing_problem LINK_PUBLIC ${Boost_LIBRARIES} ${YAML_CPP_LIBRARIES})
target_link_libraries(templated_vs_legacy_subgraph_test LINK_PUBLIC ${Boost_LIBRARIES} ${YAML_CPP_LIBRARIES})

# Gurobi
# TODO: Only link if necessary.
if (GUROBI_FOUND)
    target_link_libraries(fpt LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARIES} ${GUROBI_CXX_LIBRARY})
    target_link_libraries(ilp LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARIES} ${GUROBI_CXX_LIBRARY})
    target_link_libraries(tests LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARIES} ${GUROBI_CXX_LIBRARY})
    target_link_libraries(legacy_finder_benchmark LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARIES} ${GUROBI_CXX_LIBRARY})
    target_link_libraries(lower_bound_benchmark LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARIES} ${GUROBI_CXX_LIBRARY})
    target_link_libraries(fpt_experiment LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARIES} ${GUROBI_CXX_LIBRARY})
    target_link_libraries(lower_bound_inspection LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARIES} ${GUROBI_CXX_LIBRARY})
endif (GUROBI_FOUND)
