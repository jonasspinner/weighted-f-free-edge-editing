cmake_minimum_required(VERSION 3.13)


project(Weighted_F-free_Edge_Editing)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")
set(CMAKE_CXX_STANDARD 17)

# Compiler flags
add_compile_options(-Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wold-style-cast -Wcast-align -Wunused -Woverloaded-virtual -Wpedantic -Wsign-conversion -Wmisleading-indentation -Wduplicated-cond -Wduplicated-branches -Wlogical-op -Wnull-dereference -Wuseless-cast -Wdouble-promotion -Wformat=2)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 -march=native -Wno-comment")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG -Wno-unused-parameter")



# Include directories

# Boost
find_package(Boost COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

# YAML
find_package(YAML-CPP REQUIRED)
include_directories(${YAML_CPP_INCLUDE_DIR})

# Gurobi
find_package(GUROBI)
if(GUROBI_FOUND)
    include_directories(SYSTEM ${GUROBI_INCLUDE_DIRS})
    add_definitions("-DGUROBI_FOUND")
endif(GUROBI_FOUND)

# NPS MWIS
# TODO: Refactor as optional dependency
# TODO: Refactor as CMAKE subfolder
set(NPS_MWIS_FOUND TRUE)
if(NPS_MWIS_FOUND)
    add_definitions(-DNPS_MWIS_FOUND)
    set(NPS_MWIS_SRC_FILES extern/nps_mwis/src/ArgPack.cpp extern/nps_mwis/src/bossa_timer.cpp extern/nps_mwis/src/Graph.cpp extern/nps_mwis/src/Solution.cpp)
    set_source_files_properties(${NPS_MWIS_SRC_FILES} PROPERTIES COMPILE_FLAGS "-w")
endif(NPS_MWIS_FOUND)

# LSSWZ MWIS
set(LSSWZ_MWIS_FOUND TRUE)
if(LSSWZ_MWIS_FOUND)
    add_definitions(-DLSSWZ_MWIS_FOUND)
    add_subdirectory(extern/lsswz_mwis/code)
endif(LSSWZ_MWIS_FOUND)


# Source files

set(TEST_SRC_FILES src/tests/EditorTests.cpp src/tests/GraphTests.cpp src/tests/test_utils.cpp)
set(FINDER_SRC_FILES src/finder/finder_utils.cpp src/finder/CenterC4P4.cpp src/finder/Center.cpp src/finder/SplitCluster.cpp src/finder/SplitGraph.cpp src/finder/CenterP3.cpp src/finder/NaiveP3.cpp src/finder/Endpoint.cpp src/finder/Naive.cpp src/finder/OuterP3.cpp)
set(SELECTOR_SRC_FILES src/selector/selector_utils.cpp)

set(LOWER_BOUND_SRC_FILES src/lower_bound/lower_bound_utils.cpp src/lower_bound/TrivialLowerBound.cpp src/lower_bound/SortedGreedyLowerBound.cpp src/lower_bound/GreedyLowerBound.cpp src/lower_bound/LocalSearchLowerBound.cpp)
if(NPS_MWIS_FOUND)
    set(LOWER_BOUND_SRC_FILES ${LOWER_BOUND_SRC_FILES} src/lower_bound/NPS_MWIS_SolverLowerBound.cpp)
endif(NPS_MWIS_FOUND)
if(LSSWZ_MWIS_FOUND)
    set(LOWER_BOUND_SRC_FILES ${LOWER_BOUND_SRC_FILES} src/lower_bound/LSSWZ_MWIS_SolverLowerBound.cpp)
endif(LSSWZ_MWIS_FOUND)

set(SRC_FILES ${FINDER_SRC_FILES} ${SELECTOR_SRC_FILES} ${LOWER_BOUND_SRC_FILES} src/Configuration.cpp src/graph/Subgraph.cpp src/options.cpp ${NPS_MWIS_SRC_FILES})


add_executable(graphedit src/main.cpp ${SRC_FILES})
add_executable(fpt apps/run_fptsolver.cpp ${SRC_FILES})
add_executable(tests apps/run_tests.cpp ${SRC_FILES} ${TEST_SRC_FILES})
add_executable(finder_benchmark apps/run_finder_benchmark.cpp ${SRC_FILES})
add_executable(lower_bound_benchmark apps/run_lower_bound_benchmark.cpp ${SRC_FILES})
add_executable(fpt_experiment apps/run_fpt_experiment.cpp ${SRC_FILES})

if(GUROBI_FOUND)
    add_executable(ilp apps/run_ilpsolver.cpp ${SRC_FILES})
    target_link_libraries(ilp LINK_PUBLIC ${Boost_LIBRARIES} ${YAML_CPP_LIBRARIES})
endif(GUROBI_FOUND)

# Boost, YAML
target_link_libraries(graphedit LINK_PUBLIC ${Boost_LIBRARIES} ${YAML_CPP_LIBRARIES})
target_link_libraries(fpt LINK_PUBLIC ${Boost_LIBRARIES} ${YAML_CPP_LIBRARIES})
target_link_libraries(tests LINK_PUBLIC ${Boost_LIBRARIES} ${YAML_CPP_LIBRARIES})
target_link_libraries(finder_benchmark LINK_PUBLIC ${Boost_LIBRARIES} ${YAML_CPP_LIBRARIES})
target_link_libraries(lower_bound_benchmark LINK_PUBLIC ${Boost_LIBRARIES} ${YAML_CPP_LIBRARIES})
target_link_libraries(fpt_experiment LINK_PUBLIC ${Boost_LIBRARIES} ${YAML_CPP_LIBRARIES})

# Gurobi
if(GUROBI_FOUND)
    target_link_libraries(graphedit LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARIES} ${GUROBI_CXX_LIBRARY})
    target_link_libraries(fpt LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARIES} ${GUROBI_CXX_LIBRARY})
    target_link_libraries(ilp LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARIES} ${GUROBI_CXX_LIBRARY})
    target_link_libraries(tests LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARIES} ${GUROBI_CXX_LIBRARY})
    target_link_libraries(finder_benchmark LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARIES} ${GUROBI_CXX_LIBRARY})
    target_link_libraries(lower_bound_benchmark LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARIES} ${GUROBI_CXX_LIBRARY})
    target_link_libraries(fpt_experiment LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARIES} ${GUROBI_CXX_LIBRARY})
endif(GUROBI_FOUND)



# Add version command
add_custom_command(TARGET graphedit PRE_BUILD COMMAND ${CMAKE_COMMAND} -P "${CMAKE_CURRENT_LIST_DIR}/cmake/FillVersionH.cmake")
add_custom_command(TARGET fpt PRE_BUILD COMMAND ${CMAKE_COMMAND} -P "${CMAKE_CURRENT_LIST_DIR}/cmake/FillVersionH.cmake")
add_custom_command(TARGET ilp PRE_BUILD COMMAND ${CMAKE_COMMAND} -P "${CMAKE_CURRENT_LIST_DIR}/cmake/FillVersionH.cmake")
add_custom_command(TARGET tests PRE_BUILD COMMAND ${CMAKE_COMMAND} -P "${CMAKE_CURRENT_LIST_DIR}/cmake/FillVersionH.cmake")
add_custom_command(TARGET finder_benchmark PRE_BUILD COMMAND ${CMAKE_COMMAND} -P "${CMAKE_CURRENT_LIST_DIR}/cmake/FillVersionH.cmake")
add_custom_command(TARGET lower_bound_benchmark PRE_BUILD COMMAND ${CMAKE_COMMAND} -P "${CMAKE_CURRENT_LIST_DIR}/cmake/FillVersionH.cmake")
add_custom_command(TARGET fpt_experiment PRE_BUILD COMMAND ${CMAKE_COMMAND} -P "${CMAKE_CURRENT_LIST_DIR}/cmake/FillVersionH.cmake")
