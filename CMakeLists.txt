cmake_minimum_required(VERSION 3.13)


project(Weighted_F-free_Edge_Editing)


set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")
set(CMAKE_CXX_STANDARD 17)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wold-style-cast -Wcast-align -Wunused -Woverloaded-virtual -Wpedantic -Wsign-conversion -Wmisleading-indentation -Wduplicated-cond -Wduplicated-branches -Wlogical-op -Wnull-dereference -Wuseless-cast -Wdouble-promotion -Wformat=2")

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
#set(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DNDEBUG")





find_package(Boost COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIR})

find_package(GUROBI)
if(GUROBI_FOUND)
    include_directories(SYSTEM ${GUROBI_INCLUDE_DIRS})
    # set(LIBS ${LIBS} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARIES} ${GUROBI_CXX_LIBRARY})
    add_definitions("-DGUROBI_FOUND")
endif(GUROBI_FOUND)


find_package(YAML-CPP REQUIRED)
include_directories(${YAML_CPP_INCLUDE_DIR})


set(TEST_SRC_FILES src/tests/EditorTests.cpp src/tests/GraphTests.cpp src/tests/test_utils.cpp)
set(FINDER_SRC_FILES src/finder/finder_utils.cpp src/finder/CenterC4P4.cpp src/finder/Center.cpp src/finder/SplitCluster.cpp src/finder/SplitGraph.cpp src/finder/CenterP3.cpp src/finder/NaiveP3.cpp)
set(SELECTOR_SRC_FILES src/selector/selector_utils.cpp)
set(LOWER_BOUND_SRC_FILES src/lower_bound/lower_bound_utils.cpp src/lower_bound/TrivialLowerBound.cpp src/lower_bound/GlobalGreedyLowerBound.cpp src/lower_bound/LocalGreedyLowerBound.cpp)
set(SRC_FILES ${FINDER_SRC_FILES} ${SELECTOR_SRC_FILES} ${LOWER_BOUND_SRC_FILES} src/Configuration.cpp src/graph/Subgraph.cpp)


add_executable(graphedit src/main.cpp ${SRC_FILES})
add_executable(fpt apps/run_fptsolver.cpp ${SRC_FILES})
add_executable(tests apps/run_tests.cpp ${SRC_FILES} ${TEST_SRC_FILES})
add_executable(finder_benchmark apps/run_finder_benchmark.cpp ${SRC_FILES})
add_executable(lower_bound_benchmark apps/run_lower_bound_benchmark.cpp ${SRC_FILES})

if(GUROBI_FOUND)
    add_executable(ilp apps/run_ilpsolver.cpp ${SRC_FILES})
    target_link_libraries(ilp LINK_PUBLIC ${Boost_LIBRARIES} ${YAML_CPP_LIBRARIES})
endif(GUROBI_FOUND)

target_link_libraries(graphedit LINK_PUBLIC ${Boost_LIBRARIES} ${YAML_CPP_LIBRARIES})
target_link_libraries(fpt LINK_PUBLIC ${Boost_LIBRARIES} ${YAML_CPP_LIBRARIES})
target_link_libraries(tests LINK_PUBLIC ${Boost_LIBRARIES} ${YAML_CPP_LIBRARIES})
target_link_libraries(finder_benchmark LINK_PUBLIC ${Boost_LIBRARIES} ${YAML_CPP_LIBRARIES})
target_link_libraries(lower_bound_benchmark LINK_PUBLIC ${Boost_LIBRARIES} ${YAML_CPP_LIBRARIES})

if(GUROBI_FOUND)
    target_link_libraries(graphedit LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARIES} ${GUROBI_CXX_LIBRARY})
    target_link_libraries(fpt LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARIES} ${GUROBI_CXX_LIBRARY})
    target_link_libraries(ilp LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARIES} ${GUROBI_CXX_LIBRARY})
    target_link_libraries(tests LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARIES} ${GUROBI_CXX_LIBRARY})
    target_link_libraries(lower_bound_benchmark LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARIES} ${GUROBI_CXX_LIBRARY})
endif(GUROBI_FOUND)


add_custom_command(TARGET graphedit PRE_BUILD COMMAND ${CMAKE_COMMAND} -P "${CMAKE_CURRENT_LIST_DIR}/FillVersionH.cmake")
add_custom_command(TARGET fpt PRE_BUILD COMMAND ${CMAKE_COMMAND} -P "${CMAKE_CURRENT_LIST_DIR}/FillVersionH.cmake")
add_custom_command(TARGET ilp PRE_BUILD COMMAND ${CMAKE_COMMAND} -P "${CMAKE_CURRENT_LIST_DIR}/FillVersionH.cmake")
add_custom_command(TARGET tests PRE_BUILD COMMAND ${CMAKE_COMMAND} -P "${CMAKE_CURRENT_LIST_DIR}/FillVersionH.cmake")
add_custom_command(TARGET finder_benchmark PRE_BUILD COMMAND ${CMAKE_COMMAND} -P "${CMAKE_CURRENT_LIST_DIR}/FillVersionH.cmake")
add_custom_command(TARGET lower_bound_benchmark PRE_BUILD COMMAND ${CMAKE_COMMAND} -P "${CMAKE_CURRENT_LIST_DIR}/FillVersionH.cmake")