cmake_minimum_required(VERSION 3.13)
set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")


project(Weighted_F-free_Edge_Editing)

set(CMAKE_CXX_STANDARD 17)


set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wold-style-cast -Wcast-align -Wunused -Woverloaded-virtual -Wpedantic -Wsign-conversion -Wmisleading-indentation -Wduplicated-cond -Wduplicated-branches -Wlogical-op -Wnull-dereference -Wuseless-cast -Wdouble-promotion -Wformat=2")
# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wold-style-cast -Wcast-align -Wunused -Woverloaded-virtual -Wpedantic -Wconversion -Wsign-conversion -Wmisleading-indentation -Wduplicated-cond -Wduplicated-branches -Wlogical-op -Wnull-dereference -Wuseless-cast -Wdouble-promotion -Wformat=2")


execute_process(
        COMMAND git rev-parse --abbrev-ref HEAD
        WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
        OUTPUT_VARIABLE GIT_BRANCH
        OUTPUT_STRIP_TRAILING_WHITESPACE
)
#add_definitions("-DGIT_BRANCH=\"${GIT_BRANCH}\"")

execute_process(
    COMMAND git log -1 --format=%H
    WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}
    OUTPUT_VARIABLE GIT_COMMIT_HASH
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
# add_definitions("-DGIT_COMMIT_HASH=\"${GIT_COMMIT_HASH}\"")
configure_file(
        ${CMAKE_SOURCE_DIR}/src/version.h.in
        ${CMAKE_BINARY_DIR}/src/version.h
)



find_package(Boost COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIR} )

find_package(GUROBI)
if(GUROBI_FOUND)
    include_directories(SYSTEM ${GUROBI_INCLUDE_DIRS})
    # set(LIBS ${LIBS} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARIES} ${GUROBI_CXX_LIBRARY})
    add_definitions("-DGUROBI_FOUND")
endif(GUROBI_FOUND)


find_package(YAML-CPP REQUIRED)
include_directories(${YAML_CPP_INCLUDE_DIR})


add_executable(graphedit src/main.cpp)
target_link_libraries(graphedit LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARIES} ${GUROBI_CXX_LIBRARY})
target_link_libraries(graphedit LINK_PUBLIC ${Boost_LIBRARIES} ${YAML_CPP_LIBRARIES})


add_executable(fpt apps/run_fptsolver.cpp)
target_link_libraries(fpt LINK_PUBLIC ${Boost_LIBRARIES} ${YAML_CPP_LIBRARIES})


if(GUROBI_FOUND)
    add_executable(ilp apps/run_ilpsolver.cpp)
    target_link_libraries(ilp LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARIES} ${GUROBI_CXX_LIBRARY})
    target_link_libraries(ilp LINK_PUBLIC ${Boost_LIBRARIES} ${YAML_CPP_LIBRARIES})
endif(GUROBI_FOUND)


add_executable(tests apps/run_tests.cpp)
target_link_libraries(tests LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARIES} ${GUROBI_CXX_LIBRARY})
target_link_libraries(tests LINK_PUBLIC ${Boost_LIBRARIES} ${YAML_CPP_LIBRARIES})

add_executable(finder_benchmark apps/run_finder_benchmark.cpp)
target_link_libraries(finder_benchmark LINK_PUBLIC ${Boost_LIBRARIES} ${YAML_CPP_LIBRARIES})

add_executable(lower_bound_benchmark apps/run_lower_bound_benchmark.cpp)
target_link_libraries(lower_bound_benchmark LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARIES} ${GUROBI_CXX_LIBRARY})
target_link_libraries(lower_bound_benchmark LINK_PUBLIC ${Boost_LIBRARIES} ${YAML_CPP_LIBRARIES})
