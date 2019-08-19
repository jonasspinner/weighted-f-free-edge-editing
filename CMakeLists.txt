cmake_minimum_required(VERSION 3.13)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/modules/")

project(Weighted_F-free_Edge_Editing)

set(CMAKE_CXX_STANDARD 17)

set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wold-style-cast -Wcast-align -Wunused -Woverloaded-virtual -Wpedantic -Wsign-conversion -Wmisleading-indentation -Wduplicated-cond -Wduplicated-branches -Wlogical-op -Wnull-dereference -Wuseless-cast -Wdouble-promotion -Wformat=2")
# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wshadow -Wnon-virtual-dtor -Wold-style-cast -Wcast-align -Wunused -Woverloaded-virtual -Wpedantic -Wconversion -Wsign-conversion -Wmisleading-indentation -Wduplicated-cond -Wduplicated-branches -Wlogical-op -Wnull-dereference -Wuseless-cast -Wdouble-promotion -Wformat=2")
# set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -Wshadow -Wnon-virtual-dtor -pedantic")

find_package(Boost COMPONENTS program_options REQUIRED)
include_directories(${Boost_INCLUDE_DIR} )

find_package(GUROBI REQUIRED)
include_directories(SYSTEM ${GUROBI_INCLUDE_DIRS})
set(LIBS ${LIBS} ${GUROBI_CXX_LIBRARY} ${GUROBI_LIBRARIES} ${GUROBI_CXX_LIBRARY})


#  GUROBI_FOUND - System has Gurobi
#  GUROBI_INCLUDE_DIRS - The Gurobi include directories
#  GUROBI_LIBRARIES - The libraries needed to use Gurobi
#message(${GUROBI_CXX_LIBRARY})
#message(${GUROBI_INCLUDE_DIRS})
#message(${GUROBI_LIBRARIES})
#message(${LIBS})

add_executable(graphedit src/main.cpp apps/run_tests.cpp)

# target_link_libraries(graphedit PRIVATE Boost::boost)
target_link_libraries(graphedit LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT} ${LIBS})
target_link_libraries(graphedit LINK_PUBLIC ${Boost_LIBRARIES} )


add_executable(fpt apps/run_fptsolver.cpp)

target_link_libraries(fpt LINK_PUBLIC ${Boost_LIBRARIES} )


add_executable(ilp apps/run_ilpsolver.cpp)

target_link_libraries(ilp LINK_PUBLIC ${CMAKE_THREAD_LIBS_INIT} ${LIBS})
target_link_libraries(ilp LINK_PUBLIC ${Boost_LIBRARIES} )


add_executable(tests apps/run_tests.cpp)

target_link_libraries(tests LINK_PUBLIC ${Boost_LIBRARIES} )
