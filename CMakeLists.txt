cmake_minimum_required(VERSION 3.10)

project(MPI_Solver_Project)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

# Общая библиотека с реализацией функций
add_library(solver solver.cpp)

# Основная программа
add_executable(main main.cpp)
target_link_libraries(main solver ${MPI_LIBRARIES})

# Тесты
add_executable(tests tests.cpp)
target_link_libraries(tests solver ${MPI_LIBRARIES})

if(MPI_COMPILE_FLAGS)
    set_target_properties(main PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
    set_target_properties(tests PROPERTIES COMPILE_FLAGS "${MPI_COMPILE_FLAGS}")
endif()

if(MPI_LINK_FLAGS)
    set_target_properties(main PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
    set_target_properties(tests PROPERTIES LINK_FLAGS "${MPI_LINK_FLAGS}")
endif()
