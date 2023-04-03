    cmake_minimum_required(VERSION 3.24)
    project(ILISI2_atelier_Shortest_path_A_start)

    set(CMAKE_CXX_STANDARD 17)

    find_package(Boost REQUIRED)

    find_package(Eigen3 3.3 REQUIRED NO_MODULE)



    add_executable(ILISI2_atelier_Shortest_path_A_start main.cpp)

    target_link_libraries(ILISI2_atelier_Shortest_path_A_start Eigen3::Eigen)

