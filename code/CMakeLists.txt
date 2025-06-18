cmake_minimum_required(VERSION 3.16)
project(Masterthesis)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Optional: add include directory if your headers are separated
# include_directories(include)

add_executable(main
        src/main.cpp
        src/Graph.cpp
        src/Quadtree.cpp
        src/Dataloader.cpp
        src/Dataloader.hpp
)

add_executable(WSPDTest
        test/WSPDTest.cpp
        src/Graph.cpp
        src/Quadtree.cpp
        src/Dataloader.cpp
        src/Dataloader.hpp
        src/ThetaSpanner.cpp
        src/HubLabels.hpp
        src/HubLabels.cpp
)

add_executable(GraphTest
        src/Graph.cpp
        src/Dataloader.cpp
        src/Dataloader.hpp
        test/GraphTest.cpp
)

add_executable(Theta
        src/Graph.cpp
        src/Dataloader.cpp
        src/Dataloader.hpp
        src/ThetaSpanner.cpp
        src/ThetaSpanner.hpp
)

