cmake_minimum_required(VERSION 3.16)
project(Masterthesis)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find and link OpenMP
find_package(OpenMP REQUIRED)

# Optional: add include directory if your headers are separated
# include_directories(include)

add_executable(main
        src/main.cpp
        src/Graph.cpp
        src/Quadtree.cpp
        src/Dataloader.cpp
        src/Dataloader.hpp
        src/DataWriter.cpp
        src/DataWriter.h
)

add_executable(sizeTest
        test/sizeTest.cpp
        src/Graph.cpp
        src/Quadtree.cpp
        src/Dataloader.cpp
        src/Dataloader.hpp
        src/ThetaSpanner.cpp
        src/HubLabels.hpp
        src/HubLabels.cpp
        src/DataWriter.cpp
        src/DataWriter.h
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
        src/DataWriter.cpp
        src/DataWriter.h
)

add_executable(GraphTest
        src/Graph.cpp
        src/Dataloader.cpp
        src/Dataloader.hpp
        src/HubLabels.hpp
        src/HubLabels.cpp
        src/Quadtree.hpp
        src/Quadtree.cpp
        test/GraphTest.cpp
        src/DataWriter.cpp
        src/DataWriter.h
)

add_executable(Theta
        src/Graph.cpp
        src/Dataloader.cpp
        src/Dataloader.hpp
        src/ThetaSpanner.cpp
        src/ThetaSpanner.hpp
        src/DataWriter.cpp
        src/DataWriter.h
)

add_executable(ThetaSpannerTest
        test/thetaSpannerTest.cpp
        src/ThetaSpanner.cpp
        src/ThetaSpanner.hpp
        src/Graph.cpp
        src/HubLabels.cpp
)


# Link OpenMP to your executable
target_link_libraries(WSPDTest PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(sizeTest PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(GraphTest PRIVATE OpenMP::OpenMP_CXX)
target_link_libraries(ThetaSpannerTest PRIVATE OpenMP::OpenMP_CXX)

