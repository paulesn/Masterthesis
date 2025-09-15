cmake_minimum_required(VERSION 3.16)
project(Masterthesis)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Find and link OpenMP
find_package(CGAL)
find_package(OpenMP)
find_package(Boost REQUIRED COMPONENTS serialization)

# Optional: add include directory if your headers are separated
# include_directories(include)
include_directories(src/daniel/theta-graph/headers)

add_executable(analyse
        src/analysis/ThetaSpannerAnalyser.cpp
        src/analysis/MultiThetaAnalysis.cpp

        src/daniel/theta-graph/headers/Progressbar.h
        src/daniel/theta-graph/headers/Timer.h
        src/daniel/theta-graph/sources/Triangulation.cpp
        src/daniel/theta-graph/headers/Triangulation.h
        src/daniel/theta-graph/headers/Geometry.h
        src/daniel/theta-graph/sources/Geometry.cpp
        src/daniel/theta-graph/headers/Structs.h

        src/io/DataWriter.cpp
        src/io/DataWriter.h
        src/io/Dataloader.cpp
        src/io/Dataloader.hpp

        src/spanner/ThetaSpanner.hpp
        src/spanner/ThetaSpanner.cpp

        src/structure/Graph.cpp
        src/structure/Graph.hpp



        src/analysis/MultiThetaAnalysis.h


)

add_executable(mixWspdThetaSpanner
        src/io/DataWriter.cpp
        src/io/DataWriter.h
        src/io/Dataloader.cpp
        src/io/Dataloader.hpp

        src/spanner/ThetaSpanner.hpp
        src/spanner/ThetaSpanner.cpp

        src/structure/Graph.cpp
        src/structure/Graph.hpp
        src/structure/Quadtree.cpp
        src/structure/Quadtree.hpp
        src/spanner/MixWspdThetaSpannerGenerator.cpp
)

# Link OpenMP to your executable
target_link_libraries(analyse OpenMP::OpenMP_CXX)
target_link_libraries(analyse CGAL::CGAL)
target_link_libraries(analyse ${Boost_LIBRARIES})

target_link_libraries(mixWspdThetaSpanner OpenMP::OpenMP_CXX)
