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


add_executable(analyse
        src/analysis/ThetaSpannerAnalyser.cpp
        src/analysis/MultiThetaAnalysis.cpp

        src/daniel/Progressbar.h
        src/daniel/Timer.h
        src/daniel/Triangulation.cpp
        src/daniel/Triangulation.h
        src/daniel/Geometry.h
        src/daniel/Geometry.cpp
        src/daniel/Structs.h

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

# Link OpenMP to your executable
target_link_libraries(analyse OpenMP::OpenMP_CXX)
target_link_libraries(analyse CGAL::CGAL)
target_link_libraries(analyse ${Boost_LIBRARIES})
