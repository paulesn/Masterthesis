cmake_minimum_required(VERSION 3.16)
project(Masterthesis)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_FLAGS "-O3 -Wall -Wextra")

# Find and link OpenMP
find_package(CGAL)
find_package(OpenMP)
find_package(Boost REQUIRED COMPONENTS serialization)

# Optional: add include directory if your headers are separated
# include_directories(include)
include_directories(src/daniel/theta-graph/headers)

add_executable(mtc
        src/io/mapToCoastlines.cc
)
add_executable(cdt6
        src/io/cdt6.cpp
)

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

add_executable(lowerboundTriang
        src/analysis/LowerBoundVergleichTriang.cpp
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
add_executable(lowerbound
        src/analysis/LowerBoundVergleich.cpp
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

add_executable(createVis
        src/spanner/createVisGraph.cpp

        src/structure/Graph.cpp
        src/structure/Graph.hpp

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
)

add_executable(cAllAnalysis
        src/analysis/CoastlineAllEdgesanalysis.cpp
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

add_executable(vsDanielAnalysis
        src/analysis/analysis_with_distance_file.cpp
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


add_executable(analyseCoastline
        src/analysis/CoastlineAnalysis.cpp
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

add_executable(analyseVis
        src/analysis/VisGraphAnalysis.cpp
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

add_executable(analysePaths
        src/analysis/PathAnalysis.cpp
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

add_executable(thetaComp
        src/ch_idea/theta_comparision.cpp
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

add_executable(ch
        src/ch_idea/CH_idea.cpp
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

add_executable(compareModify
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

        src/analysis/compareAlgoModifiedThetaGraph.cpp

        src/daniel/theta-graph/headers/Progressbar.h
        src/daniel/theta-graph/headers/Timer.h
        src/daniel/theta-graph/sources/Triangulation.cpp
        src/daniel/theta-graph/headers/Triangulation.h
        src/daniel/theta-graph/headers/Geometry.h
        src/daniel/theta-graph/sources/Geometry.cpp
        src/daniel/theta-graph/headers/Structs.h

        src/analysis/MultiThetaAnalysis.cpp
)

add_executable(Lowerbound
        src/daniel/Lowerbound_Graph.cpp
        src/daniel/Lowerbound_Graph.h
        src/structure/Graph.cpp
        src/structure/Graph.hpp
)

# Link OpenMP to your executable
target_link_libraries(analyse OpenMP::OpenMP_CXX)
target_link_libraries(analyse CGAL::CGAL)
target_link_libraries(analyse ${Boost_LIBRARIES})

target_link_libraries(analyseVis OpenMP::OpenMP_CXX)
target_link_libraries(analyseVis CGAL::CGAL)
target_link_libraries(analyseVis ${Boost_LIBRARIES})

target_link_libraries(mixWspdThetaSpanner OpenMP::OpenMP_CXX)
target_link_libraries(mixWspdThetaSpanner CGAL::CGAL)
target_link_libraries(mixWspdThetaSpanner ${Boost_LIBRARIES})

target_link_libraries(ch OpenMP::OpenMP_CXX)
target_link_libraries(ch CGAL::CGAL)
target_link_libraries(ch ${Boost_LIBRARIES})

target_link_libraries(analysePaths OpenMP::OpenMP_CXX)
target_link_libraries(analysePaths CGAL::CGAL)
target_link_libraries(analysePaths ${Boost_LIBRARIES})

target_link_libraries(thetaComp OpenMP::OpenMP_CXX)
target_link_libraries(thetaComp CGAL::CGAL)
target_link_libraries(thetaComp ${Boost_LIBRARIES})

target_link_libraries(analyseCoastline OpenMP::OpenMP_CXX)
target_link_libraries(analyseCoastline CGAL::CGAL)
target_link_libraries(analyseCoastline ${Boost_LIBRARIES})

target_link_libraries(cdt6 OpenMP::OpenMP_CXX)
target_link_libraries(cdt6 CGAL::CGAL)
target_link_libraries(cdt6 ${Boost_LIBRARIES})

target_link_libraries(lowerbound OpenMP::OpenMP_CXX)
target_link_libraries(lowerbound CGAL::CGAL)
target_link_libraries(lowerbound ${Boost_LIBRARIES})

target_link_libraries(lowerboundTriang OpenMP::OpenMP_CXX)
target_link_libraries(lowerboundTriang CGAL::CGAL)
target_link_libraries(lowerboundTriang ${Boost_LIBRARIES})

target_link_libraries(cAllAnalysis OpenMP::OpenMP_CXX)
target_link_libraries(cAllAnalysis CGAL::CGAL)
target_link_libraries(cAllAnalysis ${Boost_LIBRARIES})

target_link_libraries(createVis OpenMP::OpenMP_CXX)
target_link_libraries(createVis CGAL::CGAL)
target_link_libraries(createVis ${Boost_LIBRARIES})

target_link_libraries(vsDanielAnalysis OpenMP::OpenMP_CXX)
target_link_libraries(vsDanielAnalysis CGAL::CGAL)
target_link_libraries(vsDanielAnalysis ${Boost_LIBRARIES})