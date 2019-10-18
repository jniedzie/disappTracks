cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(event)

list(APPEND CMAKE_PREFIX_PATH $ENV{ROOTSYS})

find_package(ROOT REQUIRED COMPONENTS Minuit Eve Geom)
include(${ROOT_USE_FILE})

SET(CMAKE_SKIP_BUILD_RPATH FALSE)
SET(CMAKE_BUILD_WITH_INSTALL_RPATH TRUE) 
SET(CMAKE_INSTALL_RPATH "/cvmfs/sft.cern.ch/lcg/releases/vdt/0.4.3-992df/x86_64-centos7-gcc8-opt/lib:/cvmfs/sft.cern.ch/lcg/releases/tbb/2019_U7-ba3eb/x86_64-centos7-gcc8-opt/lib:/cvmfs/sft.cern.ch/lcg/releases/Davix/0.7.3-d94fa/x86_64-centos7-gcc8-opt/lib64/")
SET(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

add_executable(getTrackHeatMap getTrackHeatMap.cpp)
add_executable(eventDisplay eventDisplay.cpp)
add_executable(prepareABCDplots prepareABCDplots.cpp)
add_executable(helixTagger helixTagger.cpp)
add_executable(analyzeTagger analyzeTagger.cpp)
add_executable(taggerScanner taggerScanner.cpp)
add_executable(taggerFitter taggerFitter.cpp)
add_executable(analyzeClusters analyzeClusters.cpp)
add_executable(runAnalysis runAnalysis.cpp)
add_executable(drawDistributions drawDistributions.cpp)
add_executable(scanCuts scanCuts.cpp)
add_executable(singleCutDetails singleCutDetails.cpp)
add_executable(cutsCorrelations cutsCorrelations.cpp)
add_executable(getFfactor getFfactor.cpp)

set_property(TARGET getTrackHeatMap PROPERTY CXX_STANDARD 17)
set_property(TARGET eventDisplay PROPERTY CXX_STANDARD 17)
set_property(TARGET prepareABCDplots PROPERTY CXX_STANDARD 17)
set_property(TARGET helixTagger PROPERTY CXX_STANDARD 17)
set_property(TARGET analyzeTagger PROPERTY CXX_STANDARD 17)
set_property(TARGET taggerScanner PROPERTY CXX_STANDARD 17)
set_property(TARGET taggerFitter PROPERTY CXX_STANDARD 17)
set_property(TARGET analyzeClusters PROPERTY CXX_STANDARD 17)
set_property(TARGET runAnalysis PROPERTY CXX_STANDARD 17)
set_property(TARGET drawDistributions PROPERTY CXX_STANDARD 17)
set_property(TARGET scanCuts PROPERTY CXX_STANDARD 17)
set_property(TARGET singleCutDetails PROPERTY CXX_STANDARD 17)
set_property(TARGET cutsCorrelations PROPERTY CXX_STANDARD 17)
set_property(TARGET getFfactor PROPERTY CXX_STANDARD 17)

include_directories("${PROJECT_SOURCE_DIR}/include")

add_library(libAnalysis src/Track.cpp
                        src/TrackCut.cpp
                        src/TrackProcessor.cpp
                        src/Jet.cpp
                        src/JetCut.cpp
                        src/JetProcessor.cpp
                        src/Lepton.cpp
                        src/LeptonCut.cpp
                        src/LeptonProcessor.cpp
                        src/Event.cpp
                        src/EventCut.cpp
                        src/EventProcessor.cpp
                        src/EventSet.cpp
                        src/ConfigManager.cpp
                        src/CutsManager.cpp
                        src/Display.cpp
                        src/Fitter.cpp
                        src/Helix.cpp
                        src/HelixProcessor.cpp
                        src/HistSet.cpp
                        src/MonitorsManager.cpp
                        src/PerformanceMonitor.cpp
                        src/Point.cpp
                        src/PointsProcessor.cpp
)

set_property(TARGET libAnalysis PROPERTY CXX_STANDARD 17)

target_link_libraries(getTrackHeatMap libAnalysis ${ROOT_LIBRARIES})
target_link_libraries(eventDisplay libAnalysis ${ROOT_LIBRARIES})
target_link_libraries(prepareABCDplots libAnalysis ${ROOT_LIBRARIES})
target_link_libraries(helixTagger libAnalysis ${ROOT_LIBRARIES})
target_link_libraries(analyzeTagger libAnalysis ${ROOT_LIBRARIES})
target_link_libraries(taggerScanner libAnalysis ${ROOT_LIBRARIES})
target_link_libraries(taggerFitter libAnalysis ${ROOT_LIBRARIES})
target_link_libraries(analyzeClusters libAnalysis ${ROOT_LIBRARIES})
target_link_libraries(runAnalysis libAnalysis ${ROOT_LIBRARIES})
target_link_libraries(drawDistributions libAnalysis ${ROOT_LIBRARIES})
target_link_libraries(scanCuts libAnalysis ${ROOT_LIBRARIES})
target_link_libraries(singleCutDetails libAnalysis ${ROOT_LIBRARIES})
target_link_libraries(cutsCorrelations libAnalysis ${ROOT_LIBRARIES})
target_link_libraries(getFfactor libAnalysis ${ROOT_LIBRARIES})
