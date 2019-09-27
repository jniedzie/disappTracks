cmake_minimum_required (VERSION 2.6)
project (DisappTracks)

find_package(ROOT REQUIRED COMPONENTS Minuit Eve Geom)
include(${ROOT_USE_FILE})

set_property(TARGET ROOT::Core PROPERTY
INTERFACE_INCLUDE_DIRECTORIES "${ROOT_INCLUDE_DIRS}")

add_library(ROOT::Flags_CXX IMPORTED INTERFACE)

# Fix for ROOT_CXX_FLAGS not actually being a CMake list
separate_arguments(ROOT_CXX_FLAGS)
set_property(TARGET ROOT::Flags_CXX APPEND PROPERTY
INTERFACE_COMPILE_OPTIONS ${ROOT_CXX_FLAGS})

# Add definitions
separate_arguments(ROOT_DEFINITIONS)
foreach(_flag ${ROOT_EXE_LINKER_FLAG_LIST})
# Remove -D or /D if present
string(REGEX REPLACE [=[^[-//]D]=] "" _flag ${_flag})
set_property(TARGET ROOT::Flags APPEND PROPERTY INTERFACE_LINK_LIBRARIES ${_flag})
endforeach()

# This also fixes a bug in the linker flags
string(REPLACE "-L " "-L" ROOT_EXE_LINKER_FLAGS "${ROOT_EXE_LINKER_FLAGS}")
separate_arguments(ROOT_EXE_LINKER_FLAGS)
set_property(TARGET ROOT::Flags_CXX APPEND PROPERTY
INTERFACE_LINK_LIBRARIES ${ROOT_EXE_LINKER_FLAGS})

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/Modules/")

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
