##
#  CMake script for the linear elastic principle program:
##


# Usually, you will not need to modify anything beyond this point...

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

set(QMAKE_CXX_FLAGS "${CMAKE_C_FLAGS} -Wno-deprecated")

FIND_PACKAGE(deal.II 8.5.0 
  HINTS ${deal.II_DIR} ${DEAL_II_DIR} ../ ../../ $ENV{DEAL_II_DIR}
  )
IF(NOT ${deal.II_FOUND})
  MESSAGE(FATAL_ERROR "\n"
    "*** Could not locate a (sufficiently recent) version of deal.II. ***\n\n"
    "You may want to either pass a flag -DDEAL_II_DIR=/path/to/deal.II to cmake\n"
    "or set an environment variable \"DEAL_II_DIR\" that contains this path."
    )
ENDIF()

DEAL_II_INITIALIZE_CACHED_VARIABLES()


 INCLUDE_DIRECTORIES(
  include
  ${CMAKE_SOURCE_DIR}/include
  )


ADD_LIBRARY(CompressedStripPacaBloch  ${CMAKE_SOURCE_DIR}/source/CompressedStripPacaBloch.cc)
DEAL_II_SETUP_TARGET(CompressedStripPacaBloch)

ADD_LIBRARY(Constituitive ${CMAKE_SOURCE_DIR}/source/Constituitive.cc)
DEAL_II_SETUP_TARGET(Constituitive)

# ADD_LIBRARY(mapping_q_eulerian_hp ${CMAKE_SOURCE_DIR}/include/mapping_q_eulerian_hp.h)
# DEAL_II_SETUP_TARGET(mapping_q_eulerian_hp)

TARGET_LINK_LIBRARIES(CompressedStripPacaBloch Constituitive)
# TARGET_LINK_LIBRARIES(CompressedStripPacaBloch mapping_q_eulerian_hp)


# ADD_EXECUTABLE(paca_to_first_bif_instability ${CMAKE_SOURCE_DIR}/source/paca_to_first_bif_instability.cc $<TARGET_OBJECTS:CompressedStripPacaBloch>)

ADD_EXECUTABLE(paca_to_first_bif_instability ${CMAKE_SOURCE_DIR}/source/paca_to_first_bif_instability.cc)
DEAL_II_SETUP_TARGET(paca_to_first_bif_instability)
TARGET_LINK_LIBRARIES(paca_to_first_bif_instability CompressedStripPacaBloch)

ADD_EXECUTABLE(paca_continue ${CMAKE_SOURCE_DIR}/source/paca_continue.cc)
DEAL_II_SETUP_TARGET(paca_continue)
TARGET_LINK_LIBRARIES(paca_continue CompressedStripPacaBloch)

ADD_EXECUTABLE(bloch_eigs_first_bif ${CMAKE_SOURCE_DIR}/source/bloch_eigs_first_bif.cc)
DEAL_II_SETUP_TARGET(bloch_eigs_first_bif)
TARGET_LINK_LIBRARIES(bloch_eigs_first_bif CompressedStripPacaBloch)

ADD_EXECUTABLE(find_first_bif_load ${CMAKE_SOURCE_DIR}/source/find_first_bif_load.cc)
DEAL_II_SETUP_TARGET(find_first_bif_load)
TARGET_LINK_LIBRARIES(find_first_bif_load CompressedStripPacaBloch)

ADD_EXECUTABLE(bloch_eigs_branch ${CMAKE_SOURCE_DIR}/source/bloch_eigs_branch.cc)
DEAL_II_SETUP_TARGET(bloch_eigs_branch)
TARGET_LINK_LIBRARIES(bloch_eigs_branch CompressedStripPacaBloch)

ADD_EXECUTABLE(compare_solutions ${CMAKE_SOURCE_DIR}/source/compare_solutions.cc)
DEAL_II_SETUP_TARGET(compare_solutions)
TARGET_LINK_LIBRARIES(compare_solutions CompressedStripPacaBloch)

ADD_EXECUTABLE(princple_info ${CMAKE_SOURCE_DIR}/source/princple_info.cc)
DEAL_II_SETUP_TARGET(princple_info)
TARGET_LINK_LIBRARIES(princple_info CompressedStripPacaBloch)


# Define a distclean target to remove every generated file:
ADD_CUSTOM_TARGET(distclean
  COMMAND ${CMAKE_COMMAND} --build ${CMAKE_BINARY_DIR} --target clean
  COMMAND ${CMAKE_COMMAND} -E remove_directory CMakeFiles
  COMMAND ${CMAKE_COMMAND} -E remove CMakeCache.txt cmake_install.cmake Makefile
  COMMENT "distclean invoked"
  )


