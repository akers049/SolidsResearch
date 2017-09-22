##
#  CMake script for the linear elastic principle program:
##


# Usually, you will not need to modify anything beyond this point...

CMAKE_MINIMUM_REQUIRED(VERSION 2.8.8)

FIND_PACKAGE(deal.II 8.5.0 QUIET
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

ADD_LIBRARY(CompressedStripPacaBloch source/CompressedStripPacaBloch.cc)
DEAL_II_SETUP_TARGET(CompressedStripPacaBloch)

ADD_EXECUTABLE(paca_to_first_bif_instability ${CMAKE_SOURCE_DIR}/source/paca_to_first_bif_instability.cc)
DEAL_II_SETUP_TARGET(paca_to_first_bif_instability)
TARGET_LINK_LIBRARIES(paca_to_first_bif_instability CompressedStripPacaBloch)

ADD_EXECUTABLE(paca_continue ${CMAKE_SOURCE_DIR}/source/paca_continue.cc)
DEAL_II_SETUP_TARGET(paca_continue)
TARGET_LINK_LIBRARIES(paca_continue CompressedStripPacaBloch)


