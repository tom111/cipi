project(cipi)

add_executable(cipi main.cpp functions.cpp ConfigFile.cpp Marginal.cpp Marginals.cpp)

SET_TARGET_PROPERTIES(cipi PROPERTIES COMPILE_FLAGS "-Wall -O3 -fomit-frame-pointer -pipe")
FIND_PACKAGE(Boost REQUIRED)
IF (Boost_FOUND)
  Message("Boost found" ${Boost_LIBRARY_DIR})
  INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIR})
  #MESSAGE("Boost was found" ${Boost_INCLUDE_DIR})
  TARGET_LINK_LIBRARIES (cipi boost_thread-mt)
ENDIF (Boost_FOUND)

INSTALL(PROGRAMS cipi DESTINATION bin)
INSTALL(FILES PARAM DESTINATION share/doc/cipi-${CIPI_VERSION})