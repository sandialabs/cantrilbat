include_directories(${CURRENT_SOURCE_DIR})
include_directories(${Cantera_INCLUDE_DIRS})


file(GLOB SOURCES *.cpp)
file(GLOB HEADERS *.h)

add_library(electrode_LIBRARIES ${SOURCES} ${HEADERS})

install(TARGETS electrode_LIBRARIES
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)
install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include
)
