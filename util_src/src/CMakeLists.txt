include_directories(${CURRENT_SOURCE_DIR})
include_directories(${Cantera_INCLUDE_DIRS})

file(GLOB SOURCES *.cpp)
file(GLOB HEADERS *.h)

add_library(util_LIBRARIES ${SOURCES} ${HEADERS})
target_link_libraries(util_LIBRARIES  ${Cantera_LIBRARIES} ${Trilinos_LIBRARIES})

install(TARGETS util_LIBRARIES
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)
install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include
)

