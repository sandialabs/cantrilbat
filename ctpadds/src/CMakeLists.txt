include_directories(${CURRENT_SOURCE_DIR})
#include_directories(${Cantera_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/util_src/src)

file(GLOB SOURCES *.cpp)
list(REMOVE_ITEM SOURCES  "${CMAKE_CURRENT_SOURCE_DIR}/ElectrodeKinetics_old.cpp")
file(GLOB HEADERS *.h)

add_library(ctpadds_LIBRARIES ${SOURCES} ${HEADERS})
target_link_libraries(ctpadds_LIBRARIES  ${Trilinos_LIBRARIES} ${util_LIBRARIES})

install(TARGETS ctpadds_LIBRARIES
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)
install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include
)
