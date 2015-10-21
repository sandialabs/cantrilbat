include_directories(${CURRENT_SOURCE_DIR})
include_directories(${Trilinos_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/1Delectrode/src)
file(GLOB SOURCES *.cpp)
file(GLOB HEADERS *.h)

add_executable(LiKCl_infPorousBat ${SOURCES} ${HEADERS})
target_link_libraries(LiKCl_infPorousBat ${Cantera_LIBRARIES} ${Trilinos_LIBRARIES} m1d m1dcantera util_src)

install(PROGRAMS LiKCl_infPorousBat
  DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
)
