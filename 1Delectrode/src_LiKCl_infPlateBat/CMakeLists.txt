include_directories(${CURRENT_SOURCE_DIR})
include_directories(${Trilinos_INCLUDE_DIRS})

file(GLOB SOURCES m1d*.cpp)
file(GLOB HEADERS *.h)

add_executable(LiKCl_infPlateBat ${SOURCES} ${HEADERS})
target_link_libraries(LiKCl_infPlateBat ${Cantera_LIBRARIES} ${Trilinos_LIBRARIES} m1d m1dcantera util_src)

install(PROGRAMS LiKCl_infPlateBat
  DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
)
