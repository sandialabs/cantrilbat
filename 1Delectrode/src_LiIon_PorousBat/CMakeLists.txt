include_directories(${CURRENT_SOURCE_DIR})
include_directories(${Trilinos_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/1Delectrode/src)
file(GLOB SOURCES *.cpp)
file(GLOB HEADERS *.h)

add_definitions(-DINTEGRATED_ELECTRODE -DDELTASHEAT_ZERO)

add_executable(LiIon_PorousBat ${SOURCES} ${HEADERS})
target_link_libraries(LiIon_PorousBat ${Cantera_LIBRARIES} ${Trilinos_LIBRARIES} m1d m1dcantera util_src)

install(PROGRAMS LiIon_PorousBat
  DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
)
