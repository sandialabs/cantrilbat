include_directories(${CURRENT_SOURCE_DIR})
include_directories(${Trilinos_INCLUDE_DIRS})
include_directories(${CURRENT_SOURCE_DIR}/../src_CanteraElectrode)

file(GLOB SOURCES *.cpp)
file(GLOB HEADERS *.h)

add_executable(Cu2S_growth ${SOURCES} ${HEADERS})
target_link_libraries(Cu2S_growth ${Cantera_LIBRARIES} ${Trilinos_LIBRARIES} m1d m1dcantera util_src)

install(PROGRAMS Cu2S_growth
  DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
)
