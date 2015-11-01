include_directories(${CURRENT_SOURCE_DIR})
include_directories(${PROJECT_SOURCE_DIR}/1Delectrode/src)
include_directories(${PROJECT_SOURCE_DIR}/1Delectrode/src_CanteraElectrode)
include_directories(${PROJECT_SOURCE_DIR}/util_src/src)
include_directories(${PROJECT_SOURCE_DIR}/Electrode/src)
include_directories(${Trilinos_INCLUDE_DIRS})

file(GLOB SOURCES *.cpp)
file(GLOB HEADERS *.h)

add_definitions(-DINTEGRATED_ELECTRODE -DDELTASHEAT_ZERO -DMECH_MODEL -Woverloaded-virtual )

add_executable(LiIon_PorousBat ${SOURCES} ${HEADERS})
target_link_libraries(LiIon_PorousBat ${Trilinos_LIBRARIES} m1d m1dcantera util_src electrode_LIBRARIES )

install(PROGRAMS LiIon_PorousBat
  DESTINATION ${CMAKE_INSTALL_PREFIX}/bin
)
