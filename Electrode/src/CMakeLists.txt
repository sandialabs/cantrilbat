include_directories(${CURRENT_SOURCE_DIR})
include_directories(${Cantera_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/util_src/src)
include_directories(${PROJECT_SOURCE_DIR}/ctpadds/src)
include_directories("/usr/netpub/sundials/sundials-2.5.0/include/")

file(GLOB SOURCES *.cpp)
list(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/Electrode_SimpleDiffCVode.cpp")
list(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/Electrode_SimpleDiff_SEI.cpp")
list(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/Electrode_PhaseStability.cpp")
list(REMOVE_ITEM SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/ctmlDiff.cpp")

file(GLOB HEADERS *.h)
list(REMOVE_ITEM SOURCES "${PROJECT_SOURCE_DIR}/Electrode/src/Electrode_PhaseStability.cpp")
list(REMOVE_ITEM SOURCES "${PROJECT_SOURCE_DIR}/Electrode/src/Electrode_SimpleDiffCVode.cpp")
list(REMOVE_ITEM SOURCES "${PROJECT_SOURCE_DIR}/Electrode/src/Electrode_SimpleDiff_SEI.cpp")
list(REMOVE_ITEM SOURCES "${PROJECT_SOURCE_DIR}/Electrode/src/ctmlDiff.cpp")

add_library(electrode_LIBRARIES ${SOURCES} ${HEADERS})

#target_link_libraries(electrode_LIBRARIES util_LIBRARIES)
#target_link_libraries(electrode_LIBRARIES ctpadds_LIBRARIES)
target_link_libraries(electrode_LIBRARIES ${Cantera_LIBRARIES} ${Trilinos_LIBRARIES} util_LIBRARIES )
target_link_libraries(electrode_LIBRARIES "-L/usr/netpub/sundials/sundials-2.5.0/lib")

# \todo 
#Add executable ctmlDiff
#

install(TARGETS electrode_LIBRARIES
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)
install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include
)
