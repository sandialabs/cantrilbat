include_directories(${CMAKE_CURRENT_SOURCE_DIR})
include_directories(${Trilinos_INCLUDE_DIRS})
include_directories(${Cantera_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/util_src/src)
include_directories(${PROJECT_SOURCE_DIR}/Electrode/src)


file(GLOB M1D_SOURCES m1d*.cpp)
list(REMOVE_ITEM M1D_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/m1d_nox_interface.cpp")
list(REMOVE_ITEM M1D_SOURCES "${CMAKE_CURRENT_SOURCE_DIR}/m1d_nox_problem.cpp")
list(APPEND M1D_SOURCES md_wrap_mpi.cpp)

file(GLOB M1D_HEADERS m1d*.h)
list(REMOVE_ITEM M1D_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/m1d_nox_interface.h")
list(REMOVE_ITEM M1D_HEADERS "${CMAKE_CURRENT_SOURCE_DIR}/m1d_nox_problem.h")
list(APPEND  M1D_HEADERS md_wrap_mpi.h)

set(BEuler_SOURCES BEulerInt.cpp)
set(BEuler_HEADERS BEulerInt.h)

add_library(m1d ${M1D_SOURCES} ${M1D_HEADERS} ${BEuler_SOURCES} ${BEuler_HEADERS}  )
target_link_libraries(m1d ${Trilinos_LIBRARIES} )
target_link_libraries(m1d ${Cantera_LIBRARIES} )
target_link_libraries(m1d ${util_LIBRARIES} )
target_link_libraries(m1d ${electrode_LIBRARIES} )


install(TARGETS m1d
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
install(FILES ${M1D_HEADERS} ${BEuler_HEADERS}
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include)

add_executable(Multi1DDomain Multi1DDomain.cpp)
target_link_libraries(Multi1DDomain m1d ${util_LIBRARIES}  ${electrode_LIBRARIES} ${Trilinos_LIBRARIES} ${Cantera_LIBRARIES}  )

install(PROGRAMS Multi1DDomain
  DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)

add_executable(tddiff tddiff.cpp)
target_link_libraries(tddiff m1d)

install(PROGRAMS tddiff
  DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
