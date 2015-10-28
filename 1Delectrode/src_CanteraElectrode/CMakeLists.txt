include_directories(${CURRENT_SOURCE_DIR})
include_directories(${Trilinos_INCLUDE_DIRS})
include_directories(${Cantera_INCLUDE_DIRS})
include_directories(${PROJECT_SOURCE_DIR}/1Delectrode/src)
include_directories(${PROJECT_SOURCE_DIR}/util_src/src)
include_directories(${PROJECT_SOURCE_DIR}/Electrode/src)

file(GLOB M1DCantera_SOURCES m1d*.cpp)
list(APPEND M1DCantera_SOURCES BEulerInt_Battery.cpp)
list(APPEND M1DCantera_SOURCES FuncElectrodeCurrent.cpp)

file(GLOB M1DCantera_HEADERS m1d*.hpp)
list(APPEND M1DCantera_HEADERS BEulerInt_Battery.h)
list(APPEND M1DCantera_HEADERS FuncElectrodeCurrent.h)

add_library(m1dcantera ${M1DCantera_SOURCES} ${M1DCantera_HEADERS})
target_link_libraries(m1dcantera ${Cantera_LIBRARIES} ${Trilinos_LIBRARIES} util_src)

install(TARGETS m1dcantera
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib)
install(FILES ${M1DCantera_HEADERS}
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include)


