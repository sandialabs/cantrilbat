include_directories(${CURRENT_SOURCE_DIR})

file(GLOB SOURCES *.cpp)
file(GLOB HEADERS *.h)

add_library(util_src ${SOURCES} ${HEADERS})

install(target util_src
  DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
)
install(FILES ${HEADERS}
  DESTINATION ${CMAKE_INSTALL_PREFIX}/include
)
