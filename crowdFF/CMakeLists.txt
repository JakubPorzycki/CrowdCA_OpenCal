cmake_minimum_required(VERSION 2.8 FATAL_ERROR)

project(crowdFF)

find_package(OpenCAL 1.0 COMPONENTS OpenCAL OpenCAL-GL REQUIRED)
find_package(OpenGL  REQUIRED)
find_package(GLUT REQUIRED)
find_package(Threads REQUIRED)

if(NOT MSVC)
  find_library(M_LIB m)
else()
  # not needed on windows
  set(M_LIB "")
endif()

FIND_PACKAGE(OpenMP)
if(OPENMP_FOUND)
    SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif()

include_directories(
	${CMAKE_CURRENT_SOURCE_DIR}/include
	${OPENGL_INCLUDE_DIRS}
	${GLUT_INCLUDE_DIR}
	${OPENCAL_INCLUDE_DIR}
)

set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin)

link_directories(${CMAKE_CURRENT_SOURCE_DIR}/source)

file(GLOB SRC
    "include/*.h"
    "source/*.cpp"
)

add_executable(crowdFF ${SRC})
target_link_libraries(crowdFF ${OPENCAL_LIBRARIES} ${M_LIB} ${OPENGL_LIBRARIES} ${GLUT_glut_LIBRARY} pthread)
