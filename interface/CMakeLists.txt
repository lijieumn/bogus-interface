project(ArgusInterface)
cmake_minimum_required(VERSION 2.8)

SET( ARGUS_ROOT ${CMAKE_CURRENT_SOURCE_DIR}/.. )
SET( ARGUS_INTERFACE_LIBNAME ${CMAKE_PROJECT_NAME}_bogus )
SET( ARGUS_LOADER_NAME       ${CMAKE_PROJECT_NAME}_bogusLoader )

####
# Libraries
####

SET( Eigen3_DIR ${ARGUS_ROOT}/extern/eigen/cmake)

# Eigen
FIND_PACKAGE(Eigen3 REQUIRED)
SET( EXTERNAL_INCLUDES ${EIGEN3_INCLUDE_DIR} )
SET( EXTERNAL_LIBRARIES "" )

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")

# OpenMP
FIND_PACKAGE(OpenMP)
if(OPENMP_FOUND)
	set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
	add_definitions(-DEIGEN_DONT_PARALLELIZE)
else()
	set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
	#set( CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-unknown-pragmas " )
endif()

# boost
find_package(Boost 1.36.0 COMPONENTS serialization REQUIRED )
SET( EXTERNAL_INCLUDES ${EXTERNAL_INCLUDES} ${Boost_INCLUDE_DIR} )
SET( EXTERNAL_LIBRARIES ${EXTERNAL_LIBRARIES} ${Boost_LIBRARIES} )


#Bogus
SET( SoBogus_INCLUDE_DIR ${ARGUS_ROOT}/extern/so-bogus )
SET( EXTERNAL_INCLUDES ${EXTERNAL_INCLUDES} ${SoBogus_INCLUDE_DIR} )
add_definitions( -DBOGUS_WITH_BOOST_SERIALIZATION )

####
# Compiler flags
####

if (NOT CMAKE_BUILD_TYPE)
  message(STATUS "No build type selected, default to Release")
  set(CMAKE_BUILD_TYPE "Release")
endif()

if ( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang" )
  set (CLANG TRUE)
  include_directories( /usr/lib/clang/3.7/include )
  set( CMAKE_EXPORT_COMPILE_COMMANDS 1 )
elseif( "${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU" )
  set (GCC TRUE)
endif()

if ( GCC OR CLANG )
				SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra -std=c++11 " )
				SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-deprecated-declarations " )
				SET(CMAKE_CXX_FLAGS_DEBUG "${CMAKE_CXX_FLAGS_DEBUG} -Werror " )

				if( GCC AND APPLE )
						  SET( EXTERNAL_LIBRARIES ${EXTERNAL_LIBRARIES} stdc++ )
				endif()

				#SET(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE}  -march=native -ffast-math -DEIGEN_NO_DEBUG ")
endif()

####
# Sources and targets
####

#Lib 
file(GLOB SRC_LIST src/*.hh src/*.cc)
file(GLOB HEADERS src/*.hh)

add_library(${ARGUS_INTERFACE_LIBNAME} SHARED ${SRC_LIST})
target_include_directories(${ARGUS_INTERFACE_LIBNAME} BEFORE
	  PUBLIC ${EXTERNAL_INCLUDES} )
target_link_libraries(${ARGUS_INTERFACE_LIBNAME} ${EXTERNAL_LIBRARIES} )

# Install
INSTALL( FILES ${HEADERS} DESTINATION include/${ARGUS_INTERFACE_LIBNAME})
INSTALL( TARGETS ${ARGUS_INTERFACE_LIBNAME}
  LIBRARY DESTINATION lib
  RUNTIME DESTINATION lib
  ARCHIVE DESTINATION lib )

# Friction problem loader
add_executable(${ARGUS_LOADER_NAME} ${CMAKE_CURRENT_SOURCE_DIR}/test/loader.cpp )
target_include_directories(${ARGUS_LOADER_NAME} BEFORE
	  PUBLIC ${EXTERNAL_INCLUDES} ${CMAKE_CURRENT_SOURCE_DIR}/src )
target_link_libraries(${ARGUS_LOADER_NAME} ${ARGUS_INTERFACE_LIBNAME} )


