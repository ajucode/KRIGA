#------------------------------------------------------------------------------
# Kriging
# author: Pedro Guarderas
# date: 30-01-2013
#------------------------------------------------------------------------------

cmake_minimum_required( VERSION 2.8 )
set( CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR}/cmake ${CMAKE_MODULE_PATH})

# LLVM and clang previous configuration
#INCLUDE (CMakeForceCompiler)
#SET (CMAKE_SYSTEM_NAME Generic )		
#	SET(CMAKE_AR "/usr/bin/llvm-ar" )
#	SET(CMAKE_LINKER "/usr/bin/llvm-ld" )
#	SET(CMAKE_RANLIB "/usr/bin/llvm-ranlib" )
#	SET(CMAKE_NM "/usr/bin/llvm-nm" )
#	SET(CMAKE_OBJDUMP "/usr/bin/llvm-objdump" )

#FIND_PACKAGE( LLVM )

#IF ( LLVM_FOUND )
#
#	SET (CMAKE_C_COMPILER             ${CLANG} )
#        SET (CMAKE_C_FLAGS                "-Wall")
#	SET (CMAKE_C_FLAGS_DEBUG          "-g")
#	SET (CMAKE_C_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
#	SET (CMAKE_C_FLAGS_RELEASE        "-O4 -DNDEBUG")
#	SET (CMAKE_C_FLAGS_RELWITHDEBINFO "-O2 -g")

#	SET (CMAKE_CXX_COMPILER             ${CLANG_PLUS} )
#        SET (CMAKE_CXX_FLAGS                "-Wall -std=c++11")
#	SET (CMAKE_CXX_FLAGS_DEBUG          "-g")
#	SET (CMAKE_CXX_FLAGS_MINSIZEREL     "-Os -DNDEBUG")
#	SET (CMAKE_CXX_FLAGS_RELEASE        "-O4 -DNDEBUG")
#	SET (CMAKE_CXX_FLAGS_RELWITHDEBINFO "-O2 -g")

	#SET(CMAKE_C_LINKER ${LLVM_LD} )
	#SET(CMAKE_CXX_LINKER ${LLVM_LD} )

	#SET(CMAKE_C_LINK_EXECUTABLE "/usr/bin/llvm-ld <OBJECTS> -o  <TARGET> <CMAKE_C_LINK_FLAGS> <LINK_FLAGS> <LINK_LIBRARIES>")
	#SET(CMAKE_CXX_LINK_EXECUTABLE "/usr/bin/llvm-ld <OBJECTS> -o  <TARGET> <CMAKE_CXX_LINK_FLAGS> <LINK_FLAGS> <LINK_LIBRARIES>")

#ELSE ( LLVM_FOUND )
#	MESSAGE ( SEND_ERROR "LLVM don't work" )

#ENDIF ( LLVM_FOUND )

#------------------------------------------------------------------------------
# Project definition
set(CPACK_PACKAGE_VERSION_MAJOR "0")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "1")

set(COMPLETE_VERSION ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH})

set(PROJECT_VERSION ${COMPLETE_VERSION})
project(kriging{PROJECT_VERSION})
message(STATUS "Kriging version: ${COMPLETE_VERSION}")

set(CMAKE_COLOR_MAKEFILE ON)
option(KRIGING_ENABLE_TESTS "Enable testing" OFF )
#CMAKE_FORCE_C_COMPILER ( /usr/bin/clang clang )
#CMAKE_FORCE_CXX_COMPILER ( /usr/bin/clang++ clang_plus )
add_definitions(-std=c++11)

#------------------------------------------------------------------------------
# Search MPI
find_package(MPI REQUIRED)
if(MPI_FOUND)
  set(CMAKE_C_COMPILER ${MPI_C_COMPILER} )
  set(CMAKE_CXX_COMPILER ${MPI_CXX_COMPILER} )
  set(CMAKE_Fortran_COMPILER ${MPI_Fortran_COMPILER} )
endif(MPI_FOUND)

#------------------------------------------------------------------------------
# Search GSL
find_package( GSL REQUIRED )
if( GSL_FOUND )
  #add_definitions("-DHAVE_GSL")
endif( GSL_FOUND )

#------------------------------------------------------------------------------
# Search Global Arrays
find_package( GA REQUIRED )
if( GA_FOUND )
#  message( SEND_ERROR "${GA_INCLUDE_DIRS}" )
#  message( SEND_ERROR "${GA_LIBRARIES}" )
#  message( SEND_ERROR "${GA_LIBRARY_DIRS}" )
#  message( SEND_ERROR "${GSL_LIBRARY_DIRS}" )
#  message( SEND_ERROR "${GSL_LIBRARIES}" )
endif( GA_FOUND )

#------------------------------------------------------------------------------
# Search OpenCL
#find_package(OpenCL REQUIRED)

#------------------------------------------------------------------------------
# Search Boost library
find_package(Boost)
if(Boost_FOUND)
#  include_directories(${Boost_INCLUDE_DIR})
  add_definitions( "-DHAS_BOOST" )
endif(Boost_FOUND)

#------------------------------------------------------------------------------
# Search for Qt4
#set(QT_MIN_VERSION 4.4.0)
#find_package(Qt4 ${QT_MIN_VERSION} REQUIRED)
#set( QT_USE_QTOPENGL TRUE )
#set( QT_USE_QTXML TRUE )
#include( ${QT_USE_FILE} )

#------------------------------------------------------------------------------
if(WIN32)
  if(MSVC)
    add_definitions("-DGUI_EXPORT=__declspec(dllimport)")
    add_definitions("-DCORE_EXPORT=__declspec(dllimport)")
  else(MSVC)
    add_definitions("\"-DGUI_EXPORT=__declspec(dllimport)\"")
    add_definitions("\"-DCORE_EXPORT=__declspec(dllimport)\"")
  endif(MSVC)
else(WIN32)
  add_definitions(-DGUI_EXPORT=)
  add_definitions(-DCORE_EXPORT=)
endif(WIN32)

if(CMAKE_BUILD_TYPE MATCHES Debug)
  add_definitions(-DkrigingDEBUG=1)
endif(CMAKE_BUILD_TYPE MATCHES Debug)

#------------------------------------------------------------------------------
set(DEFAULT_BIN_SUBDIR     bin)
set(DEFAULT_LIB_SUBDIR     lib)
set(DEFAULT_INCLUDE_SUBDIR include)
set(DEFAULT_DATA_SUBDIR    share)

#------------------------------------------------------------------------------
# Including directories
#include_directories(
#  ${CMAKE_CURRENT_BINARY_DIR}
#  ${CMAKE_CURRENT_SOURCE_DIR}
#  ${GSL_INCLUDE_DIR}
#  ${MPI_INCLUDE_PATH}
#  ${OPENCL_INCLUDE_DIRS}
#  ${QT_INCLUDE_DIR}
#  .
#)

#------------------------------------------------------------------------------
# Subdirectories
add_subdirectory( src )
if( KRIGING_ENABLE_TESTS )
  message(STATUS "Enabling tests")
  enable_testing()
  add_subdirectory( test )
endif( KRIGING_ENABLE_TESTS )

