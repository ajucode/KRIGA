#---------------------------------------------------------------------------------------------------
# Sources
set( KRIGA_SRCS
  kriga.cpp )  

set( KRIGA_HDRS
  distance.hpp
  kernels.hpp
  kriga.hpp
)

set( KRIGA_UIS )

set( KRIGA_MOC_HDRS )

set( KRIGA_RCCS )

#---------------------------------------------------------------------------------------------------
# Including directories
include_directories(
  ${Boost_INCLUDE_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_CURRENT_SOURCE_DIR}
  ${GA_INCLUDE_DIRS}
  ${GSL_INCLUDE_DIRS}
  ${GSLCBLAS_INCLUDE_DIRS}
  ${MPI_INCLUDE_PATH}
  ..
  . 
)

#link_directories(
#  ${KRIGA_LIBRARIES}
#  ${MPI_LIBRARIES} 
#  . 
#)

#---------------------------------------------------------------------------------------------------  
# Creating libray kriga
# ${KRIGA_MOC_SRCS} ${KRIGA_UIS_H} ${KRIGA_RCC_SRCS} 
add_library( Kriga STATIC ${KRIGA_SRCS} ${KRIGA_HDRS} )

set_target_properties( Kriga PROPERTIES 
  LINKER_LANGUAGE CXX )

target_link_libraries( Kriga
  ${GA_LIBRARIES}
  ${GSL_LIBRARIES}
  ${MPI_LIBRARIES} )

#---------------------------------------------------------------------------------------------------
# Install
install( TARGETS Kriga 
#  RUNTINE DESTINATION ${CMAKE_INSTALL_PREFIX}
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}/lib
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_PREFIX}/include
)
