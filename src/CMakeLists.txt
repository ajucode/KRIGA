#---------------------------------------------------------------------------------------------------
# Sources
set( KRIGING_SRCS
  kriging.cpp )  

set( KRIGING_HDRS
  distance.hpp
  kernels.hpp
  kriging.hpp
)

set( KRIGING_UIS )

set( KRIGING_MOC_HDRS )

set( KRIGING_RCCS )

#---------------------------------------------------------------------------------------------------
# Wrapping Qt files
#qt4_wrap_ui( KRIGING_UIS_H  ${KRIGING_UIS} )
#qt4_wrap_cpp( KRIGING_MOC_SRCS  ${KRIGING_MOC_HDRS} )
#qt4_add_resources( KRIGING_RCC_SRCS ${KRIGING_RCCS} )

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
#  ${QT_INCLUDE_DIR}
  ..
  . 
)

#link_directories(
#  ${KRIGING_LIBRARIES}
#  ${MPI_LIBRARIES} 
#  . 
#)

#---------------------------------------------------------------------------------------------------  
# Creating libray kriging
  #${KRIGING_MOC_SRCS} ${KRIGING_UIS_H} ${KRIGING_RCC_SRCS} 
add_library( kriging STATIC ${KRIGING_SRCS} ${KRIGING_HDRS} )

set_target_properties( kriging PROPERTIES 
  LINKER_LANGUAGE CXX )

#  ${QT_LIBRARIES} 
target_link_libraries( kriging
  ${GA_LIBRARIES}
  ${GSL_LIBRARIES}
  ${MPI_LIBRARIES} )

#---------------------------------------------------------------------------------------------------
# Install
install( TARGETS kriging 
#  RUNTINE DESTINATION ${CMAKE_INSTALL_PREFIX}
  LIBRARY DESTINATION ${CMAKE_INSTALL_PREFIX}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_PREFIX}
  PUBLIC_HEADER DESTINATION ${CMAKE_INSTALL_PREFIX} )
