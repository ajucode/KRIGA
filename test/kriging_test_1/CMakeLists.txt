#------------------------------------------------------------------------------
enable_testing()

include_directories(
  ${GA_INCLUDE_DIRS}
  ../../src
  .
)

set( KRIG_T1_SRCS
  main.cpp
)

add_executable( krig_test_1
  ${KRIG_T1_SRCS}
  ${MPIEXEC} )

#set_target_properties( krig_test_1 PROPERTIES )

target_link_libraries( krig_test_1
  kriging
  ${GA_LIBRARIES}
  ${GSL_LIBRARIES}
)


add_test( krig_test_1 ${CMAKE_CURRENT_BINARY_DIR}/krig_test_1 )