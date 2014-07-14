#------------------------------------------------------------------------------
if( KRIGING_ENABLE_TESTS )
#  enable_testing()
  add_subdirectory( 
#    matrix_test_1
    kriging_test_1
    #sorting_test_1 
    #twontree_test_1  
)
endif( KRIGING_ENABLE_TESTS )