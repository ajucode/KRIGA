#___________________________________________________________________________________________________
#
# Try to find the Global Arrays library
# See
# 
# author: Pedro Guarderas
# date: 21-03-2013
# Provides:
#	GA_INCLUDE_DIR
#	GA_LIBRARY
#

set( GA_FOUND OFF )


#if( UNIX )
find_program( GA_CONFIG_EXECUTABLE ga-config
  $ENV{GA_DIR}/bin
  ${GA_DIR}/bin
  /usr/bin/
  /usr/local/bin )

if( GA_CONFIG_EXECUTABLE ) 
  set( GA_FOUND ON )
  execute_process(
    COMMAND sh "${GA_CONFIG_EXECUTABLE}" --cflags --cppflags --cflags
    OUTPUT_VARIABLE GA_CFLAGS
    RESULT_VARIABLE RET
    ERROR_QUIET )
  
  if( RET EQUAL 0 )

    string( STRIP "${GA_CFLAGS}" GA_CFLAGS )
    separate_arguments( GA_CFLAGS )
    
    # parse definitions from cflags; drop -D* from CFLAGS
    string( REGEX MATCHALL "-D[^;]+"
      GA_DEFINITIONS  "${GA_CFLAGS}" )
    string( REGEX REPLACE "-D[^;]+;" ""
      GA_CFLAGS "${GA_CFLAGS}" )
    
    # parse include dirs from cflags; drop -I prefix
    string( REGEX MATCHALL "-I[^;]+"
      GA_INCLUDE_DIRS "${GA_CFLAGS}" )
    string( REPLACE "-I" ""
      GA_INCLUDE_DIRS "${GA_INCLUDE_DIRS}")
    string( REGEX REPLACE "-I[^;]+;" ""
      GA_CFLAGS "${GA_CFLAGS}")

  else( RET EQUAL 0 )
    set( GA_FOUND FALSE )
  endif( RET EQUAL 0 )
  
  execute_process(
    COMMAND sh "${GA_CONFIG_EXECUTABLE}" --libs --ldflags -fflags
    OUTPUT_VARIABLE GA_LIBRARIES
    RESULT_VARIABLE RET
    ERROR_QUIET )

  if( RET EQUAL 0 )
    set( GA_LIBRARIES "-lga++ ${GA_LIBRARIES}" )
    string(STRIP "${GA_LIBRARIES}" GA_LIBRARIES )
    separate_arguments( GA_LIBRARIES )
    
    # extract linkdirs (-L) for rpath (i.e., LINK_DIRECTORIES)
    string( REGEX MATCHALL "-L[^;]+"
      GA_LIBRARY_DIRS "${GA_LIBRARIES}" )
    string( REPLACE "-L" ""
      GA_LIBRARY_DIRS "${GA_LIBRARY_DIRS}" )

  else( RET EQUAL 0 )
    set( GA_FOUND FALSE )
  endif( RET EQUAL 0 )

  mark_as_advanced( GA_CFLAGS )
  
else( GA_CONFIG_EXECUTABLE )
  message( STATUS "FindGA: ga-config not found.")
endif( GA_CONFIG_EXECUTABLE )

#endif( UNIX )

mark_as_advanced(
  GA_CONFIG_EXECUTABLE
  GA_INCLUDE
  GA_LIBRARY
  GA_CBLAS_LIBRARY )