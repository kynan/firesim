SET( Boost_ADDITIONAL_VERSIONS "1.37" "1.37.0" "1.38" "1.38.0" CACHE STRING "Additional Boost library versions to search for" )
SET( BOOST_ROOT "" CACHE PATH "Path to the boost installation" )
SET( BOOST_INCLUDEDIR "" CACHE PATH "Path to the boost headers" )
SET( BOOST_LIBRARYDIR "" CACHE PATH "Path to the boost library binaries" )

# Check for presence of Boost, especially regex component
FIND_PACKAGE( Boost 1.37.0 REQUIRED COMPONENTS regex )
IF( NOT Boost_FOUND )
    MESSAGE( FATAL_ERROR "Boost library not found!" )
ENDIF( NOT Boost_FOUND )

ADD_LIBRARY( confparser ConfBlock.cpp ConfParser.cpp )
TARGET_LINK_LIBRARIES( confparser ${Boost_LIBRARIES} )

ADD_EXECUTABLE( parsertest parsertest.cpp )
TARGET_LINK_LIBRARIES( parsertest confparser )