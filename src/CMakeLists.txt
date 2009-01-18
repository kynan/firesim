# Compiler flags for custom build type
SET( CMAKE_CXX_FLAGS_DEBUGNSMAGO "${CMAKE_CXX_FLAGS_DEBUG} -DNSMAGO" )
SET( CMAKE_C_FLAGS_DEBUGNSMAGO "${CMAKE_C_FLAGS_DEBUG} -DNSMAGO" )
SET( CMAKE_CXX_FLAGS_RELEASENSMAGO "${CMAKE_CXX_FLAGS_RELEASE} -DNSMAGO" )
SET( CMAKE_C_FLAGS_RELEASENSMAGO "${CMAKE_C_FLAGS_RELEASE} -DNSMAGO" )

# Module to automatically create doc target
INCLUDE( TargetDoc OPTIONAL)

SET( BOOST_ROOT "" CACHE PATH "Path to the boost installation" )
SET( BOOST_INCLUDEDIR "" CACHE PATH "Path to the boost headers" )
SET( BOOST_LIBRARYDIR "" CACHE PATH "Path to the boost library binaries" )

# Check for presence of Boost, especially regex component
FIND_PACKAGE( Boost 1.18.0 REQUIRED COMPONENTS regex )
IF( NOT Boost_FOUND )
    MESSAGE( FATAL_ERROR "Boost library not found!" )
ENDIF( NOT Boost_FOUND )

ADD_EXECUTABLE( lbm_reference main.cpp confparser/ConfParser.cpp confparser/ConfBlock.cpp )
TARGET_LINK_LIBRARIES( lbm_reference ${Boost_LIBRARIES} )
