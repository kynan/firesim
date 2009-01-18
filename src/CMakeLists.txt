# Compiler flags for custom build type
SET( CMAKE_CXX_FLAGS_DEBUGNSMAGO "${CMAKE_CXX_FLAGS_DEBUG} -DNSMAGO" )
SET( CMAKE_C_FLAGS_DEBUGNSMAGO "${CMAKE_C_FLAGS_DEBUG} -DNSMAGO" )
SET( CMAKE_CXX_FLAGS_RELEASENSMAGO "${CMAKE_CXX_FLAGS_RELEASE} -DNSMAGO" )
SET( CMAKE_C_FLAGS_RELEASENSMAGO "${CMAKE_C_FLAGS_RELEASE} -DNSMAGO" )

# Module to automatically create doc target
INCLUDE( TargetDoc OPTIONAL)

ADD_EXECUTABLE( lbm_reference main.cpp confparser/ConfParser.cpp confparser/ConfBlock.cpp )
