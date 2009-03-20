
CMAKE_MINIMUM_REQUIRED( VERSION 2.4 )
CMAKE_POLICY(SET CMP0003 NEW)

PROJECT( lbm_reference )

SET( CMAKE_C_FLAGS_DEBUG "-Wall -Wextra -Winline -Wshadow -ansi -pedantic -fmessage-length=0 -g -DDEBUG" )
SET( CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG -fmessage-length=0" )
SET( CMAKE_CXX_FLAGS_DEBUG "-Wall -Wextra -Winline -Wshadow -ansi -pedantic -fmessage-length=0 -g -DDEBUG" )
SET( CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG -fmessage-length=0" )

SET( CMAKE_VERBOSE_MAKEFILE ON )

ADD_SUBDIRECTORY( src )
