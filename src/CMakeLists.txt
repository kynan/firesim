# Compiler flags for custom build type
SET( CMAKE_CXX_FLAGS_DEBUGNSMAGO "${CMAKE_CXX_FLAGS_DEBUG} -DNSMAGO" )
SET( CMAKE_C_FLAGS_DEBUGNSMAGO "${CMAKE_C_FLAGS_DEBUG} -DNSMAGO" )
SET( CMAKE_CXX_FLAGS_RELEASENSMAGO "${CMAKE_CXX_FLAGS_RELEASE} -DNSMAGO" )
SET( CMAKE_C_FLAGS_RELEASENSMAGO "${CMAKE_C_FLAGS_RELEASE} -DNSMAGO" )

SET( PLATFORM "cip" CACHE STRING "Platform to build for" )
#SET( BIN_ROOT_DIR "${CMAKE_SOURCE_DIR}" CACHE PATH "Folder where to place binaries" )

# Module to automatically create doc target
INCLUDE( TargetDoc.cmake OPTIONAL)

ADD_SUBDIRECTORY( confparser )

MACRO(ADD_FILE_DEPENDANCY file)
   SET(${file}_deps ${${file}_deps} ${ARGN})
   SET_SOURCE_FILES_PROPERTIES( ${file}
      PROPERTIES OBJECT_DEPENDS "${${file}_deps}" )
ENDMACRO(ADD_FILE_DEPENDANCY)

ADD_EXECUTABLE( lbm_reference main_lbm_reference.cpp )
SET_TARGET_PROPERTIES( lbm_reference PROPERTIES COMPILE_FLAGS "-DLBM_ONLY" )
TARGET_LINK_LIBRARIES( lbm_reference confparser )
ADD_EXECUTABLE( firesim main_particle_system.cpp
                        particles/ParticleSystem.cpp
                        particles/Particle.cpp
                        particles/Emitter.cpp )
IF( ${PLATFORM} STREQUAL "cluster" )
  SET_TARGET_PROPERTIES( firesim PROPERTIES LINK_FLAGS "-L/central/irrlicht/lib/Linux" )
ENDIF( ${PLATFORM} STREQUAL "cluster" )
TARGET_LINK_LIBRARIES( firesim confparser Irrlicht GL Xxf86vm Xext X11 )
SET_TARGET_PROPERTIES( lbm_reference firesim PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_SOURCE_DIR}/bin/${PLATFORM}/${CMAKE_BUILD_TYPE}" )
