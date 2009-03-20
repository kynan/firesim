# Compiler flags for custom build type
SET( CMAKE_CXX_FLAGS_DEBUGNSMAGO "${CMAKE_CXX_FLAGS_DEBUG} -DNSMAGO" )
SET( CMAKE_C_FLAGS_DEBUGNSMAGO "${CMAKE_C_FLAGS_DEBUG} -DNSMAGO" )
SET( CMAKE_CXX_FLAGS_RELEASENSMAGO "${CMAKE_CXX_FLAGS_RELEASE} -DNSMAGO" )
SET( CMAKE_C_FLAGS_RELEASENSMAGO "${CMAKE_C_FLAGS_RELEASE} -DNSMAGO" )

# Module to automatically create doc target
INCLUDE( TargetDoc.cmake OPTIONAL)

ADD_SUBDIRECTORY( confparser )

MACRO(ADD_FILE_DEPENDANCY file)
   SET(${file}_deps ${${file}_deps} ${ARGN})
   SET_SOURCE_FILES_PROPERTIES( ${file}
      PROPERTIES OBJECT_DEPENDS "${${file}_deps}" )
ENDMACRO(ADD_FILE_DEPENDANCY)

#ADD_FILE_DEPENDANCY( main.cpp lbm/LBM.h lbm/LBM_def.h particles/ParticleSystem.h particles/Particle.h particles/Emitter.h )
ADD_EXECUTABLE( lbm_reference main.cpp )
TARGET_LINK_LIBRARIES( lbm_reference confparser Irrlicht GL Xxf86vm Xext X11 )