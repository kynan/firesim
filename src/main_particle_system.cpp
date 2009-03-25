//! \file main_particle_system.cpp
//! \date   Mar 24, 2009
//! \author Florian Rathgeber

#include <iostream>

#include "particles/ParticleSystem.h"

//! Main function

//! Runs a fire simulation according to the configuration file given.
//! \note The commandline for the executable is:
//!       <b>exename configFileName</b>
//! \param[in] configFileName \b string Path to the configuration file

int main ( int argc, char** argv ) {

  if ( argc < 2 || ( argc >= 2 &&
        ( strcmp( argv[1], "--help" ) == 0 || strcmp( argv[1], "-h" ) == 0 ) ) ) {
    std::cerr << "usage: " << argv[0] << " <configFileName>" << std::endl;
    exit( -1 );
  }

  particles::ParticleSystem myFire( argv[1] );
  myFire.run();

  return 0;
}
