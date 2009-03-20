//! \file main.cpp
//! \date   Jan 2, 2009
//! \author Florian Rathgeber

#include <iostream>

#include "particles/ParticleSystem.h"

#ifndef Real
  #define Real float
#endif

//! Main function

//! Initializes a lid driven cavity setup with no slip boundaries at all sides
//! and a fixed velocity of 0.08 in x-direction at the lid (top w.r.t. to z)
//! \note The commandline for the executable is:
//!       <b>exename configFileName</b>
//! \param[in] configFileName \b string Path to the configuration file

int main ( int argc, char** argv ) {

  if ( argc < 2 || ( argc >= 2 &&
        ( strcmp( argv[1], "--help" ) == 0 || strcmp( argv[1], "-h" ) == 0 ) ) ) {
    std::cerr << "usage: " << argv[0] << " <configFileName>" << std::endl;
    exit( -1 );
  }

  particles::ParticleSystem<Real> myFire( argv[1] );
  myFire.run();

  return 0;
}
