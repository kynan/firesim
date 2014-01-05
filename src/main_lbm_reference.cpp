//! \file main_lbm_reference.cpp
//! \date   Mar 24, 2009
//! \author Florian Rathgeber

#ifndef RealType
#define RealType float
#endif

#include <iostream>

#include "lbm/LBM.h"
#include "lbm/LBM_def.h"

//! Main function

//! Runs a Lattice Boltzmann simulation according to the configuration file given.
//! \note The commandline for the executable is:
//!       <b>exename configFileName</b>
//! \param[in] configFileName \b string Path to the configuration file

int main ( int argc, char** argv ) {

  if ( argc < 2 || ( argc >= 2 &&
        ( strcmp( argv[1], "--help" ) == 0 || strcmp( argv[1], "-h" ) == 0 ) ) ) {
    std::cerr << "usage: " << argv[0] << " <configFileName>" << std::endl;
    exit( -1 );
  }

  lbm::LBM<RealType> myLBM( argv[1] );
  myLBM.run();

  return 0;
}
