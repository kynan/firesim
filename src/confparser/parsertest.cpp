//! \file parsertest.cpp
//! Test setup for the ConfParser

//! \date   Jan 19, 2009
//! \author Florian Rathgeber

#include <cstdlib>
#include <iostream>
#include <string>

#include "ConfParser.h"
#include "ConfBlock.h"

using namespace std;
using namespace confparser;

int main( int argc, char** argv ) {

  if ( argc < 2 ) {
    cerr << "Usage: " << argv[0] << " <config_file> [<config_file> ...]" << endl;
    return EXIT_FAILURE;
  }

  ConfParser p;

  for ( int i = 1; i < argc; ++i ) {
    try {
      ConfBlock* b = p.parse( argv[i] );
      b->writeConfigFile( string( argv[i] ) + ".out" );
    } catch ( BadSyntax e ) {
      cerr << e.what() << endl;
    }
  }

  return EXIT_SUCCESS;
}
