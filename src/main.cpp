/*
 * main.cpp
 *
 *  Created on: Jan 2, 2009
 *      Author: Florian Rathgeber
 */

#include <cstdlib>

#include "LBM.h"

int main ( int argc, char** argv ) {

  lbm::LBM myLBM( 40, 40, 40 );
  myLBM.run( 1.9, 100 );

  return EXIT_SUCCESS;
}
