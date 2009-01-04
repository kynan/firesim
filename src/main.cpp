/*
 * main.cpp
 *
 *  Created on: Jan 2, 2009
 *      Author: Florian Rathgeber
 */

#include <cstdlib>

#include "LBM.h"

int main ( int argc, char** argv ) {

  int sizeX = 70;
  int sizeY = 80;
  int sizeZ = 90;

  std::vector< Vec3<int> > bound;
  // Create boundary cells
  // Bottom
  for ( int y = 1; y < sizeY - 1; ++y )
    for ( int x = 1; x < sizeX - 1; ++x )
      bound.push_back( Vec3<int>( x, y, 1 ) );
  // Front and back side
  for ( int z = 2; z < sizeZ - 2; ++z )
    for ( int x = 1; x < sizeX - 1; ++ x ) {
      bound.push_back( Vec3<int>( x, 1, z ) );
      bound.push_back( Vec3<int>( x, 38, z ) );
    }
  // Left and right side
  for ( int z = 2; z < sizeZ - 2; ++z )
    for ( int y = 2; y < sizeY - 2; ++y ) {
      bound.push_back( Vec3<int>( 1, y, z ) );
      bound.push_back( Vec3<int>( 38, y, z ) );
    }

  std::vector< Vec3<int> > lid;
  // Create lid cells
  for ( int y = 1; y < sizeY - 1; ++y )
      for ( int x = 1; x < sizeX - 1; ++x )
        bound.push_back( Vec3<int>( x, y, 38 ) );

  std::vector< Vec3<double> > vel;
  // Create lid velocities
  for ( int i = 0; i < (sizeX - 2) * (sizeY - 2); ++i )
    vel.push_back( Vec3<double>( 0.08, 0., 0. ) );


  lbm::LBM myLBM( sizeX, sizeY, sizeZ, bound, lid, vel );
  myLBM.run( 1.9, 100 );

  return EXIT_SUCCESS;
}
