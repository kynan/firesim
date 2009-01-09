/*
 * main.cpp
 *
 *  Created on: Jan 2, 2009
 *      Author: Florian Rathgeber
 */

#include <cstdlib>
#include <iostream>

#include "LBM.h"

int main ( int argc, char** argv ) {

  if ( argc < 9 ) {
    std::cerr << "usage: " << argv[0];
    std::cerr << " sizeX sizeY sizeZ omega cSmago timesteps vtkStep vtkFile" << std::endl;
    exit( EXIT_FAILURE );
  }
  int sizeX = atoi( argv[1] );
  int sizeY = atoi( argv[2] );
  int sizeZ = atoi( argv[3] );
  double omega = atof( argv[4] );
  double cSmago = atof( argv[5] );
  int timesteps = atoi( argv[6] );
  int vtkStep = atoi( argv[7] );
  char* vtkFile = argv[8];

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
      bound.push_back( Vec3<int>( x, sizeY - 2, z ) );
    }
  // Left and right side
  for ( int z = 2; z < sizeZ - 2; ++z )
    for ( int y = 2; y < sizeY - 2; ++y ) {
      bound.push_back( Vec3<int>( 1, y, z ) );
      bound.push_back( Vec3<int>( sizeX - 2, y, z ) );
    }

  std::vector< Vec3<int> > lid;
  // Create lid cells
  for ( int y = 1; y < sizeY - 1; ++y )
      for ( int x = 1; x < sizeX - 1; ++x )
        lid.push_back( Vec3<int>( x, y, sizeZ - 2 ) );

  std::vector< Vec3<double> > vel;
  // Create lid velocities
  for ( int i = 0; i < (sizeX - 2) * (sizeY - 2); ++i )
    vel.push_back( Vec3<double>( 0.08, 0., 0. ) );


  lbm::LBM myLBM( sizeX, sizeY, sizeZ, bound, lid, vel );
  myLBM.run( omega, cSmago, timesteps, vtkStep, vtkFile );

  return EXIT_SUCCESS;
}
