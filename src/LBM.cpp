/*
 * LBM.cpp
 *
 *  Created on: Jan 2, 2009
 *      Author: Florian Rathgeber
 */

#include "LBM.h"

namespace lbm {

// Constructors and destructors

LBM::LBM() {}

LBM::LBM( int sizeX,
          int sizeY,
          int sizeZ,
          std::vector< Vec3<int> > &boundaryCells,
          std::vector< Vec3<int> > &velocityCells,
          std::vector< Vec3<double> > &velocities )
    : grid0_( new dfField( sizeX, sizeY, sizeZ, 0. ) ),
      grid1_( new dfField( sizeX, sizeY, sizeZ, 0. ) ),
      u_( sizeX, sizeY, sizeZ, 0. ),
      rho_( sizeX, sizeY, sizeZ, 1. ),
      boundaryCells_( boundaryCells ) {}

LBM::~LBM() {
  assert ( grid0_ != 0 && grid1_ != 0 );
  delete grid0_;
  delete grid1_;
}

// Set methods

void LBM::run( double omega, int maxSteps ) {

  // loop over maxSteps time steps
  for ( int step = 0; step < maxSteps; ++step ) {
    // loop over all cells but the ghost layer
    for ( int z = 1; z < grid0_->getSizeZ() - 1; z++ ) {
      for ( int y = 1; y < grid0_->getSizeY() - 1; y++ ) {
        for ( int x = 1; x < grid0_->getSizeX() - 1; x++ ) {

          // Perform actual collision and streaming step
          collideStream( x, y, z, omega );

        } // x
      } // y
    } // z
  } // step

  // Treat no-slip boundary conditions (walls)
  treatBoundary();

  // Treat velocity cells
  treatVelocities();

  // exchange grids for current and previous time step
  dfField *gridTmp = grid0_;
  grid0_ = grid1_;
  grid1_ = gridTmp;
}

// Internal helper functions

inline void LBM::collideStream( int x, int y, int z, double omega ) {

  // calculate rho and u
  double rho = (*grid0_)( x, y, z, 0 ); // df in center
  double ux = 0.;
  double uy = 0.;
  double uz = 0.;
  // loop over all velocity directions but center
  for ( int f = 1; f < Dim; ++f ) {
    double fi = (*grid0_)( x, y, z, f );
    rho += fi;
    ux += ex[f] * fi;
    uy += ey[f] * fi;
    uz += ez[f] * fi;
  }
  rho_( x, y, z ) = rho;
  u_( x, y, z, 0 ) = ux;
  u_( x, y, z, 1 ) = uy;
  u_( x, y, z, 2 ) = uz;

  // collision step: calculate equilibrium distribution values and
  // perform collision (weighting with current distribution values)
  // streaming step: stream distribution values to neighboring cells
  double fc = rho - 1.5 * ( ux * ux + uy * uy + uz * uz );
  double omegai = 1 - omega;
  // treat center value specially
  (*grid1_)( x, y, z, 0 ) = omegai * (*grid0_)( x, y, z, 0 )
                        + omega *  w[0] * fc;
  // loop over all velocity directions but center
  for ( int f = 1; f < Dim; ++f ) {
    double eiu = ex[f] * ux + ey[f] * uy + ez[f] * uz;
    (*grid1_)( x + ex[f], y + ey[f], z + ez[f], f )
      =   omegai * (*grid0_)( x, y, z, f )
        + omega  * ( fc +  3 * eiu + 4.5 * eiu * eiu);
  }
}

inline void LBM::treatBoundary() {

  // Iterate over all no-slip boundary cells
  std::vector< Vec3<int> >::iterator iter;
  for( iter = boundaryCells_.begin(); iter != boundaryCells_.end(); iter++ ) {

    // Fetch coordinates of current boundary cell
    int x = (*iter)[0];
    int y = (*iter)[1];
    int z = (*iter)[2];

    // Go over all distribution values and stream to inverse distribution
    // value of adjacent cell in inverse direction (bounce back)
    for ( int i = 1; i < Dim; ++i ) {
      (*grid1_)( x + erx[i], y + erx[i], z + erx[i], finv[i] ) = (*grid1_)( x, y, z, i );
    }
  }
}

inline void LBM::treatVelocities() {

  // Iterate over all velocity boundary cells
  for( int i = 0; i < velocityCells_.size(); ++i ) {

    // Fetch coordinates of current boundary cell
    int x = velocityCells_[i][0];
    int y = velocityCells_[i][1];
    int z = velocityCells_[i][2];
    // Fetch velocity of moving wall
    double ux = velocities_[i][0];
    double uy = velocities_[i][1];
    double uz = velocities_[i][2];
    // Fetch density of current cell
    double rho = 6 * rho_( x, y, z );

    // Go over all distribution values, stream to inverse distribution value of
    // adjacent cell in inverse direction (bounce back) and modify by velocity
    // of moving wall
    for ( int i = 1; i < Dim; ++i ) {
      (*grid1_)( x + erx[i], y + erx[i], z + erx[i], finv[i] )
        = (*grid1_)( x, y, z, i )
          - rho * w[i] * ( ex[i] * ux + ey[i] * uy + ez[i] * uz );
    }
  }
}

} // namespace lbm
