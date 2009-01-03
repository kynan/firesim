/*
 * LBM.cpp
 *
 *  Created on: Jan 2, 2009
 *      Author: Florian Rathgeber
 */

#include "LBM.h"

namespace lbm {

// Constructors and destructors

LBM::LBM() {
}

LBM::LBM( int sizeX, int sizeY, int sizeZ )
    : grid0_( new dfField( sizeX, sizeY, sizeZ, 0. ) ),
      grid1_( new dfField( sizeX, sizeY, sizeZ, 0. ) ),
      u_( sizeX, sizeY, sizeZ, 0. ),
      rho_( sizeX, sizeY, sizeZ, 1. ) {
}

LBM::~LBM() {
  assert ( grid0_ != 0 && grid1_ != 0 );
  delete grid0_;
  delete grid1_;
}

// Set methods

void LBM::run( double omega, int maxSteps ) {

  // loop over maxSteps time steps
  for ( int step = 0; step < maxSteps; ++step ) {
    // loop over all non-boundary cells
    for ( int z = 1; z < grid0_->getSizeZ() - 1; z++ ) {
      for ( int y = 1; y < grid0_->getSizeY() - 1; y++ ) {
        for ( int x = 1; x < grid0_->getSizeX() - 1; x++ ) {

          // calculate rho and u
          double rho = (*grid0_)( x, y, z, 0 ); // df in center
          double ux = 0.;
          double uy = 0.;
          double uz = 0.;
          // loop over all velocity directions but center
          for ( int f = 1; f < dim; ++f ) {
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
          for ( int f = 1; f < dim; ++f ) {
            double eiu = ex[f] * ux + ey[f] * uy + ez[f] * uz;
            (*grid1_)( x + ex[f], y + ey[f], z + ez[f], f )
              =   omegai * (*grid0_)( x, y, z, f )
                + omega  * ( fc +  3 * eiu + 4.5 * eiu * eiu);
          }
        } // x
      } // y
    } // z
  } // step

  // exchange grids for current and previous time step
  dfField *gridTmp = grid0_;
  grid0_ = grid1_;
  grid1_ = gridTmp;
}

}
