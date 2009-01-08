/*
 * LBM.cpp
 *
 *  Created on: Jan 2, 2009
 *      Author: Florian Rathgeber
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <endian.h>

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
    : grid0_( new dfField( sizeX, sizeY, sizeZ ) ),
      grid1_( new dfField( sizeX, sizeY, sizeZ ) ),
      u_( sizeX, sizeY, sizeZ, 0. ),
      rho_( sizeX, sizeY, sizeZ, 1. ),
      boundaryCells_( boundaryCells ),
      velocityCells_( velocityCells ),
      velocities_( velocities ) {

  // initialize distribution functions with equilibrium
  for ( int z = 0; z < sizeZ; ++z )
    for ( int y = 0; y < sizeY; ++y )
      for ( int x = 0; x < sizeX; ++x )
        for ( int i = 0; i < Dim; ++i )
          (*grid0_)( x, y, z, i ) = w[i];
  for ( int z = 0; z < sizeZ; ++z )
    for ( int y = 0; y < sizeY; ++y )
      for ( int x = 0; x < sizeX; ++x )
        for ( int i = 0; i < Dim; ++i )
          (*grid1_)( x, y, z, i ) = w[i];

}

LBM::~LBM() {
  assert ( grid0_ != 0 && grid1_ != 0 );
  delete grid0_;
  delete grid1_;
}

// Set methods

void LBM::run( double omega,
               double cSmagorinsky,
               int maxSteps,
               int vtkStep,
               std::string vtkFileName ) {

  // relaxation parameter
  double tau = 1. / omega;
  // lattice viscosity
  double nu = ( 2. * tau - 1. ) * (1./6.);
  // squared Smagorinsky constant
  double cSqr = cSmagorinsky * cSmagorinsky;

  // loop over maxSteps time steps
  for ( int step = 0; step < maxSteps; ++step ) {

    std::cout << "Time step " << step << " of " << maxSteps << std::endl;

    // loop over all cells but the ghost layer
    for ( int z = 1; z < grid0_->getSizeZ() - 1; z++ ) {
      for ( int y = 1; y < grid0_->getSizeY() - 1; y++ ) {
        for ( int x = 1; x < grid0_->getSizeX() - 1; x++ ) {

          // Perform actual collision and streaming step
          collideStream( x, y, z, omega );
          // Perform collision and streaming step with Smagorinsky turbulence
          // correction
//          collideStreamSmagorinsky( x, y, z, nu, cSqr );

        } // x
      } // y
    } // z

    // Treat no-slip boundary conditions (walls)
    treatBoundary();

    // Treat velocity cells
    treatVelocities();

    // exchange grids for current and previous time step
    dfField *gridTmp = grid0_;
    grid0_ = grid1_;
    grid1_ = gridTmp;

    if ( step % vtkStep == 0 ) writeVtkFile( step, vtkFileName );

  } // step
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
  // DEBUG assertions
  assert ( rho > 0.8 && rho < 1.2 );
  assert ( fabs(ux) < 2. );
  assert ( fabs(uy) < 2. );
  assert ( fabs(uz) < 2. );
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
        + omega  * w[f] * ( fc +  3 * eiu + 4.5 * eiu * eiu);
  }
}

inline void LBM::collideStreamSmagorinsky( int x, int y, int z, double nu, double cSqr ) {

  // Calculate rho and u
  double rho = (*grid0_)( x, y, z, 0 ); // df in center
  double ux = 0.;
  double uy = 0.;
  double uz = 0.;
  // Loop over all velocity directions but center
  for ( int f = 1; f < Dim; ++f ) {
    double fi = (*grid0_)( x, y, z, f );
    rho += fi;
    ux += ex[f] * fi;
    uy += ey[f] * fi;
    uz += ez[f] * fi;
  }
  // DEBUG assertions
  assert ( rho > 0.8 && rho < 1.2 );
  assert ( fabs(ux) < 2. );
  assert ( fabs(uy) < 2. );
  assert ( fabs(uz) < 2. );
  rho_( x, y, z ) = rho;
  u_( x, y, z, 0 ) = ux;
  u_( x, y, z, 1 ) = uy;
  u_( x, y, z, 2 ) = uz;

  // Collision step: calculate equilibrium distribution values and
  // perform collision (weighting with current distribution values)
  // streaming step: stream distribution values to neighboring cells
  double fc = rho - 1.5 * ( ux * ux + uy * uy + uz * uz );
  double feq[19];

  // Calculate equilibrium distribution functions
  feq[0]  = (1./3.)  *   fc; // C
  feq[1]  = (1./18.) * ( fc + 3 *   uy        + 4.5 *   uy        *   uy ); // N
  feq[2]  = (1./18.) * ( fc + 3 *   ux        + 4.5 *   ux        *   ux ); // E
  feq[3]  = (1./18.) * ( fc - 3 *   uy        + 4.5 *   uy        *   uy ); // S
  feq[4]  = (1./18.) * ( fc - 3 *   ux        + 4.5 *   ux        *   ux ); // W
  feq[5]  = (1./18.) * ( fc + 3 *   uz        + 4.5 *   uz        *   uz ); // T
  feq[6]  = (1./18.) * ( fc - 3 *   uz        + 4.5 *   uz        *   uz ); // B
  feq[7]  = (1./36.) * ( fc + 3 * ( ux + uy ) + 4.5 * ( ux + uy ) * ( ux + uy ) ); // NE
  feq[8]  = (1./36.) * ( fc + 3 * ( ux - uy ) + 4.5 * ( ux - uy ) * ( ux - uy ) ); // SE
  feq[9]  = (1./36.) * ( fc - 3 * ( ux + uy ) + 4.5 * ( ux + uy ) * ( ux + uy ) ); // SW
  feq[10] = (1./36.) * ( fc - 3 * ( ux - uy ) + 4.5 * ( ux - uy ) * ( ux - uy ) ); // NW
  feq[11] = (1./36.) * ( fc + 3 * ( uy + uz ) + 4.5 * ( uy + uz ) * ( uy + uz ) ); // TN
  feq[12] = (1./36.) * ( fc + 3 * ( ux + uz ) + 4.5 * ( ux + uz ) * ( ux + uz ) ); // TE
  feq[13] = (1./36.) * ( fc - 3 * ( uy - uz ) + 4.5 * ( uy - uz ) * ( uy - uz ) ); // TS
  feq[14] = (1./36.) * ( fc - 3 * ( ux - uz ) + 4.5 * ( ux - uz ) * ( ux - uz ) ); // TW
  feq[15] = (1./36.) * ( fc + 3 * ( uy - uz ) + 4.5 * ( uy - uz ) * ( uy - uz ) ); // BN
  feq[16] = (1./36.) * ( fc + 3 * ( ux - uz ) + 4.5 * ( ux - uz ) * ( ux - uz ) ); // BE
  feq[17] = (1./36.) * ( fc - 3 * ( uy + uz ) + 4.5 * ( uy + uz ) * ( uy + uz ) ); // BS
  feq[18] = (1./36.) * ( fc - 3 * ( ux + uz ) + 4.5 * ( ux + uz ) * ( ux + uz ) ); // BW

  // Calculate non-equilibrium stress tensor
  double qo = 0.;
  for ( int i = 0; i < 3; ++i ) {
    double qadd = 0.;
    for ( int f = 1; f < 19; ++f ) {
      qadd += ep[i][f] * ( (*grid0_)( x, y, z, f ) - feq[f] );
    }
    qo += qadd * qadd;
  }
  qo *= 2.;
  for ( int i = 4; i < 6; ++i ) {
    double qadd = 0.;
    for ( int f = 7; f < 19; ++f ) {
      qadd += ep[i][f] * ( (*grid0_)( x, y, z, f ) - feq[f] );
    }
    qo += qadd * qadd;
  }
  qo = sqrt( qo );

  // Calculate local stress tensor
  double s = ( sqrt( nu * nu + 18. * cSqr * qo ) - nu ) / ( 6. * cSqr);
  // Calculate turbulence modified inverse lattice viscosity
  double omega = 1. / ( 3. * ( nu + cSqr * s ) + .5 );
  double omegai = 1. - omega;

  // Loop over all velocity directions and stream collided distribution value
  // to neighboring cells
  for ( int f = 0; f < Dim; ++f ) {
    (*grid1_)( x + ex[f], y + ey[f], z + ez[f], f )
      =   omegai * (*grid0_)( x, y, z, f ) + omega * feq[f];
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
      (*grid1_)( x - ex[i], y - ey[i], z - ez[i], finv[i] ) = (*grid1_)( x, y, z, i );
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
    // Set velocity of this cell
    u_( x, y, z, 0 ) = ux;
    u_( x, y, z, 1 ) = uy;
    u_( x, y, z, 2 ) = uz;

    // Go over all distribution values, stream to inverse distribution value of
    // adjacent cell in inverse direction (bounce back) and modify by velocity
    // of moving wall
    for ( int i = 1; i < Dim; ++i ) {
      int op = finv[i];
      (*grid1_)( x - ex[i], y - ey[i], z - ez[i], op )
        = (*grid1_)( x, y, z, i )
          + rho * w[i] * ( ex[op] * ux + ey[op] * uy + ez[op] * uz );
    }
  }
}

void LBM::writeVtkFile( int timestep, std::string vtkFileName ) {

  // Open file for writing
  std::ostringstream oss;
  oss << vtkFileName << "." << timestep << ".vtk";
  std::cout << "Writing file '" << oss.str() << "' for time step " << timestep << std::endl;
  std::ofstream vtkFile( oss.str().c_str(), std::ios::binary | std::ios::out );

  // Get size of domain without ghost layers
  int sizeX = grid0_->getSizeX() - 2;
  int sizeY = grid0_->getSizeY() - 2;
  int sizeZ = grid0_->getSizeZ() - 2;

  // Write file header
  vtkFile << "# vtk DataFile Version 2.0\n";
  vtkFile << "VTK output file for time step " << timestep << "\n\n";
  vtkFile << "BINARY\n\n";
  vtkFile << "DATASET STRUCTURED_POINTS\n";
  vtkFile << "DIMENSIONS " << sizeX << " " << sizeY << " " << sizeZ << "\n";
  vtkFile << "ORIGIN 0.0 0.0 0.0\n";
  vtkFile << "SPACING 1.0 1.0 1.0\n\n";
  vtkFile << "POINT_DATA " << sizeX * sizeY * sizeZ << "\n\n";

  // Write density field
  vtkFile << "SCALARS density double\n";
  vtkFile << "LOOKUP_TABLE default\n";
  for ( int z = 1; z <= sizeZ; ++z ) {
    for ( int y = 1; y <= sizeY; ++y ) {
      for ( int x = 1; x <= sizeX; ++x ) {
        // evil hack because vtk requires binary data to be in big endian
        uint64_t dump = htobe64( *reinterpret_cast<uint64_t *>( &rho_( x, y, z ) ) );
        vtkFile.write( reinterpret_cast<char *>( &dump ), sizeof(double) );
      }
    }
  }

  // Write velocity vector field
  vtkFile << "VECTORS velocity double\n";
  for ( int z = 1; z <= sizeZ; ++z ) {
    for ( int y = 1; y <= sizeY; ++y ) {
      for ( int x = 1; x <= sizeX; ++x ) {
        // evil hack because vtk requires binary data to be in big endian
        uint64_t dump0 = htobe64( *reinterpret_cast<uint64_t *>( &u_( x, y, z, 0 ) ) );
        uint64_t dump1 = htobe64( *reinterpret_cast<uint64_t *>( &u_( x, y, z, 1 ) ) );
        uint64_t dump2 = htobe64( *reinterpret_cast<uint64_t *>( &u_( x, y, z, 2 ) ) );
        vtkFile.write( reinterpret_cast<char *>( &dump0 ), sizeof(double) );
        vtkFile.write( reinterpret_cast<char *>( &dump1 ), sizeof(double) );
        vtkFile.write( reinterpret_cast<char *>( &dump2 ), sizeof(double) );
      }
    }
  }

}

} // namespace lbm
