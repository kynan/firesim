//! \file LBM_def.h
//! Implementation of the LBM class

//! \date   Jan 16, 2009
//! \author Florian Rathgeber

#ifndef LBM_DEF_H_
#define LBM_DEF_H_

namespace lbm {

// ============================ //
// Constructors and destructors //
// ============================ //

template<typename T>
LBM<T>::LBM( std::string configFileName ) {
  setup( configFileName );
}

template<typename T>
LBM<T>::~LBM() {
  assert ( grid0_ != 0 && grid1_ != 0 );
  delete grid0_;
  delete grid1_;
}

template<typename T>
void LBM<T>::run() {

#ifndef NSMAGO
  // relaxation parameter
  T tau = 1. / omega_;
  // lattice viscosity
  T nu = (2. * tau - 1.) * (1. / 6.);
  // squared Smagorinsky constant
  T cSqr = cSmagorinsky_ * cSmagorinsky_;
#endif

  T totTime = 0.;

  // Get effective domain size (without ghost layer)
  int sizeX = grid0_->getSizeX() - 2;
  int sizeY = grid0_->getSizeY() - 2;
  int sizeZ = grid0_->getSizeZ() - 2;

  int numCells = sizeX * sizeY * sizeZ;

#ifdef NSMAGO
  std::cout << "Starting LBM without Smagorinsky turbulence correction" << std::endl;
#else
  std::cout << "Starting LBM with Smagorinsky turbulence correction"
      << std::endl;
#endif

  // loop over maxSteps time steps
  for (int step = 0; step < maxSteps_; ++step) {

    struct timeval start, end;
    T stepTime = 0.;

    gettimeofday(&start, NULL);

    // loop over all cells but the ghost layer
    for (int z = 1; z <= sizeZ; z++) {
      for (int y = 1; y <= sizeY; y++) {
        for (int x = 1; x <= sizeX; x++) {

#ifdef NSMAGO
          // Perform actual collision and streaming step
          collideStream( x, y, z, omega_ );
#else
          // Perform collision and streaming step with Smagorinsky turbulence
          // correction
          collideStreamSmagorinsky(x, y, z, nu, cSqr);
#endif

        } // x
      } // y
    } // z

    gettimeofday(&end, NULL);
    T scTime = getTime(start, end);

    // Treat no-slip boundary conditions (walls)
    treatBoundary();

    gettimeofday(&start, NULL);
    T bTime = getTime(end, start);

    // Treat velocity cells
    treatVelocities();

    gettimeofday(&end, NULL);
    T vTime = getTime(start, end);

    // exchange grids for current and previous time step
    Grid<T, Dim> *gridTmp = grid0_;
    grid0_ = grid1_;
    grid1_ = gridTmp;

    stepTime = scTime + bTime + vTime;
    totTime += stepTime;

    std::cout << "Time step " << step << " of " << maxSteps_ << " took ";
    std::cout << scTime << " + " << bTime << " + " << vTime << " = ";
    std::cout << stepTime << "secs -> " << numCells / (stepTime * 1000000);
    std::cout << " MLUP/s" << std::endl;

    if ( vtkStep_ != 0 && step % vtkStep_ == 0 )
      writeVtkFile(step);

  } // step

  std::cout << "LBM finished! Processed " << maxSteps_ << " timesteps in ";
  std::cout << totTime << " secs!" << std::endl;
  std::cout << "Average speed of " << (maxSteps_ * numCells) / (totTime
      * 1000000);
  std::cout << " MLUP/s" << std::endl;
}

// ========================= //
// Internal helper functions //
// ========================= //

template<typename T>
void LBM<T>::setup( std::string configFileName ) {
  try {

    ConfParser p;
    ConfBlock* base = p.parse( configFileName );
    std::cout << "Parsed configuration file " << configFileName << std::endl;

    // Read the parameters from the config file

    ConfBlockIterator bit = base->findRec( "domain" );
    if ( bit == base->end( "domain" ) ) {
      throw "No domain size given.";
    }
    ConfBlock paramBlock = *bit;
    int sizeX = paramBlock.getParam<int>( "sizeX" );
    int sizeY = paramBlock.getParam<int>( "sizeY" );
    int sizeZ = paramBlock.getParam<int>( "sizeZ" );
    std::cout << "Read domain specification:" << std::endl;
    std::cout << "sizeX : " << sizeX << std::endl;
    std::cout << "sizeY : " << sizeY << std::endl;
    std::cout << "sizeZ : " << sizeZ << std::endl;

    bit = base->findRec( "parameters" );
    if ( bit == base->end( "parameters" ) ) {
      throw "No parameters given.";
    }
    paramBlock = *bit;
    omega_ = paramBlock.getParam<T>( "omega" );
    cSmagorinsky_ = paramBlock.getParam<T>( "cSmagorinsky" );
    maxSteps_ = paramBlock.getParam<int>( "maxSteps" );
    std::cout << "Read parameter specification:" << std::endl;
    std::cout << "omega                : " << omega_ << std::endl;
    std::cout << "Smagorinsky constant : " << cSmagorinsky_ << std::endl;
    std::cout << "Number of steps      : " << maxSteps_ << std::endl;

    bit = base->findRec( "vtk" );
    if ( bit != base->end( "vtk" ) ) {
      paramBlock = *bit;
      vtkStep_ = paramBlock.getParam<int>( "vtkStep" );
      vtkFileName_ = paramBlock.getParam<std::string>( "vtkFileName" );
      std::cout << "VTK output specification:" << std::endl;
      std::cout << "VTK step modulo           : " << vtkStep_ << std::endl;
      std::cout << "VTK output file base name : " << vtkFileName_ << std::endl;
    } else {
      std::cout << "No vtk block given in configuration file, no output will be created." << std::endl;
      vtkStep_ = 0;
    }

    std::cout << "Set up the lattices..." << std::endl;

    // Set up the lattices accordingly
    grid0_ = new Grid<T,Dim>( sizeX, sizeY, sizeZ );
    grid1_ = new Grid<T,Dim>( sizeX, sizeY, sizeZ );
    u_.init( sizeX, sizeY, sizeZ, 0. );
    rho_.init( sizeX, sizeY, sizeZ, 1. );
    flag_.init( sizeX, sizeY, sizeZ, UNDEFINED );

    std::cout << "Set up the boundaries..." << std::endl;

    // Set up the boundaries
    bit = base->findRec( "boundaries" );
    if ( bit == base->end( "boundaries" ) ) {
      throw "No boundary description given.";
    }
    paramBlock = *bit;

    std::cout << "South..." << std::endl;

    // South
    bit = paramBlock.findRec( "south" );
    if ( bit != paramBlock.end( "south" ) ) {
      ConfBlock* bd = (*bit).descend();
      for ( ConfBlockIterator it = bd->begin(); it != bd->end(); ++it ) {
        ConfBlock b = *it;
        int xStart = b.getParam<int>( "xStart" );
        int xEnd   = b.getParam<int>( "xEnd" );
        int zStart = b.getParam<int>( "zStart" );
        int zEnd   = b.getParam<int>( "zEnd" );
        assert( xStart > 0 && xEnd < sizeX - 1 && zStart > 1 && zEnd < sizeZ - 2 );
        if ( b.getName() == "noslip" ) {
          for ( int z = zStart; z <= zEnd; ++z ) {
            for ( int x = xStart; x <= xEnd; ++x ) {
              assert( flag_( x, 1, z ) == UNDEFINED );
              flag_( x, 1, z ) = NOSLIP;
              boundaryCells_.push_back( Vec3<int>( x, 1, z ) );
            }
          }
        } else if ( b.getName() == "velocity" ) {
          T u_x = b.getParam<T>( "u_x" );
          T u_y = b.getParam<T>( "u_y" );
          T u_z = b.getParam<T>( "u_z" );
          for ( int z = zStart; z <= zEnd; ++z ) {
            for ( int x = xStart; x <= xEnd; ++x ) {
              assert( flag_( x, 1, z ) == UNDEFINED );
              flag_( x, 1, z ) = VELOCITY;
              velocityCells_.push_back( Vec3<int>( x, 1, z ) );
              velocities_.push_back( Vec3<T>( u_x, u_y, u_z ) );
            }
          }
        } else {
          std::cerr << "Unknown boundary type " << b.getName() << " (south)" << std::endl;
          continue;
        }

      }
    }
    // Fill unspecified cells with no slip boundary
    for ( int z = 2; z < sizeZ - 2; ++z ) {
      for ( int x = 1; x <= sizeX - 2; ++x ) {
        if ( flag_( x, 1, z ) == UNDEFINED ) {
          flag_( x, 1, z ) = NOSLIP;
          boundaryCells_.push_back( Vec3<int>( x, 1, z ) );
        }
      }
    }

    std::cout << "North..." << std::endl;

    // North
    bit = paramBlock.findRec( "north" );
    if ( bit != paramBlock.end( "north" ) ) {
      ConfBlock* bd = (*bit).descend();
        for ( ConfBlockIterator it = bd->begin(); it != bd->end(); ++it ) {
        ConfBlock b = *it;
        int xStart = b.getParam<int>( "xStart" );
        int xEnd   = b.getParam<int>( "xEnd" );
        int zStart = b.getParam<int>( "zStart" );
        int zEnd   = b.getParam<int>( "zEnd" );
        assert( xStart > 0 && xEnd < sizeX - 1 && zStart > 1 && zEnd < sizeZ - 2 );
        if ( b.getName() == "noslip" ) {
          for ( int z = zStart; z <= zEnd; ++z ) {
            for ( int x = xStart; x <= xEnd; ++x ) {
              assert( flag_( x, sizeY - 2, z ) == UNDEFINED );
              flag_( x, sizeY - 2, z ) = NOSLIP;
              boundaryCells_.push_back( Vec3<int>( x, sizeY - 2, z ) );
            }
          }
        } else if ( b.getName() == "velocity" ) {
          T u_x = b.getParam<T>( "u_x" );
          T u_y = b.getParam<T>( "u_y" );
          T u_z = b.getParam<T>( "u_z" );
          for ( int z = zStart; z <= zEnd; ++z ) {
            for ( int x = xStart; x <= xEnd; ++x ) {
              assert( flag_( x, sizeY - 2, z ) == UNDEFINED );
              flag_( x, sizeY - 2, z ) = VELOCITY;
              velocityCells_.push_back( Vec3<int>( x, sizeY - 2, z ) );
              velocities_.push_back( Vec3<T>( u_x, u_y, u_z ) );
            }
          }
        } else {
          std::cerr << "Unknown boundary type " << b.getName() << " (north)" << std::endl;
          continue;
        }

      }
    }
    // Fill unspecified cells with no slip boundary
    for ( int z = 2; z < sizeZ - 2; ++z ) {
      for ( int x = 1; x <= sizeX - 2; ++x ) {
        if ( flag_( x, sizeY - 2, z ) == UNDEFINED ) {
          flag_( x, sizeY - 2, z ) = NOSLIP;
          boundaryCells_.push_back( Vec3<int>( x, sizeY - 2, z ) );
        }
      }
    }

    std::cout << "West..." << std::endl;

    // West
    bit = paramBlock.findRec( "west" );
    if ( bit != paramBlock.end( "west" ) ) {
      ConfBlock* bd = (*bit).descend();
        for ( ConfBlockIterator it = bd->begin(); it != bd->end(); ++it ) {
        ConfBlock b = *it;
        int yStart = b.getParam<int>( "yStart" );
        int yEnd   = b.getParam<int>( "yEnd" );
        int zStart = b.getParam<int>( "zStart" );
        int zEnd   = b.getParam<int>( "zEnd" );
        assert( yStart > 1 && yEnd < sizeY - 2 && zStart > 1 && zEnd < sizeZ - 2 );
        if ( b.getName() == "noslip" ) {
          for ( int z = zStart; z <= zEnd; ++z ) {
            for ( int y = yStart; y <= yEnd; ++y ) {
              assert( flag_( 1, y, z ) == UNDEFINED );
              flag_( 1, y, z ) = NOSLIP;
              boundaryCells_.push_back( Vec3<int>( 1, y, z ) );
            }
          }
        } else if ( b.getName() == "velocity" ) {
          T u_x = b.getParam<T>( "u_x" );
          T u_y = b.getParam<T>( "u_y" );
          T u_z = b.getParam<T>( "u_z" );
          for ( int z = zStart; z <= zEnd; ++z ) {
            for ( int y = yStart; y <= yEnd; ++y ) {
              assert( flag_( 1, y, z ) == UNDEFINED );
              flag_( 1, y, z ) = VELOCITY;
              velocityCells_.push_back( Vec3<int>( 1, y, z ) );
              velocities_.push_back( Vec3<T>( u_x, u_y, u_z ) );
            }
          }
        } else {
          std::cerr << "Unknown boundary type " << b.getName() << " (west)" << std::endl;
          continue;
        }

      }
    }
    // Fill unspecified cells with no slip boundary
    for ( int z = 2; z < sizeZ - 2; ++z ) {
      for ( int y = 2; y < sizeY - 2; ++y ) {
        if ( flag_( 1, y, z ) == UNDEFINED ) {
          flag_( 1, y, z ) = NOSLIP;
          boundaryCells_.push_back( Vec3<int>( 1, y, z ) );
        }
      }
    }

    std::cout << "East..." << std::endl;

    // East
    bit = paramBlock.findRec( "east" );
    if ( bit != paramBlock.end( "east" ) ) {
      ConfBlock* bd = (*bit).descend();
        for ( ConfBlockIterator it = bd->begin(); it != bd->end(); ++it ) {
        ConfBlock b = *it;
        int yStart = b.getParam<int>( "yStart" );
        int yEnd   = b.getParam<int>( "yEnd" );
        int zStart = b.getParam<int>( "zStart" );
        int zEnd   = b.getParam<int>( "zEnd" );
        assert( yStart > 1 && yEnd < sizeY - 2 && zStart > 1 && zEnd < sizeZ - 2 );
        if ( b.getName() == "noslip" ) {
          for ( int z = zStart; z <= zEnd; ++z ) {
            for ( int y = yStart; y <= yEnd; ++y ) {
              assert( flag_( sizeX - 2, y, z ) == UNDEFINED );
              flag_( sizeX - 2, y, z ) = NOSLIP;
              boundaryCells_.push_back( Vec3<int>( sizeX - 2, y, z ) );
            }
          }
        } else if ( b.getName() == "velocity" ) {
          T u_x = b.getParam<T>( "u_x" );
          T u_y = b.getParam<T>( "u_y" );
          T u_z = b.getParam<T>( "u_z" );
          for ( int z = zStart; z <= zEnd; ++z ) {
            for ( int y = yStart; y <= yEnd; ++y ) {
              assert( flag_( sizeX - 2, y, z ) == UNDEFINED );
              flag_( sizeX - 2, y, z ) = VELOCITY;
              velocityCells_.push_back( Vec3<int>( sizeX - 2, y, z ) );
              velocities_.push_back( Vec3<T>( u_x, u_y, u_z ) );
            }
          }
        } else {
          std::cerr << "Unknown boundary type " << b.getName() << " (east)" << std::endl;
          continue;
        }

      }
    }
    // Fill unspecified cells with no slip boundary
    for ( int z = 2; z < sizeZ - 2; ++z ) {
      for ( int y = 2; y < sizeY - 2; ++y ) {
        if ( flag_( sizeX - 2, y, z ) == UNDEFINED ) {
          flag_( sizeX - 2, y, z ) = NOSLIP;
          boundaryCells_.push_back( Vec3<int>( sizeX - 2, y, z ) );
        }
      }
    }

    std::cout << "Bottom..." << std::endl;

    // Bottom
    bit = paramBlock.findRec( "bottom" );
    if ( bit != paramBlock.end( "bottom" ) ) {
      ConfBlock* bd = (*bit).descend();
        for ( ConfBlockIterator it = bd->begin(); it != bd->end(); ++it ) {
        ConfBlock b = *it;
        int xStart = b.getParam<int>( "xStart" );
        int xEnd   = b.getParam<int>( "xEnd" );
        int yStart = b.getParam<int>( "yStart" );
        int yEnd   = b.getParam<int>( "yEnd" );
        assert( xStart > 0 && xEnd < sizeX - 1 && yStart > 0 && yEnd < sizeY - 1 );
        if ( b.getName() == "noslip" ) {
          for ( int y = yStart; y <= yEnd; ++y ) {
            for ( int x = xStart; x <= xEnd; ++x ) {
              assert( flag_( x, y, 1 ) == UNDEFINED );
              flag_( x, y, 1 ) = NOSLIP;
              boundaryCells_.push_back( Vec3<int>( x, y, 1 ) );
            }
          }
        } else if ( b.getName() == "velocity" ) {
          T u_x = b.getParam<T>( "u_x" );
          T u_y = b.getParam<T>( "u_y" );
          T u_z = b.getParam<T>( "u_z" );
          for ( int y = yStart; y <= yEnd; ++y ) {
            for ( int x = xStart; x <= xEnd; ++x ) {
              assert( flag_( x, y, 1 ) == UNDEFINED );
              flag_( x, y, 1 ) = VELOCITY;
              velocityCells_.push_back( Vec3<int>( x, y, 1 ) );
              velocities_.push_back( Vec3<T>( u_x, u_y, u_z ) );
            }
          }
        } else {
          std::cerr << "Unknown boundary type " << b.getName() << " (bottom)" << std::endl;
          continue;
        }

      }
    }
    // Fill unspecified cells with no slip boundary
    for ( int y = 1; y <= sizeY - 2; ++y ) {
      for ( int x = 1; x <= sizeX - 2; ++x ) {
        if ( flag_( x, y, 1 ) == UNDEFINED ) {
          flag_( x, y, 1 ) = NOSLIP;
          boundaryCells_.push_back( Vec3<int>( x, y, 1 ) );
        }
      }
    }

    std::cout << "Top..." << std::endl;

    // Top
    bit = paramBlock.findRec( "top" );
    if ( bit != paramBlock.end( "top" ) ) {
      ConfBlock* bd = (*bit).descend();
        for ( ConfBlockIterator it = bd->begin(); it != bd->end(); ++it ) {
        ConfBlock b = *it;
        int xStart = b.getParam<int>( "xStart" );
        int xEnd   = b.getParam<int>( "xEnd" );
        int yStart = b.getParam<int>( "yStart" );
        int yEnd   = b.getParam<int>( "yEnd" );
        assert( xStart > 0 && xEnd < sizeX - 1 && yStart > 0 && yEnd < sizeY - 1 );
        if ( b.getName() == "noslip" ) {
          for ( int y = yStart; y <= yEnd; ++y ) {
            for ( int x = xStart; x <= xEnd; ++x ) {
              assert( flag_( x, y, sizeZ - 2 ) == UNDEFINED );
              flag_( x, y, sizeZ - 2 ) = NOSLIP;
              boundaryCells_.push_back( Vec3<int>( x, y, sizeZ - 2 ) );
            }
          }
        } else if ( b.getName() == "velocity" ) {
          T u_x = b.getParam<T>( "u_x" );
          T u_y = b.getParam<T>( "u_y" );
          T u_z = b.getParam<T>( "u_z" );
          for ( int y = yStart; y <= yEnd; ++y ) {
            for ( int x = xStart; x <= xEnd; ++x ) {
              assert( flag_( x, y, sizeZ - 2 ) == UNDEFINED );
              flag_( x, y, sizeZ - 2 ) = VELOCITY;
              velocityCells_.push_back( Vec3<int>( x, y, sizeZ - 2 ) );
              velocities_.push_back( Vec3<T>( u_x, u_y, u_z ) );
            }
          }
        } else {
          std::cerr << "Unknown boundary type " << b.getName() << " (top)" << std::endl;
          continue;
        }

      }
    }
    // Fill unspecified cells with no slip boundary
    for ( int y = 1; y <= sizeY - 2; ++y ) {
      for ( int x = 1; x <= sizeX - 2; ++x ) {
        if ( flag_( x, y, sizeZ - 2 ) == UNDEFINED ) {
          flag_( x, y, sizeZ - 2 ) = NOSLIP;
          boundaryCells_.push_back( Vec3<int>( x, y, sizeZ - 2 ) );
        }
      }
    }

    std::cout << "Flag fluid cells..." << std::endl;

    // Non-boundary cells are fluid cells
    for ( int z = 2; z < sizeZ - 2; ++z )
      for ( int y = 2; y < sizeY - 2; ++y )
        for ( int x = 2; x < sizeX - 2; ++x )
          flag_( x, y, z ) = FLUID;

    std::cout << "Initialize distribution functions with equilibrium..." << std::endl;

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

//    delete base;

  } catch ( BadSyntax e ) {
    std::cerr << e.what() << std::endl;
    exit( -1 );
  } catch ( ParameterNotFound e ) {
    std::cerr << e.what() << std::endl;
    exit( -1 );
  } catch ( const char* e ) {
    std::cerr << e << std::endl;
    exit( -1 );
  }

  std::cout << "Setup finished!" << std::endl;
}

template<typename T>
inline T LBM<T>::getTime( timeval &start, timeval &end ) {
  return (T) ( end.tv_sec - start.tv_sec )
          + (T) ( end.tv_usec - start.tv_usec ) / 1000000.;
}

template<typename T>
inline void LBM<T>::collideStream( int x, int y, int z, T omega ) {

  // calculate rho and u
  T rho = (*grid0_)( x, y, z, 0 ); // df in center
  T ux = 0.;
  T uy = 0.;
  T uz = 0.;
  // loop over all velocity directions but center
  for ( int f = 1; f < Dim; ++f ) {
    T fi = (*grid0_)( x, y, z, f );
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
  T fc = rho - 1.5 * ( ux * ux + uy * uy + uz * uz );
  T omegai = 1 - omega;
  // treat center value specially
  (*grid1_)( x, y, z, 0 ) = omegai * (*grid0_)( x, y, z, 0 )
                        + omega *  w[0] * fc;
  // loop over all velocity directions but center
  for ( int f = 1; f < Dim; ++f ) {
    T eiu = ex[f] * ux + ey[f] * uy + ez[f] * uz;
    (*grid1_)( x + ex[f], y + ey[f], z + ez[f], f )
      =   omegai * (*grid0_)( x, y, z, f )
        + omega  * w[f] * ( fc +  3 * eiu + 4.5 * eiu * eiu);
  }
}

template<typename T>
inline void LBM<T>::collideStreamSmagorinsky( int x, int y, int z, T nu, T cSqr ) {

  // Calculate rho and u
  T rho = (*grid0_)( x, y, z, 0 ); // df in center
  T ux = 0.;
  T uy = 0.;
  T uz = 0.;
  // Loop over all velocity directions but center
  for ( int f = 1; f < Dim; ++f ) {
    T fi = (*grid0_)( x, y, z, f );
    rho += fi;
    ux += ex[f] * fi;
    uy += ey[f] * fi;
    uz += ez[f] * fi;
  }
  // DEBUG assertions
  assert ( rho > 0.5 && rho < 1.5 );
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
  T fc = rho - 1.5 * ( ux * ux + uy * uy + uz * uz );
  T feq[19];

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
  T qo = 0.;
  for ( int i = 0; i < 3; ++i ) {
    T qadd = 0.;
    for ( int f = 1; f < 19; ++f ) {
      qadd += ep[i][f] * ( (*grid0_)( x, y, z, f ) - feq[f] );
    }
    qo += qadd * qadd;
  }
  qo *= 2.;
  for ( int i = 4; i < 6; ++i ) {
    T qadd = 0.;
    for ( int f = 7; f < 19; ++f ) {
      qadd += ep[i][f] * ( (*grid0_)( x, y, z, f ) - feq[f] );
    }
    qo += qadd * qadd;
  }
  qo = sqrt( qo );

  // Calculate local stress tensor
  T s = ( sqrt( nu * nu + 18. * cSqr * qo ) - nu ) / ( 6. * cSqr);
  // Calculate turbulence modified inverse lattice viscosity
  T omega = 1. / ( 3. * ( nu + cSqr * s ) + .5 );
  T omegai = 1. - omega;

  // Loop over all velocity directions and stream collided distribution value
  // to neighboring cells
  for ( int f = 0; f < Dim; ++f ) {
    (*grid1_)( x + ex[f], y + ey[f], z + ez[f], f )
      =   omegai * (*grid0_)( x, y, z, f ) + omega * feq[f];
  }
}

template<typename T>
inline void LBM<T>::treatBoundary() {

  // Iterate over all no-slip boundary cells
  std::vector< Vec3<int> >::iterator iter;
  for( iter = boundaryCells_.begin(); iter != boundaryCells_.end(); iter++ ) {

    // Fetch coordinates of current boundary cell
    int x = (*iter)[0];
    int y = (*iter)[1];
    int z = (*iter)[2];

    // Go over all distribution values and stream to inverse distribution
    // value of adjacent cell in inverse direction (bounce back)
    for ( int f = 1; f < Dim; ++f ) {
      (*grid1_)( x - ex[f], y - ey[f], z - ez[f], finv[f] ) = (*grid1_)( x, y, z, f );
    }
  }
}

template<typename T>
inline void LBM<T>::treatVelocities() {

  // Iterate over all velocity boundary cells
  for( int i = 0; i < (int) velocityCells_.size(); ++i ) {

    // Fetch coordinates of current boundary cell
    int x = velocityCells_[i][0];
    int y = velocityCells_[i][1];
    int z = velocityCells_[i][2];
    // Fetch velocity of moving wall
    T ux = velocities_[i][0];
    T uy = velocities_[i][1];
    T uz = velocities_[i][2];
    // Fetch density of current cell
    T rho = 6 * rho_( x, y, z );
    // Set velocity of this cell
    u_( x, y, z, 0 ) = ux;
    u_( x, y, z, 1 ) = uy;
    u_( x, y, z, 2 ) = uz;

    // Go over all distribution values, stream to inverse distribution value of
    // adjacent cell in inverse direction (bounce back) and modify by velocity
    // of moving wall
    for ( int f = 1; f < Dim; ++f ) {
      int op = finv[f];
      (*grid1_)( x - ex[f], y - ey[f], z - ez[f], op )
        = (*grid1_)( x, y, z, f )
          + rho * w[f] * ( ex[op] * ux + ey[op] * uy + ez[op] * uz );
    }
  }
}

template<>
void LBM<double>::writeVtkFile( int timestep ) {

  // Open file for writing
  std::ostringstream oss;
  oss << vtkFileName_ << "." << timestep << ".vtk";
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

  // Write flag field
  vtkFile << "SCALARS flags int\n";
  vtkFile << "LOOKUP_TABLE default\n";
  for ( int z = 1; z <= sizeZ; ++z ) {
    for ( int y = 1; y <= sizeY; ++y ) {
      for ( int x = 1; x <= sizeX; ++x ) {
        // evil hack because vtk requires binary data to be in big endian
        uint32_t dump = htobe32( *reinterpret_cast<uint32_t *>( &flag_( x, y, z ) ) );
        vtkFile.write( reinterpret_cast<char *>( &dump ), sizeof(int) );
      }
    }
  }

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

template<>
void LBM<float>::writeVtkFile( int timestep ) {

  // Open file for writing
  std::ostringstream oss;
  oss << vtkFileName_ << "." << timestep << ".vtk";
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

  // Write flag field
  vtkFile << "SCALARS flags int\n";
  vtkFile << "LOOKUP_TABLE default\n";
  for ( int z = 1; z <= sizeZ; ++z ) {
    for ( int y = 1; y <= sizeY; ++y ) {
      for ( int x = 1; x <= sizeX; ++x ) {
        // evil hack because vtk requires binary data to be in big endian
        uint32_t dump = htobe32( *reinterpret_cast<uint32_t *>( &flag_( x, y, z ) ) );
        vtkFile.write( reinterpret_cast<char *>( &dump ), sizeof(int) );
      }
    }
  }

  // Write density field
  vtkFile << "SCALARS density float\n";
  vtkFile << "LOOKUP_TABLE default\n";
  for ( int z = 1; z <= sizeZ; ++z ) {
    for ( int y = 1; y <= sizeY; ++y ) {
      for ( int x = 1; x <= sizeX; ++x ) {
        // evil hack because vtk requires binary data to be in big endian
        uint32_t dump = htobe32( *reinterpret_cast<uint32_t *>( &rho_( x, y, z ) ) );
        vtkFile.write( reinterpret_cast<char *>( &dump ), sizeof(float) );
      }
    }
  }

  // Write velocity vector field
  vtkFile << "VECTORS velocity float\n";
  for ( int z = 1; z <= sizeZ; ++z ) {
    for ( int y = 1; y <= sizeY; ++y ) {
      for ( int x = 1; x <= sizeX; ++x ) {
        // evil hack because vtk requires binary data to be in big endian
        uint32_t dump0 = htobe32( *reinterpret_cast<uint32_t *>( &u_( x, y, z, 0 ) ) );
        uint32_t dump1 = htobe32( *reinterpret_cast<uint32_t *>( &u_( x, y, z, 1 ) ) );
        uint32_t dump2 = htobe32( *reinterpret_cast<uint32_t *>( &u_( x, y, z, 2 ) ) );
        vtkFile.write( reinterpret_cast<char *>( &dump0 ), sizeof(float) );
        vtkFile.write( reinterpret_cast<char *>( &dump1 ), sizeof(float) );
        vtkFile.write( reinterpret_cast<char *>( &dump2 ), sizeof(float) );
      }
    }
  }

}

} // namespace lbm

#endif /* LBM_DEF_H_ */
