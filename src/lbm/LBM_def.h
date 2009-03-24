//! \file LBM_def.h
//! Implementation of the LBM class

//! \date   Jan 16, 2009
//! \author Florian Rathgeber

#ifndef LBM_DEF_H_
#define LBM_DEF_H_

#ifndef htobe32

  #if __BYTE_ORDER == __LITTLE_ENDIAN
    #include <byteswap.h>
    #define htobe32(x) bswap_32(x)
    #define htobe64(x) bswap_64(x)
  #else
    #define htobe32(x) (x)
    #define htobe64(x) (x)
  #endif

#endif

namespace lbm {

// ============================ //
// Constructors and destructors //
// ============================ //

template<typename T>
LBM<T>::LBM( const std::string configFileName ) : curStep_( 0 ) {

  try {
    ConfParser p;
    ConfBlock base = p.parse( configFileName );
    std::cout << "Parsed configuration file " << configFileName << std::endl;
    setup( base );
  } catch ( std::exception e ) {
    std::cerr << e.what() << std::endl;
    exit( -1 );
  }
}

template<typename T>
LBM<T>::LBM( ConfBlock& base ) : curStep_( 0 ) {
  setup( base );
}

template<typename T>
LBM<T>::~LBM() {
  assert ( grid0_ != 0 && grid1_ != 0 );
  delete grid0_;
  delete grid1_;
}

template<typename T>
void LBM<T>::run() {

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
    totTime += runStep();
  }

  std::cout << "LBM finished! Processed " << maxSteps_ << " timesteps in ";
  std::cout << totTime << " secs!" << std::endl;
  std::cout << "Average speed of " << (maxSteps_ * numCells) / (totTime
      * 1000000);
  std::cout << " MLUP/s" << std::endl;
}

template<typename T>
double LBM<T>::runStep() {

  // Get effective domain size (without ghost layer)
  int sizeX = grid0_->getSizeX() - 2;
  int sizeY = grid0_->getSizeY() - 2;
  int sizeZ = grid0_->getSizeZ() - 2;

  int numCells = sizeX * sizeY * sizeZ;

  struct timeval start, end;
  T stepTime = 0.;

  gettimeofday(&start, NULL);

  // loop over all cells but the ghost layer
  for (int z = 1; z <= sizeZ - 1; z++) {
    for (int y = 1; y <= sizeY - 1; y++) {
      for (int x = 1; x <= sizeX - 1; x++) {

#ifdef NSMAGO
        // Perform actual collision and streaming step
        collideStream( x, y, z );
#else
        // Perform collision and streaming step with Smagorinsky turbulence
        // correction
        collideStreamSmagorinsky( x, y, z );
#endif

      } // x
    } // y
  } // z

  gettimeofday(&end, NULL);
  T scTime = getTime(start, end);

  // Treat no-slip boundary conditions (walls)
  treatNoslip();

  gettimeofday(&start, NULL);
  T nTime = getTime(end, start);

  // Treat velocity cells
  treatVelocity();

  gettimeofday(&end, NULL);
  T vTime = getTime(start, end);

  // Treat inflow boundary conditions
  treatInflow();

  gettimeofday(&start, NULL);
  T iTime = getTime(end, start);

  // Treat outflow boundary conditions
  treatOutflow();

  gettimeofday(&end, NULL);
  T oTime = getTime(start, end);

  // Treat pressure cells
  treatPressure();

  gettimeofday(&start, NULL);
  T pTime = getTime(end, start);

  // exchange grids for current and previous time step
  Grid<T, Dim> *gridTmp = grid0_;
  grid0_ = grid1_;
  grid1_ = gridTmp;

  stepTime = scTime + nTime + vTime + iTime + oTime + pTime;

#ifdef DEBUG
  std::cout << "Time step " << curStep_ << " of " << maxSteps_ << " took ";
  std::cout << scTime << " + " << nTime << " + " << vTime << " + " << iTime;
  std::cout << " + " << oTime << " + " << pTime << " = " << stepTime;
  std::cout << "secs -> " << numCells / (stepTime * 1000000) << " MLUP/s";
  std::cout << std::endl;
#endif

  if ( vtkStep_ != 0 && curStep_ % vtkStep_ == 0 )
    writeVtkFile();

  ++curStep_;
  return stepTime;
}

template<typename T>
inline Vec3<T> LBM<T>::getVelocity( T x, T y, T z ) {
  int xf = (int) x;
  int yf = (int) y;
  int zf = (int) z;
  int xc = xf + 1;
  int yc = yf + 1;
  int zc = zf + 1;
  T xd = x - xf;
  T yd = y - yf;
  T zd = z - zf;
  T u_x = (1. - xd) * (
        (1. - yd) * ( u_(xf, yf, zf, 0) * (1. - zd) + u_(xf, yf, zc, 0) * zd )
      + yd        * ( u_(xf, yc, zf, 0) * (1. - zd) + u_(xf, yc, zc, 0) * zd )
                      ) + xd * (
        (1. - yd) * ( u_(xc, yf, zf, 0) * (1. - zd) + u_(xc, yf, zc, 0) * zd )
      + yd        * ( u_(xc, yc, zf, 0) * (1. - zd) + u_(xc, yc, zc, 0) * zd )
                               );
  T u_y = (1. - xd) * (
        (1. - yd) * ( u_(xf, yf, zf, 1) * (1. - zd) + u_(xf, yf, zc, 1) * zd )
      + yd        * ( u_(xf, yc, zf, 1) * (1. - zd) + u_(xf, yc, zc, 1) * zd )
                      ) + xd * (
        (1. - yd) * ( u_(xc, yf, zf, 1) * (1. - zd) + u_(xc, yf, zc, 1) * zd )
      + yd        * ( u_(xc, yc, zf, 1) * (1. - zd) + u_(xc, yc, zc, 1) * zd )
                               );
  T u_z = (1. - xd) * (
        (1. - yd) * ( u_(xf, yf, zf, 2) * (1. - zd) + u_(xf, yf, zc, 2) * zd )
      + yd        * ( u_(xf, yc, zf, 2) * (1. - zd) + u_(xf, yc, zc, 2) * zd )
                      ) + xd * (
        (1. - yd) * ( u_(xc, yf, zf, 2) * (1. - zd) + u_(xc, yf, zc, 2) * zd )
      + yd        * ( u_(xc, yc, zf, 2) * (1. - zd) + u_(xc, yc, zc, 2) * zd )
                               );
//   std::cout << "Velocity at <" << x << "," << y << "," << z << "> in timestep " << curStep_ << ": <";
//   std::cout << u_x << "," << u_y << "," << u_z << ">" << std::endl;
  return Vec3<T>( u_x, u_y, u_z );
}

// ========================= //
// Internal helper functions //
// ========================= //

template<typename T>
void LBM<T>::setup( ConfBlock& base ) {
  try {

    // Read the parameters from the config file

    std::cout << "Setting up LBM..." << std::endl;

    ConfBlock* paramBlock = base.find( "domain" );
    if ( paramBlock == NULL ) {
      throw "No domain size given.";
    }
    int sizeX = paramBlock->getParam<int>( "sizeX" );
    int sizeY = paramBlock->getParam<int>( "sizeY" );
    int sizeZ = paramBlock->getParam<int>( "sizeZ" );
    std::cout << "Read domain specification:" << std::endl;
    std::cout << "sizeX : " << sizeX << std::endl;
    std::cout << "sizeY : " << sizeY << std::endl;
    std::cout << "sizeZ : " << sizeZ << std::endl;

    paramBlock = base.find( "parameters" );
    if ( paramBlock == NULL ) {
      throw "No parameters given.";
    }
    omega_ = paramBlock->getParam<T>( "omega" );
#ifndef NSMAGO
    T cSmagorinsky = paramBlock->getParam<T>( "cSmagorinsky" );
#endif
    maxSteps_ = paramBlock->getParam<int>( "maxSteps" );
    std::cout << "Read parameter specification:" << std::endl;
    std::cout << "omega                : " << omega_ << std::endl;
#ifndef NSMAGO
    std::cout << "Smagorinsky constant : " << cSmagorinsky << std::endl;
#endif
    std::cout << "Number of steps      : " << maxSteps_ << std::endl;

#ifndef NSMAGO
    // lattice viscosity
    nu_ = (2. / omega_ - 1.) * (1. / 6.);
    // squared Smagorinsky constant
    cSmagoSqr_ = cSmagorinsky * cSmagorinsky;
#endif

    paramBlock = base.find( "vtk" );
    if ( paramBlock != NULL ) {
      vtkStep_ = paramBlock->getParam<int>( "vtkStep" );
      vtkFileName_ = paramBlock->getParam<std::string>( "vtkFileName" );
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
    paramBlock = base.find( "boundaries" );
    if ( paramBlock == NULL ) {
      throw "No boundary description given.";
    }

    std::cout << "Bottom..." << std::endl;

    // Bottom
    ConfBlock* boundBlock = paramBlock->find( "bottom" );
    if ( boundBlock != NULL ) setupBoundary( *boundBlock, -1, 1, -1 );
    // Fill unspecified cells with no slip boundary
    for ( int z = 2; z < sizeZ - 2; ++z ) {
      for ( int x = 1; x <= sizeX - 2; ++x ) {
        if ( flag_( x, 1, z ) == UNDEFINED ) {
          flag_( x, 1, z ) = NOSLIP;
          noslipCells_.push_back( Vec3<int>( x, 1, z ) );
        }
      }
    }

    std::cout << "Top..." << std::endl;

    // Top
    boundBlock = paramBlock->find( "top" );
    if ( boundBlock != NULL ) setupBoundary( *boundBlock, -1, sizeY - 2, -1 );
    // Fill unspecified cells with no slip boundary
    for ( int z = 2; z < sizeZ - 2; ++z ) {
      for ( int x = 1; x <= sizeX - 2; ++x ) {
        if ( flag_( x, sizeY - 2, z ) == UNDEFINED ) {
          flag_( x, sizeY - 2, z ) = NOSLIP;
          noslipCells_.push_back( Vec3<int>( x, sizeY - 2, z ) );
        }
      }
    }

    std::cout << "Left..." << std::endl;

    // Left
    boundBlock = paramBlock->find( "left" );
    if ( boundBlock != NULL ) setupBoundary( *boundBlock, 1, -1, -1 );
    // Fill unspecified cells with no slip boundary
    for ( int z = 2; z < sizeZ - 2; ++z ) {
      for ( int y = 2; y < sizeY - 2; ++y ) {
        if ( flag_( 1, y, z ) == UNDEFINED ) {
          flag_( 1, y, z ) = NOSLIP;
          noslipCells_.push_back( Vec3<int>( 1, y, z ) );
        }
      }
    }

    std::cout << "Right..." << std::endl;

    // Right
    boundBlock = paramBlock->find( "right" );
    if ( boundBlock != NULL ) setupBoundary( *boundBlock, sizeX - 2, -1, -1 );
    // Fill unspecified cells with no slip boundary
    for ( int z = 2; z < sizeZ - 2; ++z ) {
      for ( int y = 2; y < sizeY - 2; ++y ) {
        if ( flag_( sizeX - 2, y, z ) == UNDEFINED ) {
          flag_( sizeX - 2, y, z ) = NOSLIP;
          noslipCells_.push_back( Vec3<int>( sizeX - 2, y, z ) );
        }
      }
    }

    std::cout << "Front..." << std::endl;

    // Front
    boundBlock = paramBlock->find( "front" );
    if ( boundBlock != NULL ) setupBoundary( *boundBlock, -1, -1, 1 );
    // Fill unspecified cells with no slip boundary
    for ( int y = 1; y <= sizeY - 2; ++y ) {
      for ( int x = 1; x <= sizeX - 2; ++x ) {
        if ( flag_( x, y, 1 ) == UNDEFINED ) {
          flag_( x, y, 1 ) = NOSLIP;
          noslipCells_.push_back( Vec3<int>( x, y, 1 ) );
        }
      }
    }

    std::cout << "Back..." << std::endl;

    // Back
    boundBlock = paramBlock->find( "back" );
    if ( boundBlock != NULL ) setupBoundary( *boundBlock, -1, -1, sizeZ - 2 );
    // Fill unspecified cells with no slip boundary
    for ( int y = 1; y <= sizeY - 2; ++y ) {
      for ( int x = 1; x <= sizeX - 2; ++x ) {
        if ( flag_( x, y, sizeZ - 2 ) == UNDEFINED ) {
          flag_( x, y, sizeZ - 2 ) = NOSLIP;
          noslipCells_.push_back( Vec3<int>( x, y, sizeZ - 2 ) );
        }
      }
    }

    // Set up the obstacles
    paramBlock = base.find( "obstacles" );
    if ( paramBlock == NULL ) {
      std::cout << "No obstacles defined" << std::endl;
    } else {

      std::cout << "Set up the obstacles..." << std::endl;

      ConfBlock::childIterPair bit = paramBlock->findAll( "cuboid_stationary" );
      for ( ConfBlock::childIter it = bit.first; it != bit.second; ++it ) {

        ConfBlock& b = it->second;

        int xStart = b.getParam<int>( "xStart" );
        int xEnd   = b.getParam<int>( "xEnd" );
        int yStart = b.getParam<int>( "yStart" );
        int yEnd   = b.getParam<int>( "yEnd" );
        int zStart = b.getParam<int>( "zStart" );
        int zEnd   = b.getParam<int>( "zEnd" );
        std::cout << "Stationary cuboid ranging from <" << xStart << ",";
        std::cout << yStart << "," << zStart << "> to <" << xEnd << "," << yEnd;
        std::cout << "," << zEnd << ">" << std::endl;
        assert( xStart > 1 && xEnd < flag_.getSizeX() - 2 &&
                yStart > 1 && yEnd < flag_.getSizeY() - 2 &&
                zStart > 1 && zEnd < flag_.getSizeZ() - 2 );
        for ( int z = zStart; z <= zEnd; ++z ) {
          for ( int y = yStart; y <= yEnd; ++y ) {
            for ( int x = xStart; x <= xEnd; ++x ) {
              assert( flag_( x, y, z ) == UNDEFINED );
              flag_( x, y, z ) = NOSLIP;
              noslipCells_.push_back( Vec3<int>( x, y, z ) );
            }
          }
        }
      }

      bit = paramBlock->findAll( "inflow" );
      for ( ConfBlock::childIter it = bit.first; it != bit.second; ++it ) {

        ConfBlock& b = it->second;

        int xStart = b.getParam<int>( "xStart" );
        int xEnd   = b.getParam<int>( "xEnd" );
        int yStart = b.getParam<int>( "yStart" );
        int yEnd   = b.getParam<int>( "yEnd" );
        int zStart = b.getParam<int>( "zStart" );
        int zEnd   = b.getParam<int>( "zEnd" );
        T u_x      = b.getParam<T>( "u_x" );
        T u_y      = b.getParam<T>( "u_y" );
        T u_z      = b.getParam<T>( "u_z" );
        assert( xStart > 1 && xEnd < flag_.getSizeX() - 2 &&
                yStart > 1 && yEnd < flag_.getSizeY() - 2 &&
                zStart > 1 && zEnd < flag_.getSizeZ() - 2 );
        for ( int z = zStart; z <= zEnd; ++z ) {
          for ( int y = yStart; y <= yEnd; ++y ) {
            for ( int x = xStart; x <= xEnd; ++x ) {
              assert( flag_( x, y, z ) == UNDEFINED );
              flag_( x, y, z ) = NOSLIP;
              inflowCells_.push_back( Vec3<int>( x, y, z ) );
              inflowVels_.push_back( Vec3<T>( u_x, u_y, u_z ) );
            }
          }
        }
      }
    }

    std::cout << "Flag fluid cells..." << std::endl;

    // Non-boundary cells and non-obstacle cells are fluid cells
    for ( int z = 2; z < sizeZ - 2; ++z )
      for ( int y = 2; y < sizeY - 2; ++y )
        for ( int x = 2; x < sizeX - 2; ++x )
          if ( flag_( x, y, z ) == UNDEFINED ) flag_( x, y, z ) = FLUID;

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

  } catch ( std::exception e ) {
    std::cerr << e.what() << std::endl;
    exit( -1 );
  } catch ( const char* e ) {
    std::cerr << e << std::endl;
    exit( -1 );
  }

  std::cout << "LBM setup finished!" << std::endl;
}

template<typename T>
inline void LBM<T>::setupBoundary( ConfBlock& block, int x, int y, int z ) {

  ConfBlock::childIterPair bit = block.findAll( "noslip" );

  for ( ConfBlock::childIter it = bit.first; it != bit.second; ++it ) {

    ConfBlock& b = it->second;

    int xStart = ( x == -1 ) ? b.getParam<int>( "xStart" ) : x;
    int xEnd   = ( x == -1 ) ? b.getParam<int>( "xEnd" )   : x;
    int yStart = ( y == -1 ) ? b.getParam<int>( "yStart" ) : y;
    int yEnd   = ( y == -1 ) ? b.getParam<int>( "yEnd" )   : y;
    int zStart = ( z == -1 ) ? b.getParam<int>( "zStart" ) : z;
    int zEnd   = ( z == -1 ) ? b.getParam<int>( "zEnd" )   : z;
    assert( xStart > 0 && xEnd < flag_.getSizeX() - 1 &&
            yStart > 0 && yEnd < flag_.getSizeY() - 1 &&
            zStart > 0 && zEnd < flag_.getSizeZ() - 1 );
    for ( z = zStart; z <= zEnd; ++z ) {
      for ( y = yStart; y <= yEnd; ++y ) {
        for ( x = xStart; x <= xEnd; ++x ) {
          assert( flag_( x, y, z ) == UNDEFINED );
          flag_( x, y, z ) = NOSLIP;
          noslipCells_.push_back( Vec3<int>( x, y, z ) );
        }
      }
    }

  }

  bit = block.findAll( "velocity" );

  for ( ConfBlock::childIter it = bit.first; it != bit.second; ++it ) {

    ConfBlock& b = it->second;

    int xStart = ( x == -1 ) ? b.getParam<int>( "xStart" ) : x;
    int xEnd   = ( x == -1 ) ? b.getParam<int>( "xEnd" )   : x;
    int yStart = ( y == -1 ) ? b.getParam<int>( "yStart" ) : y;
    int yEnd   = ( y == -1 ) ? b.getParam<int>( "yEnd" )   : y;
    int zStart = ( z == -1 ) ? b.getParam<int>( "zStart" ) : z;
    int zEnd   = ( z == -1 ) ? b.getParam<int>( "zEnd" )   : z;
    T u_x = b.getParam<T>( "u_x" );
    T u_y = b.getParam<T>( "u_y" );
    T u_z = b.getParam<T>( "u_z" );
    assert( xStart > 0 && xEnd < flag_.getSizeX() - 1 &&
            yStart > 0 && yEnd < flag_.getSizeY() - 1 &&
            zStart > 0 && zEnd < flag_.getSizeZ() - 1 );
    for ( z = zStart; z <= zEnd; ++z ) {
      for ( y = yStart; y <= yEnd; ++y ) {
        for ( x = xStart; x <= xEnd; ++x ) {
          assert( flag_( x, y, z ) == UNDEFINED );
          flag_( x, y, z ) = VELOCITY;
          velocityCells_.push_back( Vec3<int>( x, y, z ) );
          velocityVels_.push_back( Vec3<T>( u_x, u_y, u_z ) );
        }
      }
    }

  }

  bit = block.findAll( "pressure" );

  for ( ConfBlock::childIter it = bit.first; it != bit.second; ++it ) {

    ConfBlock& b = it->second;

    int xStart = ( x == -1 ) ? b.getParam<int>( "xStart" ) : x;
    int xEnd   = ( x == -1 ) ? b.getParam<int>( "xEnd" )   : x;
    int yStart = ( y == -1 ) ? b.getParam<int>( "yStart" ) : y;
    int yEnd   = ( y == -1 ) ? b.getParam<int>( "yEnd" )   : y;
    int zStart = ( z == -1 ) ? b.getParam<int>( "zStart" ) : z;
    int zEnd   = ( z == -1 ) ? b.getParam<int>( "zEnd" )   : z;
    assert( xStart > 0 && xEnd < flag_.getSizeX() - 1 &&
            yStart > 0 && yEnd < flag_.getSizeY() - 1 &&
            zStart > 0 && zEnd < flag_.getSizeZ() - 1 );
    for ( z = zStart; z <= zEnd; ++z ) {
      for ( y = yStart; y <= yEnd; ++y ) {
        for ( x = xStart; x <= xEnd; ++x ) {
          assert( flag_( x, y, z ) == UNDEFINED );
          flag_( x, y, z ) = PRESSURE;
          pressureCells_.push_back( Vec3<int>( x, y, z ) );
        }
      }
    }

  }

  bit = block.findAll( "inflow" );

  for ( ConfBlock::childIter it = bit.first; it != bit.second; ++it ) {

    ConfBlock& b = it->second;

    int xStart = ( x == -1 ) ? b.getParam<int>( "xStart" ) : x;
    int xEnd   = ( x == -1 ) ? b.getParam<int>( "xEnd" )   : x;
    int yStart = ( y == -1 ) ? b.getParam<int>( "yStart" ) : y;
    int yEnd   = ( y == -1 ) ? b.getParam<int>( "yEnd" )   : y;
    int zStart = ( z == -1 ) ? b.getParam<int>( "zStart" ) : z;
    int zEnd   = ( z == -1 ) ? b.getParam<int>( "zEnd" )   : z;
    T u_x = b.getParam<T>( "u_x" );
    T u_y = b.getParam<T>( "u_y" );
    T u_z = b.getParam<T>( "u_z" );
    assert( xStart > 0 && xEnd < flag_.getSizeX() - 1 &&
            yStart > 0 && yEnd < flag_.getSizeY() - 1 &&
            zStart > 0 && zEnd < flag_.getSizeZ() - 1 );
    for ( z = zStart; z <= zEnd; ++z ) {
      for ( y = yStart; y <= yEnd; ++y ) {
        for ( x = xStart; x <= xEnd; ++x ) {
          assert( flag_( x, y, z ) == UNDEFINED );
          flag_( x, y, z ) = INFLOW;
          inflowCells_.push_back( Vec3<int>( x, y, z ) );
          inflowVels_.push_back( Vec3<T>( u_x, u_y, u_z ) );
        }
      }
    }

  }

  bit = block.findAll( "outflow" );

  for ( ConfBlock::childIter it = bit.first; it != bit.second; ++it ) {

    ConfBlock& b = it->second;

    int xStart = ( x == -1 ) ? b.getParam<int>( "xStart" ) : x;
    int xEnd   = ( x == -1 ) ? b.getParam<int>( "xEnd" )   : x;
    int yStart = ( y == -1 ) ? b.getParam<int>( "yStart" ) : y;
    int yEnd   = ( y == -1 ) ? b.getParam<int>( "yEnd" )   : y;
    int zStart = ( z == -1 ) ? b.getParam<int>( "zStart" ) : z;
    int zEnd   = ( z == -1 ) ? b.getParam<int>( "zEnd" )   : z;
    int xDir   = ( x == -1 ) ? 0 : ( ( x > 1 ) ? -1 : 1 );
    int yDir   = ( y == -1 ) ? 0 : ( ( y > 1 ) ? -1 : 1 );
    int zDir   = ( z == -1 ) ? 0 : ( ( z > 1 ) ? -1 : 1 );
    assert( xStart > 0 && xEnd < flag_.getSizeX() - 1 &&
            yStart > 0 && yEnd < flag_.getSizeY() - 1 &&
            zStart > 0 && zEnd < flag_.getSizeZ() - 1 );
    for ( z = zStart; z <= zEnd; ++z ) {
      for ( y = yStart; y <= yEnd; ++y ) {
        for ( x = xStart; x <= xEnd; ++x ) {
          assert( flag_( x, y, z ) == UNDEFINED );
          flag_( x, y, z ) = OUTFLOW;
          outflowCells_.push_back( Vec3<int>( x, y, z ) );
          outflowDirs_.push_back( Vec3<int>( xDir, yDir, zDir ) );
        }
      }
    }

  }
}

template<typename T>
inline T LBM<T>::getTime( timeval &start, timeval &end ) {
  return (T) ( end.tv_sec - start.tv_sec )
          + (T) ( end.tv_usec - start.tv_usec ) / 1000000.;
}

template<typename T>
inline void LBM<T>::collideStream( int x, int y, int z ) {

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
  T omegai = 1 - omega_;
  // treat center value specially
  (*grid1_)( x, y, z, 0 ) = omegai * (*grid0_)( x, y, z, 0 )
                        + omega_ *  w[0] * fc;
  // loop over all velocity directions but center
  for ( int f = 1; f < Dim; ++f ) {
    T eiu = ex[f] * ux + ey[f] * uy + ez[f] * uz;
    (*grid1_)( x + ex[f], y + ey[f], z + ez[f], f )
      =   omegai * (*grid0_)( x, y, z, f )
        + omega_  * w[f] * ( fc +  3 * eiu + 4.5 * eiu * eiu);
  }
}

template<typename T>
inline void LBM<T>::collideStreamSmagorinsky( int x, int y, int z ) {

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
  T s = ( sqrt( nu_ * nu_ + 18. * cSmagoSqr_ * qo ) - nu_ ) / ( 6. * cSmagoSqr_);
  // Calculate turbulence modified inverse lattice viscosity
  T omega = 1. / ( 3. * ( nu_ + cSmagoSqr_ * s ) + .5 );
  T omegai = 1. - omega;

  // Loop over all velocity directions and stream collided distribution value
  // to neighboring cells
  for ( int f = 0; f < Dim; ++f ) {
    (*grid1_)( x + ex[f], y + ey[f], z + ez[f], f )
      =   omegai * (*grid0_)( x, y, z, f ) + omega * feq[f];
  }
}

template<typename T>
inline void LBM<T>::treatNoslip() {

  // Iterate over all no-slip boundary cells
  std::vector< Vec3<int> >::iterator iter;
  for( iter = noslipCells_.begin(); iter != noslipCells_.end(); iter++ ) {

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
inline void LBM<T>::treatVelocity() {

  // Iterate over all velocity boundary cells
  for( int i = 0; i < (int) velocityCells_.size(); ++i ) {

    // Fetch coordinates of current boundary cell
    int x = velocityCells_[i][0];
    int y = velocityCells_[i][1];
    int z = velocityCells_[i][2];
    // Fetch velocity of moving wall
    T ux = velocityVels_[i][0];
    T uy = velocityVels_[i][1];
    T uz = velocityVels_[i][2];
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

template<typename T>
inline void LBM<T>::treatInflow() {

  // Iterate over all inflow boundary cells
  for( int i = 0; i < (int) inflowCells_.size(); ++i ) {

    // Fetch coordinates of current boundary cell
    int x = inflowCells_[i][0];
    int y = inflowCells_[i][1];
    int z = inflowCells_[i][2];
    // Fetch inflow velocity
    T ux = inflowVels_[i][0];
    T uy = inflowVels_[i][1];
    T uz = inflowVels_[i][2];
    // Set velocity of this cell
    u_( x, y, z, 0 ) = ux;
    u_( x, y, z, 1 ) = uy;
    u_( x, y, z, 2 ) = uz;

    // Calculate equilibrium distribution functions with fixed density of 1.0
    // and set as distribution values of the inflow cell
    T fc = 1. - 1.5 * ( ux * ux + uy * uy + uz * uz );
    (*grid1_)( x, y, z, 0 )  = (1./3.)  *   fc; // C
    (*grid1_)( x, y, z, 1 )  = (1./18.) * ( fc + 3 *   uy        + 4.5 *   uy        *   uy ); // N
    (*grid1_)( x, y, z, 2 )  = (1./18.) * ( fc + 3 *   ux        + 4.5 *   ux        *   ux ); // E
    (*grid1_)( x, y, z, 3 )  = (1./18.) * ( fc - 3 *   uy        + 4.5 *   uy        *   uy ); // S
    (*grid1_)( x, y, z, 4 )  = (1./18.) * ( fc - 3 *   ux        + 4.5 *   ux        *   ux ); // W
    (*grid1_)( x, y, z, 5 )  = (1./18.) * ( fc + 3 *   uz        + 4.5 *   uz        *   uz ); // T
    (*grid1_)( x, y, z, 6 )  = (1./18.) * ( fc - 3 *   uz        + 4.5 *   uz        *   uz ); // B
    (*grid1_)( x, y, z, 7 )  = (1./36.) * ( fc + 3 * ( ux + uy ) + 4.5 * ( ux + uy ) * ( ux + uy ) ); // NE
    (*grid1_)( x, y, z, 8 )  = (1./36.) * ( fc + 3 * ( ux - uy ) + 4.5 * ( ux - uy ) * ( ux - uy ) ); // SE
    (*grid1_)( x, y, z, 9 )  = (1./36.) * ( fc - 3 * ( ux + uy ) + 4.5 * ( ux + uy ) * ( ux + uy ) ); // SW
    (*grid1_)( x, y, z, 10 ) = (1./36.) * ( fc - 3 * ( ux - uy ) + 4.5 * ( ux - uy ) * ( ux - uy ) ); // NW
    (*grid1_)( x, y, z, 11 ) = (1./36.) * ( fc + 3 * ( uy + uz ) + 4.5 * ( uy + uz ) * ( uy + uz ) ); // TN
    (*grid1_)( x, y, z, 12 ) = (1./36.) * ( fc + 3 * ( ux + uz ) + 4.5 * ( ux + uz ) * ( ux + uz ) ); // TE
    (*grid1_)( x, y, z, 13 ) = (1./36.) * ( fc - 3 * ( uy - uz ) + 4.5 * ( uy - uz ) * ( uy - uz ) ); // TS
    (*grid1_)( x, y, z, 14 ) = (1./36.) * ( fc - 3 * ( ux - uz ) + 4.5 * ( ux - uz ) * ( ux - uz ) ); // TW
    (*grid1_)( x, y, z, 15 ) = (1./36.) * ( fc + 3 * ( uy - uz ) + 4.5 * ( uy - uz ) * ( uy - uz ) ); // BN
    (*grid1_)( x, y, z, 16 ) = (1./36.) * ( fc + 3 * ( ux - uz ) + 4.5 * ( ux - uz ) * ( ux - uz ) ); // BE
    (*grid1_)( x, y, z, 17 ) = (1./36.) * ( fc - 3 * ( uy + uz ) + 4.5 * ( uy + uz ) * ( uy + uz ) ); // BS
    (*grid1_)( x, y, z, 18 ) = (1./36.) * ( fc - 3 * ( ux + uz ) + 4.5 * ( ux + uz ) * ( ux + uz ) ); // BW
  }
}

template<typename T>
inline void LBM<T>::treatOutflow() {

  // Iterate over all outflow boundary cells
  for( int i = 0; i < (int) outflowCells_.size(); ++i ) {

    // Fetch coordinates of current boundary cell
    int x = outflowCells_[i][0];
    int y = outflowCells_[i][1];
    int z = outflowCells_[i][2];
    // Fetch outflow direction
    int xDir = outflowDirs_[i][0];
    int yDir = outflowDirs_[i][1];
    int zDir = outflowDirs_[i][2];

    // Go over all distribution values, and copy the distribution values from
    // the neighboring cell in outflow direction
    for ( int f = 1; f < Dim; ++f ) {
      (*grid1_)( x - xDir, y - yDir, z - zDir, f ) = (*grid1_)( x, y, z, f );
    }
  }
}

template<typename T>
inline void LBM<T>::treatPressure() {

  // Iterate over all pressure boundary cells
  std::vector< Vec3<int> >::iterator iter;
  for( iter = pressureCells_.begin(); iter != pressureCells_.end(); iter++ ) {

    // Fetch coordinates of current boundary cell
    int x = (*iter)[0];
    int y = (*iter)[1];
    int z = (*iter)[2];
    // Fetch velocity of current boundary cell
    T ux = u_( x, y, z, 0 );
    T uy = u_( x, y, z, 1 );
    T uz = u_( x, y, z, 2 );

    // Calculate pressure corrected equilibrium distribution functions for
    // atmospheric pressure and set as distribution values of the pressure cell
    T fc = 2. - 3. * ( ux * ux + uy * uy + uz * uz );
    (*grid1_)( x, y, z, 0 )  = (1./3.)  *   fc                                    - (*grid1_)( x, y, z, finv[0] ); // C
    (*grid1_)( x, y, z, 1 )  = (1./18.) * ( fc + 9. *   uy        *   uy )        - (*grid1_)( x, y, z, finv[1] ); // N
    (*grid1_)( x, y, z, 2 )  = (1./18.) * ( fc + 9. *   ux        *   ux )        - (*grid1_)( x, y, z, finv[2] ); // E
    (*grid1_)( x, y, z, 3 )  = (1./18.) * ( fc + 9. *   uy        *   uy )        - (*grid1_)( x, y, z, finv[3] ); // S
    (*grid1_)( x, y, z, 4 )  = (1./18.) * ( fc + 9. *   ux        *   ux )        - (*grid1_)( x, y, z, finv[4] ); // W
    (*grid1_)( x, y, z, 5 )  = (1./18.) * ( fc + 9. *   uz        *   uz )        - (*grid1_)( x, y, z, finv[5] ); // T
    (*grid1_)( x, y, z, 6 )  = (1./18.) * ( fc + 9. *   uz        *   uz )        - (*grid1_)( x, y, z, finv[6] ); // B
    (*grid1_)( x, y, z, 7 )  = (1./36.) * ( fc + 9. * ( ux + uy ) * ( ux + uy ) ) - (*grid1_)( x, y, z, finv[7] ); // NE
    (*grid1_)( x, y, z, 8 )  = (1./36.) * ( fc + 9. * ( ux - uy ) * ( ux - uy ) ) - (*grid1_)( x, y, z, finv[8] ); // SE
    (*grid1_)( x, y, z, 9 )  = (1./36.) * ( fc + 9. * ( ux + uy ) * ( ux + uy ) ) - (*grid1_)( x, y, z, finv[9] ); // SW
    (*grid1_)( x, y, z, 10 ) = (1./36.) * ( fc + 9. * ( ux - uy ) * ( ux - uy ) ) - (*grid1_)( x, y, z, finv[10] ); // NW
    (*grid1_)( x, y, z, 11 ) = (1./36.) * ( fc + 9. * ( uy + uz ) * ( uy + uz ) ) - (*grid1_)( x, y, z, finv[11] ); // TN
    (*grid1_)( x, y, z, 12 ) = (1./36.) * ( fc + 9. * ( ux + uz ) * ( ux + uz ) ) - (*grid1_)( x, y, z, finv[12] ); // TE
    (*grid1_)( x, y, z, 13 ) = (1./36.) * ( fc + 9. * ( uy - uz ) * ( uy - uz ) ) - (*grid1_)( x, y, z, finv[13] ); // TS
    (*grid1_)( x, y, z, 14 ) = (1./36.) * ( fc + 9. * ( ux - uz ) * ( ux - uz ) ) - (*grid1_)( x, y, z, finv[14] ); // TW
    (*grid1_)( x, y, z, 15 ) = (1./36.) * ( fc + 9. * ( uy - uz ) * ( uy - uz ) ) - (*grid1_)( x, y, z, finv[15] ); // BN
    (*grid1_)( x, y, z, 16 ) = (1./36.) * ( fc + 9. * ( ux - uz ) * ( ux - uz ) ) - (*grid1_)( x, y, z, finv[16] ); // BE
    (*grid1_)( x, y, z, 17 ) = (1./36.) * ( fc + 9. * ( uy + uz ) * ( uy + uz ) ) - (*grid1_)( x, y, z, finv[17] ); // BS
    (*grid1_)( x, y, z, 18 ) = (1./36.) * ( fc + 9. * ( ux + uz ) * ( ux + uz ) ) - (*grid1_)( x, y, z, finv[18] ); // BW
  }
}

template<>
void LBM<double>::writeVtkFile() {

  // Open file for writing
  std::ostringstream oss;
  oss << vtkFileName_ << "." << curStep_ << ".vtk";
  std::cout << "Writing file '" << oss.str() << "' for time step " << curStep_ << std::endl;
  std::ofstream vtkFile( oss.str().c_str(), std::ios::binary | std::ios::out );

  // Get size of domain without ghost layers
  int sizeX = grid0_->getSizeX() - 2;
  int sizeY = grid0_->getSizeY() - 2;
  int sizeZ = grid0_->getSizeZ() - 2;

  // Write file header
  vtkFile << "# vtk DataFile Version 2.0\n";
  vtkFile << "VTK output file for time step " << curStep_ << "\n\n";
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
void LBM<float>::writeVtkFile() {

  // Open file for writing
  std::ostringstream oss;
  oss << vtkFileName_ << "." << curStep_ << ".vtk";
  std::cout << "Writing file '" << oss.str() << "' for time step " << curStep_ << std::endl;
  std::ofstream vtkFile( oss.str().c_str(), std::ios::binary | std::ios::out );

  // Get size of domain without ghost layers
  int sizeX = grid0_->getSizeX() - 2;
  int sizeY = grid0_->getSizeY() - 2;
  int sizeZ = grid0_->getSizeZ() - 2;

  // Write file header
  vtkFile << "# vtk DataFile Version 2.0\n";
  vtkFile << "VTK output file for time step " << curStep_ << "\n\n";
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
