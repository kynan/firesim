//! \file LBM.h
//! \date   Jan 2, 2009
//! \author Florian Rathgeber

#ifndef LBM_H_
#define LBM_H_

#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <string>
#include <endian.h>
#include <sys/time.h>

#include "D3Q19.h"
#include "../Vec.h"
#include "../Grid.h"
#include "../confparser/ConfParser.h"
#include "../confparser/ConfBlock.h"

using namespace confparser;

//! Common namespace for all LBM classes

namespace lbm {

enum Flag {
  UNDEFINED = 0,
  FLUID     = 1,
  NOSLIP    = 2,
  VELOCITY  = 3,
  INFLOW    = 4,
  OUTFLOW   = 5,
  PRESSURE  = 6
};

//! Lattice Boltzmann Method fluid solver

//! Uses the BGK collision model and Smagorinsky turbulence correction.
//!
//! The following boundary conditions are implemented:
//! - \em no-slip  Bounce back condition (fixed wall)
//! - \em velocity Bounce back condition (moving wall)
//! - \em inflow   Inflow velocity with density of 1.0
//! - \em outflow  Enforces zero-gradient pressure at the boundary
//! - \em pressure Enforces atmospheric pressure at the boundary
//!
//! The constructor requires a configuration file in a ConfParser compatible
//! format to be given, which must specify the following hierarchy of blocks:
//! - \em domain Specifies the domain size
//!  - \b sizeX \em int Number of cells in x-dimension
//!  - \b sizeY \em int Number of cells in y-dimension
//!  - \b sizeZ \em int Number of cells in z-dimension
//! - \em parameters Specifies the simulation parameters
//!  - \b omega        \em float Inverse lattice velocity
//!  - \b cSmagorinsky \em float Smagorinsky turbulence constant
//!  - \b maxSteps     \em int   Number of time steps to simulate
//! - \em vtk Specifies vtk output (optional, no output will be generated if
//!           this block is not present)
//!  - \b vtkFileName \em string Base name and path for the vtk file to write out
//!                           (relative to location of executable, will be
//!                           appended by current timestep and .vtk extension)
//!  - \b vtkStep     \em int    Number of timesteps in between 2 vtk file outputs
//! - \em boundaries Specifies the boundary conditions, contains a subblock for
//!                  each side of the domain (\em bottom, \em top, \em north,
//!                  \em south, \em east, \em west) which in turn can contain
//!                  1 or more of the following subblocks to specify the type of
//!                  boundary condition.
//!  - \em noslip   Bounce back condition (fixed wall)
//!  - \em velocity Bounce back condition (moving wall)
//!   - \b u_x \e float Velocity of moving wall in x-direction
//!   - \b u_y \e float Velocity of moving wall in y-direction
//!   - \b u_z \e float Velocity of moving wall in y-direction
//!  - \em inflow   Inflow velocity with density of 1.0
//!   - \b u_x \e float Inflow velocity in x-direction
//!   - \b u_y \e float Inflow velocity in y-direction
//!   - \b u_z \e float Inflow velocity in z-direction
//!  - \em outflow  Enforces zero-gradient pressure at the boundary
//!  - \em pressure Enforces atmospheric pressure at the boundary
//!  .
//!  Every boundary block needs to specify the coordinates of the rectangular
//!  patch, where the fixed coordinate is implicitly set and needs not to be
//!  specified (i.e. \e z for \e bottom and \e top, \e y for \e north and \e
//!  south and \e x for \e east and \e west.
//!  - \b xStart \e int Start cell in x-direction (not required for \e east and \e west)
//!  - \b xEnd   \e int End cell in x-direction (not required for \e east and \e west)
//!  - \b yStart \e int Start cell in y-direction (not required for \e north and \e south)
//!  - \b yEnd   \e int End cell in y-direction (not required for \e north and \e south)
//!  - \b zStart \e int Start cell in z-direction (not required for \e bottom and \e top)
//!  - \b zEnd   \e int End cell in z-direction (not required for \e bottom and \e top)
//!  .
//!  Boundary cells that are not specified are implicitly set to \e noslip. An
//!  error is thrown if a boundary cell is specified twice. \e bottom and \e top
//!  have priority over \e north and \e south which in turn have priority over
//!  \e east and \e west, i.e. the edges of the rectangle are added to the
//!  boundary regions according to these priorities given.
//!
//! Example configuration file:
//! \code
//! domain {
//!   sizeX 62;
//!   sizeY 62;
//!   sizeZ 62;
//! }
//!
//! parameters {
//!   omega 1.9;
//!   cSmagorinsky 0.03;
//!   maxSteps 1201;
//! }
//!
//! vtk {
//!   vtkFileName vtk/outflow_60_60_60_1.9_0.03_1201_25;
//!   vtkStep 25;
//! }
//!
//! boundaries {
//!   top {
//!     outflow {
//!       xStart  1;
//!       xEnd   60;
//!       yStart  1;
//!       yEnd   60;
//!     }
//!   }
//!   bottom {
//!     inflow {
//!       xStart 26;
//!       xEnd 35;
//!       yStart 26;
//!       yEnd 35;
//!       u_x 0.0;
//!       u_y 0.0;
//!       u_z 0.1;
//!     }
//!   }
//!   west {
//!     inflow {
//!       zStart 2;
//!       zEnd 59;
//!       yStart 2;
//!       yEnd 59;
//!       u_x 0.001;
//!       u_y 0.0;
//!       u_z 0.0;
//!     }
//!   }
//!   east {
//!     outflow {
//!       zStart 2;
//!       zEnd 59;
//!       yStart 2;
//!       yEnd 59;
//!     }
//!   }
//!   north {
//!     outflow {
//!       xStart 1;
//!       xEnd 60;
//!       zStart 2;
//!       zEnd 59;
//!     }
//!   }
//!   south {
//!     outflow {
//!       xStart 1;
//!       xEnd 60;
//!       zStart 2;
//!       zEnd 59;
//!     }
//!   }
//! }
//! \endcode
//! \note define the preprocessor symbol NSMAGO to disable turbulence correction

template<typename T>
class LBM {

public:

  // ============================ //
  // Constructors and destructors //
  // ============================ //

  //! Constructor

  //! Initializes the geometry of the domain and the boundaries according to the
  //! configuration file given.
  //! \note There is no default constructor available!
  //! \param[in] configFileName Path to configuration file

  LBM( std::string configFileName );

  //! Destructor

  //! Deletes the grids of the distribution functions

  virtual ~LBM();

  //! Run the solver

  //! The output will be written to a series of VTK files (legacy format)

  void run();

protected:

  // ========================= //
  // Internal helper functions //
  // ========================= //

  //! Setup the solver by processing configuration file

  //! Initializes the geometry of the domain and the boundaries according to the
  //! configuration file given.
  //! \param[in] configFileName Path to configuration file

  void setup( std::string configFileName );

  //! Process a boundary block

  //! \param[in] block Boundary block to process (either \e bottom, \e top, \e
  //!                  north, \e south, \e east, \e west)

  inline void setupBoundary( ConfBlock& block, int x, int y, int z );

  //! Get the time difference between two measurements

  //! \param[in] start Start time of measurement as return by \e gettimeofday
  //! \param[in] end   End time of measurement as return by \e gettimeofday
  //! \returns Time difference \e end - \e start in seconds

  inline T getTime( timeval &start, timeval &end );

  //! Perform a collide-stream step without turbulence correction

  //! \param[in] x Cell coordinate for dimension x
  //! \param[in] y Cell coordinate for dimension y
  //! \param[in] z Cell coordinate for dimension z
  //! \param[in] omega Inverse lattice velocity

  inline void collideStream( int x, int y, int z, T omega );

  //! Perform a collide-stream step with turbulence correction

  //! \param[in] x Cell coordinate for dimension x
  //! \param[in] y Cell coordinate for dimension y
  //! \param[in] z Cell coordinate for dimension z
  //! \param[in] nu Lattice viscosity
  //! \param[in] cSqr Squared Smagorinsky constant

  inline void collideStreamSmagorinsky( int x, int y, int z, T nu, T cSqr );

  //! Treat the no-slip boundary cells

  inline void treatNoslip();

  //! Treat the boundary cells with fixed velocity

  inline void treatVelocity();

  //! Treat the inflow boundary cells

  inline void treatInflow();

  //! Treat the outflow boundary cells

  inline void treatOutflow();

  //! Treat the outflow boundary cells with fixed athmospheric pressure

  inline void treatPressure();

  //! Write out the VTK file for a given timestep

  //! \note Output is in binary VTK legacy file format
  //! \param[in] timestep Current simulation step

  void writeVtkFile( int timestep );

  // ============ //
  // Data members //
  // ============ //

  //! Inverse lattice velocity

  T omega_;

  //! Smagorinsky turbulence constant

  T cSmagorinsky_;

  //! Number of collide-stream steps to perform

  int maxSteps_;

  //! VTK files will be written out at multiples of this step size

  int vtkStep_;

  //! Path and base name of VTK files to write out. Will be appended by
  //! \em .currentStep.vtk

  std::string vtkFileName_;

  //! Distribution function field (switched with grid1_ after each time step)

  Grid<T,Dim>* grid0_;

  //! Distribution function field (switched with grid0_ after each time step)

  Grid<T,Dim>* grid1_;

  //! Velocity field

  Grid<T,3> u_;

  //! Density field

  Grid<T,1> rho_;

  //! Flag field

  Grid<Flag,1> flag_;

  //! List with coordinates of all noslip cells

  std::vector< Vec3<int> > noslipCells_;

  //! List with coordinates of all velocity cells

  std::vector< Vec3<int> > velocityCells_;

  //! List with velocities for all velocity cells

  std::vector< Vec3<T> > velocityVels_;

  //! List with coordinates of all inflow cells

  std::vector< Vec3<int> > inflowCells_;

  //! List with velocities for all inflow cells

  std::vector< Vec3<T> > inflowVels_;

  //! List with coordinates of all outflow cells

  std::vector< Vec3<int> > outflowCells_;

  //! List with outflow directions for all outflow cells

  std::vector< Vec3<int> > outflowDirs_;

  //! List with coordinates of all pressure cells

  std::vector< Vec3<int> > pressureCells_;

};

//! Specialization of the the VTK file writer for template type double
//! \note This is necessary as the type needs to be written to the VTK file in
//!       ASCII and there is no way to get an ASCII representation of the \e
//!       type of a template paramter

template<>
void LBM<double>::writeVtkFile( int timestep );

//! Specialization of the the VTK file writer for template type float
//! \note This is necessary as the type needs to be written to the VTK file in
//!       ASCII and there is no way to get an ASCII representation of the \e
//!       type of a template paramter

template<>
void LBM<float>::writeVtkFile( int timestep );

} // namespace lbm

// Include the definitions
#include "LBM_def.h"

#endif /* LBM_H_ */
