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

template<typename T>
struct Sphere {

  Sphere( T _x, T _y, T _z, T _r, T _u_x, T _u_y, T _u_z ) :
      x(_x), y(_y), z(_z), r(_r), u_x(_u_x), u_y(_u_y), u_z(_u_z) {}

  void move() {
    x += u_x;
    y += u_y;
    z += u_z;
  }

  T x;
  T y;
  T z;
  T r;
  T u_x;
  T u_y;
  T u_z;
};

template<typename T>
struct PerformanceData {
  PerformanceData() : scTime( 0. ), nTime( 0. ), vTime( 0. ), iTime( 0. ), oTime( 0. ), pTime( 0. ), cTime( 0. ), sTime( 0. ), mTime( 0. ) {}
  std::string fileName;
  T scTime;
  T nTime;
  T vTime;
  T iTime;
  T oTime;
  T pTime;
  T cTime;
  T sTime;
  T mTime;
};

//! Enum describing possible states of a cell

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

//! Uses the BGK collision model and Smagorinsky turbulence correction. The
//! coordinate system used is left-handed, i.e. screen coordinates
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
//!                  each side of the domain (\em bottom, \em top, \em back,
//!                  \em front, \em left, \em right) which in turn can contain
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
//!  specified (i.e. \e y for \e bottom and \e top, \e z for \e back and \e
//!  front and \e x for \e left and \e right.
//!  - \b xStart \e int Start cell in x-direction (not required for \e left and \e right)
//!  - \b xEnd   \e int End cell in x-direction (not required for \e left and \e right)
//!  - \b yStart \e int Start cell in y-direction (not required for \e bottom and \e top)
//!  - \b yEnd   \e int End cell in y-direction (not required for \e bottom and \e top)
//!  - \b zStart \e int Start cell in z-direction (not required for \e back and \e front)
//!  - \b zEnd   \e int End cell in z-direction (not required for \e back and \e front)
//!  .
//!  Boundary cells that are not specified are implicitly set to \e noslip. An
//!  error is thrown if a boundary cell is specified twice. \e bottom and \e top
//!  have priority over \e back and \e front which in turn have priority over
//!  \e left and \e right, i.e. the edges of the rectangle are added to the
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
//!       zStart  1;
//!       zEnd   60;
//!     }
//!   }
//!   bottom {
//!     inflow {
//!       xStart 26;
//!       xEnd 35;
//!       zStart 26;
//!       zEnd 35;
//!       u_x 0.0;
//!       u_y 0.0;
//!       u_z 0.1;
//!     }
//!   }
//!   left {
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
//!   right {
//!     outflow {
//!       zStart 2;
//!       zEnd 59;
//!       yStart 2;
//!       yEnd 59;
//!     }
//!   }
//!   back {
//!     outflow {
//!       xStart 1;
//!       xEnd 60;
//!       yStart 2;
//!       yEnd 59;
//!     }
//!   }
//!   front {
//!     outflow {
//!       xStart 1;
//!       xEnd 60;
//!       yStart 2;
//!       yEnd 59;
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

  //! Default constructor

  //! Does no initialization.

  LBM() : curStep_( 0 ) {}

  //! Constructor

  //! Initializes the geometry of the domain and the boundaries according to the
  //! configuration file given.
  //! \note There is no default constructor available!
  //! \param[in] configFileName Path to configuration file

  LBM( const std::string configFileName );

  //! Constructor

  //! Initializes the geometry of the domain and the boundaries according to the
  //! configuration block given.
  //! \note There is no default constructor available!
  //! \param[in] base Root configuration block

  LBM( ConfBlock& base );

  //! Destructor

  //! Deletes the grids of the distribution functions

  virtual ~LBM();

  //! Run the solver

  //! The output will be written to a series of VTK files (legacy format)

  void run();

  //! Perform a single LBM step

  //! The output will be written to a VTK file if appropriate (legacy format)

  double runStep();

  //! Get tri-linearly interpolated velocity at given position

  //! \param[in] x x-coordinate of position to get velocity
  //! \param[in] y y-coordinate of position to get velocity
  //! \param[in] z z-coordinate of position to get velocity

  inline Vec3<T> getVelocity( T x, T y, T z );

  //! Get tri-linearly interpolated velocity at given position

  //! \param[in] v 3-component vector (x,y,z) of position to get velocity

  inline Vec3<T> getVelocity( const Vec3<T>& v ) {
    return getVelocity( v[0], v[1], v[2] );
  }

  //! Setup the solver by processing configuration file

  //! Initializes the geometry of the domain and the boundaries according to the
  //! configuration given.
  //! \param[in] base Root block of parsed configuration file

  void setup( ConfBlock& base );

protected:

  // ========================= //
  // Internal helper functions //
  // ========================= //

  //! Process a boundary block

  //! \param[in] block Boundary block to process (either \e bottom, \e top, \e
  //!                  back, \e front, \e right, \e left)
  //! \param[in] x     Either fixed value for x-coordinate or -1 to read range
  //!                  from configuration
  //! \param[in] y     Either fixed value for y-coordinate or -1 to read range
  //!                  from configuration
  //! \param[in] z     Either fixed value for z-coordinate or -1 to read range
  //!                  from configuration

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

  inline void collideStream( int x, int y, int z );

#ifndef NSMAGO
  //! Perform a collide-stream step with turbulence correction

  //! \param[in] x Cell coordinate for dimension x
  //! \param[in] y Cell coordinate for dimension y
  //! \param[in] z Cell coordinate for dimension z

  inline void collideStreamSmagorinsky( int x, int y, int z );
#endif

  //! Treat the no-slip boundary cells

  inline void treatNoslip();

  inline void treatStaircase();

  //! Treat the boundary cells with fixed velocity

  inline void treatVelocity();

  //! Treat the inflow boundary cells

  inline void treatInflow();

  //! Treat the outflow boundary cells

  inline void treatOutflow();

  //! Treat the outflow boundary cells with fixed athmospheric pressure

  inline void treatPressure();

  //! Treat the curved boundary cells

  inline void treatCurved();

  inline void moveSphere();

  //! Write out the VTK file for a given timestep

  //! \note Output is in binary VTK legacy file format

  void writeVtkFile();

  void writePerformanceSummary();

  std::string calcMLUP( T time, int cells );

  // ============ //
  // Data members //
  // ============ //

  //! Inverse lattice velocity

  T omega_;

#ifndef NSMAGO
  //! Lattice viscosity

  T nu_;

  //! Squared Smagorinsky turbulence constant

  T cSmagoSqr_;
#endif

  //! Number of collide-stream steps to perform

  int maxSteps_;

  //! Current time step

  int curStep_;

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

  //! Vector holding all sphere obstacles

  std::vector< Sphere<T> > sphereObstacles_;

  //! List with coordinates of all noslip cells

  std::vector< Vec3<int> > noslipCells_;
  std::vector< Vec3<int> > staircaseCells_;

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

  std::vector< int > outflowDFs_;

  //! List with coordinates of all pressure cells

  std::vector< Vec3<int> > pressureCells_;

  //! List with outflow directions for all pressure cells

  std::vector< int > pressureDFs_;

  //! List with coordinates of all curved boundary cells

  std::vector< Vec3<int> > curvedCells_;

  //! List with fluid fractions for all lattice links of curved boundary cells

  std::vector< std::vector<T> > curvedDeltas_;

  PerformanceData<T> prof_;

};

//! Specialization of the the VTK file writer for template type double
//! \note This is necessary as the type needs to be written to the VTK file in
//!       ASCII and there is no way to get an ASCII representation of the \e
//!       type of a template paramter

template<>
void LBM<double>::writeVtkFile();

//! Specialization of the the VTK file writer for template type float
//! \note This is necessary as the type needs to be written to the VTK file in
//!       ASCII and there is no way to get an ASCII representation of the \e
//!       type of a template paramter

template<>
void LBM<float>::writeVtkFile();

} // namespace lbm

#endif /* LBM_H_ */
