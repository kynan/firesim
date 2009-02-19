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
  PRESSURE  = 5
};

//! Lattice Boltzmann Method fluid solver

//! Uses the BGK collision model and Smagorinsky turbulence correction
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

  void setup( std::string configFileName );

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

  //! Perform a collide-stream step without turbulence correction

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
