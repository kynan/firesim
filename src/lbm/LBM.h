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

//! Common namespace for all LBM classes

namespace lbm {

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

    //! Initializes the geometry of the domain.
    //! \note There is no default constructor available!
    //! \param[in] sizeX Number of cells in x dimension, including ghost layer
    //! \param[in] sizeY Number of cells in y dimension, including ghost layer
    //! \param[in] sizeZ Number of cells in z dimension, including ghost layer
    //! \param[in] boundaryCells Vector of coordinate triples specifying the
    //!                          location of the no-slip boundary
    //! \param[in] velocityCells Vector of coordinate triples specifying the
    //!                          location of fixed-velocity cells
    //! \param[in] velocities    Vector of 3-dimensional velocities corresponding
    //!                          to the velocity cells (vector sizes must match!)

    LBM( int sizeX,
         int sizeY,
         int sizeZ,
         std::vector< Vec3<int> > &boundaryCells,
         std::vector< Vec3<int> > &velocityCells,
         std::vector< Vec3<T> > &velocities );

    //! Destructor

    //! Deletes the grids of the distribution functions

    virtual ~LBM();

    //! Run the solver with a given parameter set

    //! The output will be written to a series of VTK files (legacy format)
    //! \param[in] omega Inverse lattice velocity
    //! \param[in] cSmagorinsky Smagorinsky constant for turbulence correction
    //! \param[in] maxSteps Number of time steps to simulate
    //! \param[in] vtkStep Number of time steps between two output file writes
    //! \param[in] vtkFileName Name of the VTK output files to write out. Should
    //!                        be prepended by the path where to put the files
    //!                        relative to the executable. Will be appended by
    //!                        \e .step.vtk where \e step is the current
    //!                        simulation time step

    void run( T omega,
              T cSmagorinsky,
               int maxSteps,
               int vtkStep,
               std::string vtkFileName );

  protected:

    // ========================= //
    // Internal helper functions //
    // ========================= //

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

    inline void treatBoundary();

    //! Treat the boundary cells with fixed velocity

    inline void treatVelocities();

    //! Write out the VTK file for a given timestep

    //! \note Output is in binary VTK legacy file format
    //! \param[in] timestep Current simulation step
    //! \param[in] vtkFileName Name of the VTK output file to write out. Must be
    //!                        prepended by the path where to put the files
    //!                        relative to the executable. Will be appended by
    //!                        \e .timestep.vtk

    void writeVtkFile( int timestep, std::string vtkFileName );

    // ============ //
    // Data members //
    // ============ //

    //! Distribution function field (switched with grid1_ after each time step)

    Grid<T,Dim>* grid0_;

    //! Distribution function field (switched with grid0_ after each time step)

    Grid<T,Dim>* grid1_;

    //! Velocity field

    Grid<T,3> u_;

    //! Density field

    Grid<T,1> rho_;

    //! List with coordinates of all boundary cells

    std::vector< Vec3<int> > boundaryCells_;

    //! List with coordinates of all velocity cells

    std::vector< Vec3<int> > velocityCells_;

    //! List with velocities for all velocity cells

    std::vector< Vec3<T> > velocities_;
  };

  //! Specialization of the the VTK file writer for template type double

  //! \note This is necessary as the type needs to be written to the VTK file in
  //!       ASCII and there is no way to get an ASCII representation of the \e
  //!       type of a template paramter

  template<>
  void LBM<double>::writeVtkFile( int timestep, std::string vtkFileName );

  //! Specialization of the the VTK file writer for template type float

  //! \note This is necessary as the type needs to be written to the VTK file in
  //!       ASCII and there is no way to get an ASCII representation of the \e
  //!       type of a template paramter

  template<>
  void LBM<float>::writeVtkFile( int timestep, std::string vtkFileName );

} // namespace lbm

// Include the definitions
#include "LBM_def.h"

#endif /* LBM_H_ */
