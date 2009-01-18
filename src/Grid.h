//! \file Grid.h
//! \date   Jan 2, 2009
//! \author Florian Rathgeber

#ifndef GRID_H_
#define GRID_H_

#include <assert.h>
#include <vector>

//! 3-dimensional grid of Cells with template specifyable number of cell
//! variables and type

//! \tparam T Type of cell variables
//! \tparam Cellsize Number of cell variables

template<typename T, int Cellsize>
class Grid {

public:

  // ============================ //
  // Constructors and destructors //
  // ============================ //

  //! Default constructor

  //! Creates a grid with zero cells in each dimension

  Grid();

  //! Constructor to specify grid dimensions

  //! \param[in] sizeX Grid size in x-dimension
  //! \param[in] sizeY Grid size in y-dimension
  //! \param[in] sizeZ Grid size in z-dimension

  Grid(int sizeX, int sizeY, int sizeZ);

  //! Constructor to specify grid dimensions and initialization value

  //! \param[in] sizeX Grid size in x-dimension
  //! \param[in] sizeY Grid size in y-dimension
  //! \param[in] sizeZ Grid size in z-dimension
  //! \param[in] val Initial value for all cell variables

  Grid(int sizeX, int sizeY, int sizeZ, T val);

  //! Destructor

  //! Clears the vector of cell variables

  virtual ~Grid();

  // ============= //
  // Set operators //
  // ============= //

  //! Initialize the grid with given dimensions and initial value

  //! \param[in] sizeX Grid size in x-dimension
  //! \param[in] sizeY Grid size in y-dimension
  //! \param[in] sizeZ Grid size in z-dimension
  //! \param[in] val Initial value for all cell variables

  void init( int sizeX, int sizeY, int sizeZ, T val );

  //! Initialize the grid with given initial value

  //! \param[in] val Initial value for all cell variables

  void init( T val );

  // ================ //
  // Access operators //
  // ================ //

  //! Get size in x-dimension

  //! \returns Grid size in x-dimension

  inline int getSizeX() {
    return sizeX_;
  }

  //! Get size in y-dimension

  //! \returns Grid size in y-dimension

  inline int getSizeY() {
    return sizeY_;
  }

  //! Get size in z-dimension

  //! \returns Grid size in z-dimension

  inline int getSizeZ() {
    return sizeZ_;
  }

  //! Get specified cell variable of specified cell

  //! \param[in] x Cell coordinate in x-dimension
  //! \param[in] y Cell coordinate in y-dimension
  //! \param[in] z Cell coordinate in z-dimension
  //! \param[in] f Cell variable index
  //! \return Reference to specified cell variable of specified cell

  inline T& operator()(int x, int y, int z, int f);

  //! Get specified const cell variable of specified cell

  //! \param[in] x Cell coordinate in x-dimension
  //! \param[in] y Cell coordinate in y-dimension
  //! \param[in] z Cell coordinate in z-dimension
  //! \param[in] f Cell variable index
  //! \return Const reference to specified cell variable of specified cell

  inline const T& operator()(int x, int y, int z, int f) const;

protected:

  // ============ //
  // Data members //
  // ============ //

  //! Number of cells in x-dimension

  int sizeX_;

  //! Number of cells in y-dimension

  int sizeY_;

  //! Number of cells in z-dimension

  int sizeZ_;

  //! Linearized, 1-dimensional representation of the 3D data grid

  std::vector<T> data_;

};


//! 3-dimensional grid of cells with template type, specialization to a single
//! cell variable

//! \tparam T Type of cell variables

template<typename T>
class Grid<T,1> {

public:

  // ============================ //
  // Constructors and destructors //
  // ============================ //

  //! Default constructor

  //! Creates a grid with zero cells in each dimension

  Grid();

  //! Constructor to specify grid dimensions

  //! \param[in] sizeX Grid size in x-dimension
  //! \param[in] sizeY Grid size in y-dimension
  //! \param[in] sizeZ Grid size in z-dimension

  Grid(int sizeX, int sizeY, int sizeZ);

  //! Constructor to specify grid dimensions and initialization value

  //! \param[in] sizeX Grid size in x-dimension
  //! \param[in] sizeY Grid size in y-dimension
  //! \param[in] sizeZ Grid size in z-dimension
  //! \param[in] val Initial value for cells

  Grid(int sizeX, int sizeY, int sizeZ, T val);

  //! Destructor

  //! Clears the vector of cell variables

  virtual ~Grid();

  // ============= //
  // Set operators //
  // ============= //

  //! Initialize the grid with given dimensions and initial value

  //! \param[in] sizeX Grid size in x-dimension
  //! \param[in] sizeY Grid size in y-dimension
  //! \param[in] sizeZ Grid size in z-dimension
  //! \param[in] val Initial value for cells

  void init( int sizeX, int sizeY, int sizeZ, T val );

  //! Initialize the grid with given initial value

  //! \param[in] val Initial value for cells

  void init( T val );

  // ================ //
  // Access operators //
  // ================ //

  //! Get size in x-dimension

  //! \returns Grid size in x-dimension

  inline int getSizeX() {
    return sizeX_;
  }

  //! Get size in y-dimension

  //! \returns Grid size in y-dimension

  inline int getSizeY() {
    return sizeY_;
  }

  //! Get size in z-dimension

  //! \returns Grid size in z-dimension

  inline int getSizeZ() {
    return sizeZ_;
  }

  //! Get cell variable of specified cell

  //! \param[in] x Cell coordinate in x-dimension
  //! \param[in] y Cell coordinate in y-dimension
  //! \param[in] z Cell coordinate in z-dimension
  //! \return Reference to specified cell variable of specified cell

  inline T& operator()(int x, int y, int z);

  //! Get specified constant cell variable of specified cell

  //! \param[in] x Cell coordinate in x-dimension
  //! \param[in] y Cell coordinate in y-dimension
  //! \param[in] z Cell coordinate in z-dimension
  //! \return Constant reference to specified cell variable of specified cell

  inline const T& operator()(int x, int y, int z) const;

protected:

  //! Number of cells in x-dimension

  int sizeX_;

  //! Number of cells in y-dimension

  int sizeY_;

  //! Number of cells in z-dimension

  int sizeZ_;

  //! Linearized, 1-dimensional representation of the 3D data grid

  std::vector<T> data_;

};

// include the definitions
#include "Grid_def.h"

#endif /* GRID_H_ */
