//! \file Grid.h
//! \date   Jan 2, 2009
//! \author Florian Rathgeber

#ifndef GRID_H_
#define GRID_H_

#include <assert.h>
#include <vector>

//! Common namespace for all LBM classes

namespace lbm {

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

  Grid() :
    sizeX_(0), sizeY_(0), sizeZ_(0), data_(0) {
  }

  //! Constructor to specify grid dimensions

  //! \param[in] sizeX Grid size in x-dimension
  //! \param[in] sizeY Grid size in y-dimension
  //! \param[in] sizeZ Grid size in z-dimension

  Grid(int sizeX, int sizeY, int sizeZ) :
    sizeX_(sizeX),
    sizeY_(sizeY),
    sizeZ_(sizeZ),
    data_(Cellsize * sizeX * sizeY * sizeZ) {
  }

  //! Constructor to specify grid dimensions and initialization value

  //! \param[in] sizeX Grid size in x-dimension
  //! \param[in] sizeY Grid size in y-dimension
  //! \param[in] sizeZ Grid size in z-dimension
  //! \param[in] val Initial value for all cell variables

  Grid(int sizeX, int sizeY, int sizeZ, T val) :
    sizeX_(sizeX),
    sizeY_(sizeY),
    sizeZ_(sizeZ),
    data_(Cellsize * sizeX * sizeY * sizeZ, val) {
  }

  //! Destructor

  //! Clears the vector of cell variables

  virtual ~Grid() {
    data_.clear();
  }

  // ============= //
  // Set operators //
  // ============= //

  //! Initialize the grid with given dimensions and initial value

  //! \param[in] sizeX Grid size in x-dimension
  //! \param[in] sizeY Grid size in y-dimension
  //! \param[in] sizeZ Grid size in z-dimension
  //! \param[in] val Initial value for all cell variables

  void init( int sizeX, int sizeY, int sizeZ, T val ) {
    data_.assign( Cellsize * sizeX * sizeY * sizeZ, val );
  }

  //! Initialize the grid with given initial value

  //! \param[in] val Initial value for all cell variables

  void init( T val ) {
    data_.assign( Cellsize * sizeX_ * sizeY_ * sizeZ_, val );
  }

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

  inline T& operator()(int x, int y, int z, int f) {
    assert( x >= 0 && x < sizeX_ );
    assert( y >= 0 && y < sizeY_ );
    assert( z >= 0 && z < sizeZ_ );
    assert( f >= 0 && f < Cellsize );

    return data_[((z * sizeY_ + y) * sizeX_ + x) * Cellsize + f];
  }

  //! Get specified const cell variable of specified cell

  //! \param[in] x Cell coordinate in x-dimension
  //! \param[in] y Cell coordinate in y-dimension
  //! \param[in] z Cell coordinate in z-dimension
  //! \param[in] f Cell variable index
  //! \return Const reference to specified cell variable of specified cell

  inline const T& operator()(int x, int y, int z, int f) const {
    assert( x >= 0 && x < sizeX_ );
    assert( y >= 0 && y < sizeY_ );
    assert( z >= 0 && z < sizeZ_ );
    assert( f >= 0 && f < Cellsize );

    return data_[((z * sizeY_ + y) * sizeX_ + x) * Cellsize + f];
  }

private:

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

  inline Grid() :
    sizeX_(0), sizeY_(0), sizeZ_(0), data_(0) {
  }

  //! Constructor to specify grid dimensions

  //! \param[in] sizeX Grid size in x-dimension
  //! \param[in] sizeY Grid size in y-dimension
  //! \param[in] sizeZ Grid size in z-dimension

  inline Grid(int sizeX, int sizeY, int sizeZ) :
    sizeX_(sizeX),
    sizeY_(sizeY),
    sizeZ_(sizeZ),
    data_(sizeX * sizeY * sizeZ) {
  }

  //! Constructor to specify grid dimensions and initialization value

  //! \param[in] sizeX Grid size in x-dimension
  //! \param[in] sizeY Grid size in y-dimension
  //! \param[in] sizeZ Grid size in z-dimension
  //! \param[in] val Initial value for cells

  inline Grid(int sizeX, int sizeY, int sizeZ, T val) :
    sizeX_(sizeX),
    sizeY_(sizeY),
    sizeZ_(sizeZ),
    data_(sizeX * sizeY * sizeZ, val) {
  }

  //! Destructor

  //! Clears the vector of cell variables

  inline virtual ~Grid() {
    data_.clear();
  }

  // ============= //
  // Set operators //
  // ============= //

  //! Initialize the grid with given dimensions and initial value

  //! \param[in] sizeX Grid size in x-dimension
  //! \param[in] sizeY Grid size in y-dimension
  //! \param[in] sizeZ Grid size in z-dimension
  //! \param[in] val Initial value for cells

  void init( int sizeX, int sizeY, int sizeZ, T val ) {
    data_.assign( sizeX * sizeY * sizeZ, val );
  }

  //! Initialize the grid with given initial value

  //! \param[in] val Initial value for cells

  void init( T val ) {
    data_.assign( sizeX_ * sizeY_ * sizeZ_, val );
  }

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

  inline T& operator()(int x, int y, int z) {
    assert( x >= 0 && x < sizeX_ );
    assert( y >= 0 && y < sizeY_ );
    assert( z >= 0 && z < sizeZ_ );

    return data_[(z * sizeY_ + y) * sizeX_ + x];
  }

  //! Get specified constant cell variable of specified cell

  //! \param[in] x Cell coordinate in x-dimension
  //! \param[in] y Cell coordinate in y-dimension
  //! \param[in] z Cell coordinate in z-dimension
  //! \param[in] f Cell variable index
  //! \return Constant reference to specified cell variable of specified cell

  inline const T& operator()(int x, int y, int z) const {
    assert( x >= 0 && x < sizeX_ );
    assert( y >= 0 && y < sizeY_ );
    assert( z >= 0 && z < sizeZ_ );

    return data_[(z * sizeY_ + y) * sizeX_ + x];
  }

private:

  //! Number of cells in x-dimension

  int sizeX_;

  //! Number of cells in y-dimension

  int sizeY_;

  //! Number of cells in z-dimension

  int sizeZ_;

  //! Linearized, 1-dimensional representation of the 3D data grid

  std::vector<T> data_;

};

} // namespace lbm

#endif /* GRID_H_ */
