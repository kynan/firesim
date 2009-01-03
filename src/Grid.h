/*
 * Grid.h
 *
 *  Created on: Jan 2, 2009
 *      Author: Florian Rathgeber
 */

#ifndef GRID_H_
#define GRID_H_

#include <assert.h>
#include <vector>

namespace lbm {

template<typename T, int Cellsize>
class Grid {

public:

  // Constructors and destructors

  inline Grid() :
    sizeX_(0), sizeY_(0), sizeZ_(0), data_(0) {
  }

  inline Grid(int xSize, int ySize, int zSize) :
    sizeX_(xSize),
    sizeY_(ySize),
    sizeZ_(zSize),
    data_(Cellsize * xSize * ySize * zSize) {
  }

  inline Grid(int xSize, int ySize, int zSize, T val) :
    sizeX_(xSize),
    sizeY_(ySize),
    sizeZ_(zSize),
    data_(Cellsize * xSize * ySize * zSize, val) {
  }

  inline virtual ~Grid() {
  }

  // Set operators

  void init( int xSize, int ySize, int zSize, T val ) {
    data_.assign( Cellsize * xSize * ySize * zSize, val );
  }

  void init( T val ) {
    data_.assign( Cellsize * sizeX_ * sizeY_ * sizeZ_, val );
  }

  // Access operators

  inline int getSizeX() {
    return sizeX_;
  }

  inline int getSizeY() {
    return sizeY_;
  }

  inline int getSizeZ() {
    return sizeZ_;
  }

  inline T& operator()(int x, int y, int z, int f) {
    assert( x >= 0 && x < sizeX_ );
    assert( y >= 0 && y < sizeY_ );
    assert( z >= 0 && z < sizeZ_ );
    assert( f >= 0 && f < Cellsize );

    return data_[((z * sizeY_ + y) * sizeX_ + x) * Cellsize + f];
  }

  inline const T& operator()(int x, int y, int z, int f) const {
    assert( x >= 0 && x < sizeX_ );
    assert( y >= 0 && y < sizeY_ );
    assert( z >= 0 && z < sizeZ_ );
    assert( f >= 0 && f < Cellsize );

    return data_[((z * sizeY_ + y) * sizeX_ + x) * Cellsize + f];
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


template<typename T>
class Grid<T,1> {

public:

  // Constructors and destructors

  inline Grid() :
    sizeX_(0), sizeY_(0), sizeZ_(0), data_(0) {
  }

  inline Grid(int xSize, int ySize, int zSize) :
    sizeX_(xSize),
    sizeY_(ySize),
    sizeZ_(zSize),
    data_(xSize * ySize * zSize) {
  }

  inline Grid(int xSize, int ySize, int zSize, T val) :
    sizeX_(xSize),
    sizeY_(ySize),
    sizeZ_(zSize),
    data_(xSize * ySize * zSize, val) {
  }

  inline virtual ~Grid() {
  }

  // Set operators

  void init( int xSize, int ySize, int zSize, T val ) {
    data_.assign( xSize * ySize * zSize, val );
  }

  void init( T val ) {
    data_.assign( sizeX_ * sizeY_ * sizeZ_, val );
  }

  // Access operators

  inline int getSizeX() {
    return sizeX_;
  }

  inline int getSizeY() {
    return sizeY_;
  }

  inline int getSizeZ() {
    return sizeZ_;
  }

  inline T& operator()(int x, int y, int z) {
    assert( x >= 0 && x < sizeX_ );
    assert( y >= 0 && y < sizeY_ );
    assert( z >= 0 && z < sizeZ_ );

    return data_[(z * sizeY_ + y) * sizeX_ + x];
  }

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
