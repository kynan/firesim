//! \file Grid_def.h
//! Brief Implementation of the Grid class

//! \date   Jan 16, 2009
//! \author Florian Rathgeber

#ifndef GRID_DEF_H_
#define GRID_DEF_H_

// Implementation of general Grid class

// ============================ //
// Constructors and destructors //
// ============================ //

template<typename T, int Cellsize>
Grid<T,Cellsize>::Grid() :
  sizeX_(0), sizeY_(0), sizeZ_(0), data_(0) {
}

template<typename T, int Cellsize>
Grid<T,Cellsize>::Grid(int sizeX, int sizeY, int sizeZ) :
  sizeX_(sizeX),
  sizeY_(sizeY),
  sizeZ_(sizeZ),
  data_(Cellsize * sizeX * sizeY * sizeZ) {
}

template<typename T, int Cellsize>
Grid<T,Cellsize>::Grid(int sizeX, int sizeY, int sizeZ, T val) :
  sizeX_(sizeX),
  sizeY_(sizeY),
  sizeZ_(sizeZ),
  data_(Cellsize * sizeX * sizeY * sizeZ, val) {
}

template<typename T, int Cellsize>
Grid<T,Cellsize>::~Grid() {
  sizeX_ = 0;
  sizeY_ = 0;
  sizeZ_ = 0;
  data_.clear();
}

// ============= //
// Set operators //
// ============= //

template<typename T, int Cellsize>
void Grid<T,Cellsize>::init( int sizeX, int sizeY, int sizeZ, T val ) {
  sizeX_ = sizeX;
  sizeY_ = sizeY;
  sizeZ_ = sizeZ;
  data_.assign( Cellsize * sizeX * sizeY * sizeZ, val );
}

template<typename T, int Cellsize>
void Grid<T,Cellsize>::init( T val ) {
  data_.assign( Cellsize * sizeX_ * sizeY_ * sizeZ_, val );
}

// ================ //
// Access operators //
// ================ //

template<typename T, int Cellsize>
inline T& Grid<T,Cellsize>::operator()(int x, int y, int z, int f) {
  assert( x >= 0 && x < sizeX_ );
  assert( y >= 0 && y < sizeY_ );
  assert( z >= 0 && z < sizeZ_ );
  assert( f >= 0 && f < Cellsize );

  return data_[((z * sizeY_ + y) * sizeX_ + x) * Cellsize + f];
}

template<typename T, int Cellsize>
inline const T& Grid<T,Cellsize>::operator()(int x, int y, int z, int f) const {
  assert( x >= 0 && x < sizeX_ );
  assert( y >= 0 && y < sizeY_ );
  assert( z >= 0 && z < sizeZ_ );
  assert( f >= 0 && f < Cellsize );

  return data_[((z * sizeY_ + y) * sizeX_ + x) * Cellsize + f];
}

// Implementation of Grid class specialized to a single cell variable

// ============================ //
// Constructors and destructors //
// ============================ //

template<typename T>
Grid<T,1>::Grid() :
  sizeX_(0), sizeY_(0), sizeZ_(0), data_(0) {
}

template<typename T>
Grid<T,1>::Grid(int sizeX, int sizeY, int sizeZ) :
  sizeX_(sizeX),
  sizeY_(sizeY),
  sizeZ_(sizeZ),
  data_(sizeX * sizeY * sizeZ) {
}

template<typename T>
Grid<T,1>::Grid(int sizeX, int sizeY, int sizeZ, T val) :
  sizeX_(sizeX),
  sizeY_(sizeY),
  sizeZ_(sizeZ),
  data_(sizeX * sizeY * sizeZ, val) {
}

template<typename T>
Grid<T,1>::~Grid() {
  sizeX_ = 0;
  sizeY_ = 0;
  sizeZ_ = 0;
  data_.clear();
}

// ============= //
// Set operators //
// ============= //

template<typename T>
void Grid<T,1>::init( int sizeX, int sizeY, int sizeZ, T val ) {
  sizeX_ = sizeX;
  sizeY_ = sizeY;
  sizeZ_ = sizeZ;
  data_.assign( sizeX * sizeY * sizeZ, val );
}

template<typename T>
void Grid<T,1>::init( T val ) {
  data_.assign( sizeX_ * sizeY_ * sizeZ_, val );
}

// ================ //
// Access operators //
// ================ //

template<typename T>
inline T& Grid<T,1>::operator()(int x, int y, int z) {
  assert( x >= 0 && x < sizeX_ );
  assert( y >= 0 && y < sizeY_ );
  assert( z >= 0 && z < sizeZ_ );

  return data_[(z * sizeY_ + y) * sizeX_ + x];
}

template<typename T>
inline const T& Grid<T,1>::operator()(int x, int y, int z) const {
  assert( x >= 0 && x < sizeX_ );
  assert( y >= 0 && y < sizeY_ );
  assert( z >= 0 && z < sizeZ_ );

  return data_[(z * sizeY_ + y) * sizeX_ + x];
}

#endif /* GRID_DEF_H_ */
