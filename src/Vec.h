/*
 * Vec.h
 *
 *  Created on: Jan 4, 2009
 *      Author: Florian Rathgeber
 */

#ifndef VEC_H_
#define VEC_H_

#include <vector>

template< typename T >
class Vec3 {

public:

  // Constructors and destructors

  Vec3() : data_( 0 ) {}

  Vec3( T x, T y, T z ) : data_( 3 ) {
    data_[0] = x;
    data_[1] = y;
    data_[2] = z;
  }

  virtual ~Vec3() {
    data_.clear();
  }

  // access operators

  inline T& operator[]( int i ) { return data_[i]; }

  inline const T& operator[]( int i ) const { return data_[i]; }

private:
  std::vector<T> data_;
};

#endif /* VEC_H_ */
