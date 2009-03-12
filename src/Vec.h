//! \file Vec.h
//! \date   Jan 4, 2009
//! \author Florian Rathgeber

#ifndef VEC_H_
#define VEC_H_

#include <vector>

//! Three component vector of type specified by template parameter.

//! \tparam T Type of the vector elements

template< typename T >
class Vec3 {

public:

  // ============================ //
  // Constructors and destructors //
  // ============================ //

  //! Default constructor

  //! Creates a zero length vector

  Vec3() : data_( 0 ) {}

  //! Copy constructor

  //! \param[in] v Vector to copy from

  Vec3( const Vec3& v ) : data_( v.data_ ) {}

  //! Constructor to set all three components

  //! \param[in] x First component
  //! \param[in] y Second component
  //! \param[in] z Third component

  Vec3( T x, T y, T z ) : data_( 3 ) {
    data_[0] = x;
    data_[1] = y;
    data_[2] = z;
  }

  //! Destructor

  //! Clears the vector

  virtual ~Vec3() {
    data_.clear();
  }

  // ================ //
  // Access operators //
  // ================ //

  //! Non-constant access operator

  //! \param[in] i Vector component to return
  //! \return Reference to specified component of template type

  inline T& operator[]( int i ) { return data_[i]; }

  //! Constant access operator

  //! \param[in] i Vector component to return
  //! \return Const reference to specified component of template type

  inline const T& operator[]( int i ) const { return data_[i]; }

  inline Vec3<T> operator+( const Vec3<T>& v ) const {
    return Vec3<T>( data_[0] + v.data_[0], data_[1] + v.data_[1], data_[2] + v.data_[2] );
  }

  inline void operator+=( const Vec3<T>& v ) {
    data_[0] += v.data_[0];
    data_[1] += v.data_[1];
    data_[2] += v.data_[2];
  }

protected:

  // ============ //
  // Data members //
  // ============ //

  //! The stl vector to hold the data

  std::vector<T> data_;
};

#endif /* VEC_H_ */
