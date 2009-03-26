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

  // ====================== //
  // Mathematical operators //
  // ====================== //

  //! Addition operator

  //! \param[in] v Vector to add to this vector
  //! \return Vector sum of this vector and given vector

  inline Vec3<T> operator+( const Vec3<T>& v ) const {
    return Vec3<T>( data_[0] + v.data_[0], data_[1] + v.data_[1], data_[2] + v.data_[2] );
  }

  //! In-place addition operator

  //! \param[in] v Vector to add to this vector in-place

  inline void operator+=( const Vec3<T>& v ) {
    data_[0] += v.data_[0];
    data_[1] += v.data_[1];
    data_[2] += v.data_[2];
  }

  //! Inner product operator

  //! \param[in] v Vector to scalar multiply with this vector
  //! \return Scalar product of this and given vector

  inline T operator*( const Vec3<T>& v ) const {
    return data_[0] * v.data_[0] + data_[1] * v.data_[1] + data_[2] * v.data_[2];
  }

  //! Multiplication with scalar operator

  //! \param[in] a Scalar to multiply this vector with
  //! \return Vector scaled by given scalar

  inline Vec3<T> operator*( const T a ) const {
    return Vec3<T>(data_[0] * a + data_[1] * a + data_[2] * a);
  }

  //! In-place multiplication with scalar operator

  //! \param[in] a Scalar to multiply this vector with in-place

  inline void operator*=( const T a ) const {
    data_[0] *= a;
    data_[1] *= a;
    data_[2] *= a;
  }

protected:

  // ============ //
  // Data members //
  // ============ //

  //! The stl vector to hold the data

  std::vector<T> data_;
};

#endif /* VEC_H_ */
