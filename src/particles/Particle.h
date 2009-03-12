//! \file Particle.h
//! \date   Mar 12, 2009
//! \author Florian Rathgeber

#ifndef PARTICLE_H_
#define PARTICLE_H_

#include "../Vec.h"

enum ParticleType { FIRE = 0, SMOKE = 1 };

namespace particles {

template<typename T>
class ParticleHandler;

template<typename T>
class Particle {

  friend class ParticleHandler<T>;

public:

  Particle () {}

  Particle ( T temp, Vec3<T> pos, int lifetime ) :
    type_( FIRE ), temp_( temp ), pos_( pos ), lifetime_( lifetime ) {}

  virtual ~Particle () {}

  T dist ( Particle& p ) {
    T dx = fabs( pos_[0] - p.pos_[0] );
    T dy = fabs( pos_[1] - p.pos_[1] );
    T dz = fabs( pos_[2] - p.pos_[2] );
    return dx * dx + dy * dy + dz * dz;
  }

protected:

  ParticleType type_;

  T temp_;

  Vec3<T> pos_;

  int lifetime_;

};

}

#endif /* PARTICLE_H_ */
