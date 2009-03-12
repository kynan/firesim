//! \file Emitter.h
//! \date   Mar 12, 2009
//! \author Florian Rathgeber

#ifndef EMITTER_H_
#define EMITTER_H_

#include <list>
#include "Particle.h"

namespace particles {

template<typename T>
class ParticleHandler;

template<typename T>
class Emitter {

  friend class ParticleHandler<T>;

public:

  Emitter ( Vec3<T> pos,
            T temp,
            int fuel,
            int emitThreshold,
            T fuelConsumption,
            T lifetimeCoeff ) :
              pos_( pos ),
              temp_( temp ),
              fuel_( fuel ),
              emitThreshold_( emitThreshold ),
              fuelConsumption_( fuelConsumption ),
              lifetimeCoeff_( lifetimeCoeff ) {}

  virtual ~Emitter () {}

  void emit() {
    if ( fuel_ > std::rand() % emitThreshold_ ) {
      particles_.push_back( Particle<T>( temp_, pos_, fuel_ * lifetimeCoeff_ ) );
      fuel_ *= fuelConsumption_;
    }
  }

protected:

  std::list< Particle<T> > particles_;

  Vec3<T> pos_;

  T temp_;

  int fuel_;
  int emitThreshold_;
  T fuelConsumption_;
  T lifetimeCoeff_;
};

}

#endif /* EMITTER_H_ */
