//! \file Emitter.h
//! \date   Mar 12, 2009
//! \author Florian Rathgeber

#ifndef EMITTER_H_
#define EMITTER_H_

#include <list>

#include "Particle.h"

//! Common namespace for all classes related to the particle system

namespace particles {

class ParticleSystem;

//! Adds new particles to the particle system

//! Particles are inserted at the emitter's position having a fixed initial
//! temperature and a lifetime depending on the emitter's current fuel value.

class Emitter {

  //! Fried declaration to allow ParticleSystem class access to protected
  //! members

  friend class ParticleSystem;

public:

  // ============================ //
  // Constructors and destructors //
  // ============================ //

  //! Constructor

  //! \param[in] pos             Position of the emitter (x,y,z)
  //! \param[in] temp            Initial temperature of emitted particles
  //! \param[in] fuel            Initial fuel value of the emitter
  //! \param[in] emitThreshold   Determines particle emission probability, the
  //!                            higher the threshold the smaller
  //! \param[in] fuelConsumption Factor by which fuel value is reduced every
  //!                            time a particle is emitted
  //! \param[in] lifetimeCoeff   Coefficient to determine particle lifetime
  //!                            based on emitters current fuel value

  Emitter ( core::vector3df pos,
            float temp,
            int fuel,
            int emitThreshold,
            float fuelConsumption,
            float lifetimeCoeff );

  //! Destructor

  //! Does nothing.

  virtual ~Emitter () {}

protected:

  // ============ //
  // Data members //
  // ============ //

  //! List of particles the emitter has already emitted

  std::list< Particle > particles_;

  //! Position of the emitter (x,y,z)

  core::vector3df pos_;

  //! Initial temperature of emitted particles

  float temp_;

  //! Current fuel value of the emitter

  int fuel_;

  //! Determines particle emission probability. The higher the threshold the
  //! smaller the emission probability for each time step

  int emitThreshold_;

  //! Factor by which fuel value is reduced every time a particle is emitted

  float fuelConsumption_;

  //! Coefficient to determine particle lifetime based on emitters current fuel
  //! value

  float lifetimeCoeff_;
};

} // namespace particles

#endif /* EMITTER_H_ */
