//! \file Particle.h
//! \date   Mar 12, 2009
//! \author Florian Rathgeber

#ifndef PARTICLE_H_
#define PARTICLE_H_

#define RAND1 .4 * rand() / (float) RAND_MAX - .2

#include <vector>
#include <irrlicht/irrlicht.h>

using namespace irr;

//! Common namespace for all classes related to the particle system

namespace particles {

//! Possible different particle types

enum ParticleType { FIRE = 0, SMOKE = 1 };

class ParticleSystem;

//! Individual particle of the particle system

//! Has a configurable number of sprites (display primitives) assigned to it.

class Particle {

  //! Fried declaration to allow ParticleSystem class access to protected
  //! members

  friend class ParticleSystem;

public:

  // ============================ //
  // Constructors and destructors //
  // ============================ //

  //! Default constructor

  //! Does nothing.

  Particle () {}

  //! Constructor

  //! \param[in] mgr        Pointer to irrlicht scene manager
  //! \param[in] id         Particle ID
  //! \param[in] position   Initial particle position
  //! \param[in] texture    Texture to assign to sprites of this particle
  //! \param[in] numSprites Number of sprites assigned to each particle instance
  //! \param[in] temp       Initial particle temperature
  //! \param[in] color      Initial particle color
  //! \param[in] lifetime   Particle lifetime (in timesteps)

  Particle ( scene::ISceneManager* mgr,
             s32 id,
             const core::vector3df &position,
             video::ITexture* texture,
             int numSprites,
             float temp,
             video::SColor& color,
             int lifetime );

  //! Destructor

  //! Does nothing.

  virtual ~Particle () {}

  // ======= //
  // Getters //
  // ======= //

  //! Get distance between this and another particle

  //! \param[in] p Particle to measure distance from
  //! \return Distance between this and given particle

  float dist ( Particle& p ) {
    float dx = pos_.X - p.pos_.X;
    float dy = pos_.Y - p.pos_.Y;
    float dz = pos_.Y - p.pos_.Y;
    return sqrt( dx * dx + dy * dy + dz * dz );
  }

  //! Get current particle position

  //! \return Current particle position as vector

  const core::vector3df& getPos() {
    return pos_;
  }

  // ======= //
  // Setters //
  // ======= //

  //! Update particle position, including assigned sprites

  //! \param[in] d Displacement vector by which particle is to be moved

  void updatePos( const core::vector3df& d ) {
    pos_ += d;
    for ( uint i = 0; i < sprites_.size(); ++i ) {
      // Sprites are moved by given displacement + a random offset
      sprites_[i]->setPosition( sprites_[i]->getPosition() + d
                                + core::vector3df( RAND1, RAND1, RAND1 ) );
    }
  }

  //! Set the color for all assigned sprites

  //! \param[in] c Color to assign to all assigned sprites

  void setColor( video::SColor c ) {
    for ( uint i = 0; i < sprites_.size(); ++i ) {
      sprites_[i]->setColor( c );
    }
  }

  //! Set the size for all assigned sprites

  //! \param[in] sz Size to assign to all assigned sprites

  void setSize( float sz ) {
    for ( uint i = 0; i < sprites_.size(); ++i ) {
      sprites_[i]->setSize( core::dimension2df( sz, sz ) );
    }
  }

  //! Delete all sprites, i.e. remove them from the scene

  void clear() {
    for ( uint i = 0; i < sprites_.size(); ++i ) {
      sprites_[i]->remove();
    }
    sprites_.clear();
  }

protected:

  // ============ //
  // Data members //
  // ============ //

  //! Vector of sprites (display primitives) assigned to this particle
  std::vector< scene::IBillboardSceneNode* > sprites_;

  //! Current particle position
  core::vector3df pos_;

  //! Particle type (either FIRE or SMOKE)
  ParticleType type_;

  //! Current particle temperature
  float temp_;

  //! Current particle lifetime, decreased every timestep
  int lifetime_;

};

}

#endif /* PARTICLE_H_ */
