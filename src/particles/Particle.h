//! \file Particle.h
//! \date   Mar 12, 2009
//! \author Florian Rathgeber

#ifndef PARTICLE_H_
#define PARTICLE_H_

#include <irrlicht/irrlicht.h>

using namespace irr;

enum ParticleType { FIRE = 0, SMOKE = 1 };

namespace particles {

template<typename T>
class ParticleSystem;

template<typename T>
class Particle {

  friend class ParticleSystem<T>;

public:

  Particle () {}

  Particle ( scene::ISceneNode* parent,
             scene::ISceneManager* mgr,
             s32 id,
             const core::vector3df &position,
             video::ITexture* texture,
             T temp,
             int lifetime ) :
               sprite_( mgr->addBillboardSceneNode( parent, core::dimension2df(10.0f,10.0f), position, id ) ),
               type_( FIRE ),
               temp_( temp ),
               lifetime_( lifetime ) {
    assert( sprite_ );
    sprite_->setMaterialTexture( 0, texture );
    sprite_->setMaterialFlag( video::EMF_LIGHTING, false );
//    sprite_->setMaterialType( video::EMT_SOLID );
    sprite_->setMaterialType( video::EMT_TRANSPARENT_ADD_COLOR );
//    sprite_->setMaterialType( video::EMT_TRANSPARENT_ALPHA_CHANNEL );
  }

  virtual ~Particle () {
//    sprite_->remove();
  }

  T dist ( Particle& p ) {
    return getPos().getDistanceFrom( p.getPos() );
  }

  void updatePos( const core::vector3df& d ) {
    sprite_->setPosition( getPos() + d );
  }

  const core::vector3df& getPos() {
    return sprite_->getPosition();
  }

protected:

  scene::IBillboardSceneNode* sprite_;

  ParticleType type_;

  T temp_;

  int lifetime_;

};

}

#endif /* PARTICLE_H_ */
