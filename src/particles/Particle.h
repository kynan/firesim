//! \file Particle.h
//! \date   Mar 12, 2009
//! \author Florian Rathgeber

#ifndef PARTICLE_H_
#define PARTICLE_H_

#define RAND1 .4 * std::rand() / (T) RAND_MAX - .2

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

  Particle ( scene::ISceneManager* mgr,
             s32 id,
             const core::vector3df &position,
             video::ITexture* texture,
             int numSprites,
             T temp,
             video::SColor& color,
             int lifetime ) :
               sprites_( numSprites ),
               pos_( position ),
               type_( FIRE ),
               temp_( temp ),
               lifetime_( lifetime ) {
    for ( int i = 0; i < numSprites; ++i ) {
      sprites_[i] = mgr->addBillboardSceneNode( NULL, core::dimension2df(10.0f,10.0f), position, id );
      sprites_[i]->setColor( color );
      sprites_[i]->setMaterialTexture( 0, texture );
      sprites_[i]->setMaterialFlag( video::EMF_LIGHTING, false );
  //    sprites_[i]->setMaterialType( video::EMT_SOLID );
//      sprites_[i]->setMaterialType( i % 2 ? video::EMT_TRANSPARENT_ADD_COLOR
//                                          : video::EMT_TRANSPARENT_ALPHA_CHANNEL );
      sprites_[i]->setMaterialType( video::EMT_TRANSPARENT_ALPHA_CHANNEL );
    }
  }

  virtual ~Particle () {
//    sprite_->remove();
  }

  T dist ( Particle& p ) {
    T dx = pos_.X - p.pos_.X;
    T dy = pos_.Y - p.pos_.Y;
    T dz = pos_.Y - p.pos_.Y;
    return sqrt( dx * dx + dy * dy + dz * dz );
  }

  void updatePos( const core::vector3df& d ) {
    pos_ += d;
    for ( int i = 0; i < sprites_.size(); ++i ) {
//      sprites_[i]->setPosition( pos_ + core::vector3df( RAND1, RAND1, RAND1 ) );
      sprites_[i]->setPosition( sprites_[i]->getPosition() + d + core::vector3df( RAND1, RAND1, RAND1 ) );
    }
  }

  const core::vector3df& getPos() {
    return pos_;
  }

  void clear() {
    for ( int i = 0; i < sprites_.size(); ++i ) {
      sprites_[i]->remove();
    }
    sprites_.clear();
  }

  void setColor( video::SColor c ) {
    for ( int i = 0; i < sprites_.size(); ++i ) {
      sprites_[i]->setColor( c );
    }
  }

  void setSize( T sz ) {
    for ( int i = 0; i < sprites_.size(); ++i ) {
      sprites_[i]->setSize( core::dimension2df( sz, sz ) );
    }
  }

protected:

  std::vector< scene::IBillboardSceneNode* > sprites_;

  core::vector3df pos_;

  ParticleType type_;

  T temp_;

  int lifetime_;

};

}

#endif /* PARTICLE_H_ */
