//! \file Particle.cpp
//! \date   Mar 12, 2009
//! \author Florian Rathgeber

#include "Particle.h"

namespace particles {

Particle::Particle ( scene::ISceneManager* mgr,
                     s32 id,
                     const core::vector3df &position,
                     std::vector< video::ITexture* >& textures,
                     int numSprites,
                     float temp,
                     video::SColor& color,
                     float size,
                     int lifetime,
                     bool dynamicLights) :
             sprites_( numSprites ),
             pos_( position ),
             type_( FIRE ),
             temp_( temp ),
             lifetime_( lifetime ) {

  for ( int i = 0; i < numSprites; ++i ) {
    // Add sprite as billboard scene node to the scene
    sprites_[i] = mgr->addBillboardSceneNode( NULL,
                                              core::dimension2df(size,size),
                                              position,
                                              id );
    sprites_[i]->setColor( color );
    sprites_[i]->setMaterialTexture( 0, textures[ rand() % textures.size() ] );
    sprites_[i]->setMaterialFlag( video::EMF_LIGHTING, false );
    sprites_[i]->setMaterialFlag( video::EMF_ZWRITE_ENABLE, false );
    sprites_[i]->setMaterialType( video::EMT_TRANSPARENT_ADD_COLOR );

    // Add a light source to the billboard
    if ( dynamicLights )
      sprites_[i]->addChild(
          mgr->addLightSceneNode( sprites_[i], position, color, temp / 50 ) );
  }
}

} // namespace particles
