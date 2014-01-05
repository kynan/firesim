//! \file ParticleSystem.h
//! \date   Mar 12, 2009
//! \author Florian Rathgeber

#ifndef PARTICLESYSTEM_H_
#define PARTICLESYSTEM_H_

#ifndef RealType
#define RealType float
#endif

#include <irrlicht/irrlicht.h>

#include "../lbm/LBM.h"
#include "Emitter.h"

using namespace irr;

//! Common namespace for all classes related to the particle system

namespace particles {

struct Sphere {

  Sphere( float x,
          float y,
          float z,
          float radius,
          float u_x = 0.0,
          float u_y = 0.0,
          float u_z = 0.0,
          scene::IMeshSceneNode* n = 0 ) :
      pos( x, y, z ), r( radius ), u( u_x, u_y, u_z ), node( n ) {}

  void move() {
    pos += u;
    if ( node )
      node->setPosition( pos );
  }

  bool isPointInside( const core::vector3df& p ) {
    return   (p.X - pos.X) * (p.X - pos.X)
           + (p.Y - pos.Y) * (p.Y - pos.Y)
           + (p.Z - pos.Z) * (p.Z - pos.Z) < r * r;
  }

  core::vector3df pos;
  float r;
  core::vector3df u;
  scene::IMeshSceneNode* node;
};

//! Particle system that handles creation, movement and visualization of
//! particles.

//! Irrlicht 3D engine is used for OpenGL visualization.

class ParticleSystem {

public:

  // ============================ //
  // Constructors and destructors //
  // ============================ //

  //! Default constructor

  //! Does nothing.

  ParticleSystem() {}

  //! Constructor that initializes the particle system according to given
  //! configuration file

  //! \param[in] configFileName Path to configuration file to parse

  ParticleSystem ( std::string configFileName );

  //! Destructor

  //! Does nothing.

  virtual ~ParticleSystem () {}

  //! Set up the particle system according to configuration options

  //! \param[in] base Root block of the parsed configuration file

  void setup( ConfBlock& base );

  //! Main simulation loop

  void run();

protected:

  // ========================= //
  // Internal helper functions //
  // ========================= //

  //! Go over all particles and update position, color and size

  inline void updateParticles();

  //! Go over all emitters and emit new particles

  inline void emitParticles();

  //! Write a povray rendering file for a timestep

  //! \param[in] step Current timestep

  void writePovray( int step );

  //! Generate a black body color table

  //! \param[in] maxTemp Maximum emitted particle temperature

  void generateBlackBodyColorTable( float maxTemp );

  float getTime( timeval &start, timeval &end ) {
    return (float) ( end.tv_sec - start.tv_sec )
            + (float) ( end.tv_usec - start.tv_usec ) / 1000000.;
  }

  // ============ //
  // Data members //
  // ============ //

  //! Lattice Boltzmann fluid solver
  lbm::LBM<float> solver_;

  //! Vector containing all emitters of the particle system
  std::vector< Emitter > emitters_;

  //! Domain size in x-direction
  int sizeX_;

  //! Domain size in y-direction
  int sizeY_;

  //! Domain size in z-direction
  int sizeZ_;

  //! Total simulation timesteps
  int maxSteps_;

  //! Total number of generated particles
  int numParticles_;

  //! Thermal conservation coefficient
  float alpha_;

  //! Thermal transferability coefficient
  float beta_;

  //! Temperature threshold for a fire particle to turn into a smoke particle
  float smokeTemp_;

  //! Ambient temperature
  float ambTemp_;

  //! Thermal expansion coefficient
  float k_;

  //! Vector of inversed gravity
  core::vector3df gravity_;

  //! Basic particle size
  float sizeBase_;

  //! Variable particle size, multiplied with a factor depending on lifetime
  float sizeVar_;

  //! Number of sprites assigned to each particle
  int numSprites_;

  //! Table of a precomputed Gaussian distribution
  std::vector<float> gaussTable_;

  //! Black body color table
  std::vector< video::SColor > bbColorTable_;

  //! Irrlicht OpenGL device
  IrrlichtDevice* device_;

  //! Irrlicht scene manager
  scene::ISceneManager* smgr_;

  //! Irrlicht OpenGL video driver
  video::IVideoDriver* drvr_;

  //! Vector of particle textures
  std::vector< video::ITexture* > textures_;

  //! Vector of axis-aligned cuboid obstacles
  std::vector< core::aabbox3df > obstacles_;

  std::vector< Sphere > spheres_;

  //! Base file name for povray output files
  std::string povFileName_;

  std::string updFileName_;
  std::string irrFileName_;

  std::string screenshots_;
  int screenStep_;

  bool dynamicLights_;
};

} // namespace particles

#endif /* PARTICLESYSTEM_H_ */
