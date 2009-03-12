//! \file ParticleHandler.h
//! \date   Mar 12, 2009
//! \author Florian Rathgeber

#ifndef PARTICLEHANDLER_H_
#define PARTICLEHANDLER_H_

#include "../lbm/LBM.h"
#include "../confparser/ConfParser.h"
#include "Emitter.h"

namespace particles {

template<typename T>
class ParticleHandler {

public:

  ParticleHandler() {}

  ParticleHandler ( std::string configFileName ) {
    try {
      ConfParser p;
      ConfBlock base = p.parse( configFileName );
      std::cout << "Parsed configuration file " << configFileName << std::endl;
      setup( base );
      solver_.setup( base );
    } catch ( std::exception e ) {
      std::cerr << e.what() << std::endl;
      exit( -1 );
    }
  }

  virtual ~ParticleHandler () {}

  void setup( ConfBlock& base ) {
    try {

      // Read the parameters from the config file

      ConfBlock* paramBlock = base.find( "domain" );
      if ( paramBlock == NULL ) throw "No domain size given.";
      sizeX_ = paramBlock->getParam<int>( "sizeX" );
      sizeY_ = paramBlock->getParam<int>( "sizeY" );
      sizeZ_ = paramBlock->getParam<int>( "sizeZ" );
      std::cout << "Read domain specification:" << std::endl;
      std::cout << "sizeX : " << sizeX_ << std::endl;
      std::cout << "sizeY : " << sizeY_ << std::endl;
      std::cout << "sizeZ : " << sizeZ_ << std::endl;

      paramBlock = base.find( "parameters" );
      if ( paramBlock == NULL ) throw "No parameter specification given.";
      alpha_ = paramBlock->getParam<T>( "alpha" );
      beta_ = paramBlock->getParam<T>( "beta" );
      smokeTemp_ = paramBlock->getParam<T>( "smokeTemp" );
      maxSteps_ = paramBlock->getParam<int>( "maxSteps" );
      std::cout << "Read parameter specification:" << std::endl;
      std::cout << "alpha (conservation coefficient)  : " << alpha_ << std::endl;
      std::cout << "beta (transferability coefficient): " << beta_ << std::endl;
      std::cout << "Temperature threshold for smoke   : " << smokeTemp_ << std::endl;
      std::cout << "Number of steps                   : " << maxSteps_ << std::endl;

      paramBlock = base.find( "particles" );
      if ( paramBlock == NULL ) throw "No particle emitters specified.";
      ConfBlock::childIterPair cip = paramBlock->findAll( "emitter" );

      for ( ConfBlock::childIter it = cip.first; it != cip.second; ++it ) {

        ConfBlock b = it->second;

        T xStart = b.getParam<T>( "xStart" );
        T xEnd   = b.getParam<T>( "xEnd" );
        T yStart = b.getParam<T>( "yStart" );
        T yEnd   = b.getParam<T>( "yEnd" );
        T zStart = b.getParam<T>( "zStart" );
        T zEnd   = b.getParam<T>( "zEnd" );
        int xRange  = b.getParam<int>( "xRange" );
        int yRange  = b.getParam<int>( "yRange" );
        int zRange  = b.getParam<int>( "zRange" );
        T dx = ( xEnd - xStart ) / ( xRange + 1 );
        T dy = ( yEnd - yStart ) / ( yRange + 1 );
        T dz = ( zEnd - zStart ) / ( zRange + 1 );
        T temp   = b.getParam<T>( "temp" );
        int fuel    = b.getParam<int>( "fuel" );
        int emitThreshold = b.getParam<int>( "emitThreshold" );
        T fuelConsumption = b.getParam<T>( "fuelConsumption" );
        T lifetimeCoeff = b.getParam<T>( "lifetimeCoeff" );
        int i = 0, j = 0, k = 0;
        for ( T z = zStart; z <= zEnd && i <= zRange; z += dz, ++i ) {
          for ( T y = yStart; y <= yEnd && j <= yRange; y += dy, ++j ) {
            for ( T x = xStart; x <= xEnd && k <= xRange ; x += dx, ++k ) {
              emitters_.push_back( Emitter<T>( Vec3<T>( x, y, z ),
                                               temp,
                                               fuel,
                                               emitThreshold,
                                               fuelConsumption,
                                               lifetimeCoeff ) );
            }
          }
        }

      }
    } catch ( std::exception e ) {
      std::cerr << e.what() << std::endl;
      exit( -1 );
    } catch ( const char* e ) {
      std::cerr << e << std::endl;
      exit( -1 );
    }

    std::cout << "Setup finished!" << std::endl;
  }

  void run() {
    for ( int step = 0; step < maxSteps_; ++step ) {
      emitParticles();
      solver_.runStep();
      updateParticles();
    }
  }

protected:

  void updateParticles() {

    // Go over all particles of all emitters
    typename std::vector< Emitter<T> >::iterator ite, ite2;
    typename std::list< Particle<T> >::iterator itp, itp2;
    for ( ite = emitters_.begin(); ite != emitters_.end(); ++ite ) {
      for ( itp = (*ite).particles_.begin(); itp != (*ite).particles_.end(); ++itp) {
        Particle<T> p = *itp;

        // Move particle according to fluid velocity at current position
        p.pos_ += solver_.getVelocity( p.pos_ );
        // Remove particle if it has left the domain
        if ( p.pos_[0] < 1 || p.pos_[0] > sizeX_ - 1 ||
             p.pos_[1] < 1 || p.pos_[1] > sizeY_ - 1 ||
             p.pos_[2] < 1 || p.pos_[2] > sizeZ_ -1 )
          (*ite).particles_.erase( itp );

        // Update temperature
        T tempExt = 0.;
        // Add up temperature contributions of all other particles weighted by distance
        for ( ite2 = emitters_.begin(); ite2 != emitters_.end(); ++ite2 ) {
          for ( itp2 = (*ite2).particles_.begin(); itp2 != (*ite2).particles_.end(); ++itp2) {
            // TODO apply gaussian filter to distance value
            if ( &(*itp) != &p ) tempExt += p.dist( *itp ) * (*itp).temp_;
          }
        }
        p.temp_ = alpha_ * p.temp_ + beta_ * tempExt;
        // If temperature has fallen below threshold, convert to smoke particle
        if ( p.temp_ < smokeTemp_ ) p.type_ = SMOKE;
      }
    }
  }

  void emitParticles() {
    // Go over all emitters
    typename std::vector< Emitter<T> >::iterator ite;
    for ( ite = emitters_.begin(); ite != emitters_.end(); ++ite ) {
      (*ite).emit();
    }
  }

  lbm::LBM<T> solver_;

  std::vector< Emitter<T> > emitters_;

  T alpha_;
  T beta_;
  T smokeTemp_;

  int sizeX_;
  int sizeY_;
  int sizeZ_;
  int maxSteps_;

};

}

#endif /* PARTICLEHANDLER_H_ */
