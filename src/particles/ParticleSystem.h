//! \file ParticleSystem.h
//! \date   Mar 12, 2009
//! \author Florian Rathgeber

#ifndef PARTICLEHANDLER_H_
#define PARTICLEHANDLER_H_

#include <irrlicht/irrlicht.h>
#include "../lbm/LBM.h"
#include "../confparser/ConfParser.h"
#include "Emitter.h"

using namespace irr;

namespace particles {

template<typename T>
class ParticleSystem {

public:

  ParticleSystem() {}

  ParticleSystem ( std::string configFileName ) : numParticles_( 0 ) {
    try {

      ConfParser p;
      ConfBlock base = p.parse( configFileName );
      std::cout << "Parsed configuration file " << configFileName << std::endl;

      device_ = createDevice( video::EDT_OPENGL, core::dimension2di(800,600), 24 );
      if (!device_) throw "Could not create OpenGL device.";
      smgr_ = device_->getSceneManager();
      drvr_ = device_->getVideoDriver();
      assert( smgr_ && drvr_ );

      setup( base );
      solver_.setup( base );

    } catch ( std::exception e ) {
      std::cerr << e.what() << std::endl;
      exit( -1 );
    } catch ( const char* e ) {
      std::cerr << e << std::endl;
      exit( -1 );
    }
  }

  virtual ~ParticleSystem () {}

  void setup( ConfBlock& base ) {
    try {

      // Read the parameters from the config file

      std::cout << "Setting up ParticleHandler..." << std::endl;

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
      k_ = paramBlock->getParam<T>( "k" );
      T g_x = paramBlock->getParam<T>( "g_x" );
      T g_y = paramBlock->getParam<T>( "g_y" );
      T g_z = paramBlock->getParam<T>( "g_z" );
      gravity_ = core::vector3df( g_x, g_y, g_z );
      smokeTemp_ = paramBlock->getParam<T>( "smokeTemp" );
      ambTemp_ = paramBlock->getParam<T>( "ambTemp" );
      maxSteps_ = paramBlock->getParam<int>( "maxSteps" );
      std::cout << "Read parameter specification:" << std::endl;
      std::cout << "alpha (conservation coefficient)   : " << alpha_ << std::endl;
      std::cout << "beta (transferability coefficient) : " << beta_ << std::endl;
      std::cout << "k (thermal expansion coefficient)  : " << k_ << std::endl;
      std::cout << "Gravity unit vector                : <" << g_x << "," << g_y << "," << g_z << ">" << std::endl;
      std::cout << "Temperature threshold for smoke (K): " << smokeTemp_ << std::endl;
      std::cout << "Ambient temperature (K)            : " << ambTemp_ << std::endl;
      std::cout << "Number of steps                    : " << maxSteps_ << std::endl;

      // Precompute Gauss function for thermal diffusion
      int maxElem = sqrt( sizeX_ * sizeX_ + sizeY_ * sizeY_ + sizeZ_ * sizeZ_ );
      gaussTable_.reserve( maxElem );
      T a = 1. / sqrt( 2. * M_PI );
      for ( int i = 0; i < maxElem; ++i ) {
        gaussTable_.push_back( a * exp(-  i * i / 2. ) );
      }

      paramBlock = base.find( "povray" );
      if ( paramBlock != NULL ) {
        povFileName_ = paramBlock->getParam<std::string>( "povFileName" );
        std::cout << "Povray output specification:" << std::endl;
        std::cout << "Povray output file base name : " << povFileName_ << std::endl;
      } else {
        std::cout << "No povray block given in configuration file, no output will be created." << std::endl;
        povFileName_ = "";
      }

      paramBlock = base.find( "particles" );
      if ( paramBlock == NULL ) throw "No particle emitters specified.";
      ConfBlock::childIterPair cip = paramBlock->findAll( "emitter" );

      T maxTemp = 0;
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
        T dx = (xRange <= 1) ? 0. : ( xEnd - xStart ) / ( xRange - 1 );
        T dy = (yRange <= 1) ? 0. : ( yEnd - yStart ) / ( yRange - 1 );
        T dz = (zRange <= 1) ? 0. : ( zEnd - zStart ) / ( zRange - 1 );
        T temp   = b.getParam<T>( "temp" );
        if ( temp > maxTemp ) maxTemp = temp;
        int fuel    = b.getParam<int>( "fuel" );
        int emitThreshold = b.getParam<int>( "emitThreshold" );
        T fuelConsumption = b.getParam<T>( "fuelConsumption" );
        T lifetimeCoeff = b.getParam<T>( "lifetimeCoeff" );
        std::cout << "Read emitter specification:" << std::endl;
        std::cout << "xStart : " << xStart << std::endl;
        std::cout << "xEnd   : " << xEnd << std::endl;
        std::cout << "yStart : " << yStart << std::endl;
        std::cout << "yEnd   : " << yEnd << std::endl;
        std::cout << "zStart : " << zStart << std::endl;
        std::cout << "zEnd   : " << zEnd << std::endl;
        std::cout << "xRange : " << xRange << std::endl;
        std::cout << "yRange : " << yRange << std::endl;
        std::cout << "zRange : " << zRange << std::endl;
        std::cout << "dx     : " << dx << std::endl;
        std::cout << "dy     : " << dy << std::endl;
        std::cout << "dz     : " << dz << std::endl;
        std::cout << "temp : " << temp << std::endl;
        std::cout << "fuel : " << fuel << std::endl;
        std::cout << "emitThreshold : " << emitThreshold << std::endl;
        std::cout << "fuelConsumption : " << fuelConsumption << std::endl;
        std::cout << "lifetimeCoeff : " << lifetimeCoeff << std::endl;

        int i, j, k;
        T x, y, z;
        for ( z = zStart, i = 0; z <= zEnd && i < zRange; z += dz, ++i ) {
          for ( y = yStart, j = 0; y <= yEnd && j < yRange; y += dy, ++j ) {
            for ( x = xStart, k = 0; x <= xEnd && k < xRange; x += dx, ++k ) {
              std::cout << "Create emiter at <" << x << "," << y << "," << z << ">" << std::endl;
              emitters_.push_back( Emitter<T>( core::vector3df( x, y, z ),
                                               temp,
                                               fuel,
                                               emitThreshold,
                                               fuelConsumption,
                                               lifetimeCoeff ) );
            }
          }
        }

      }

      paramBlock = base.find( "irrlicht" );
      if ( paramBlock == NULL ) throw "No irrlicht configuration specified.";
      cip = paramBlock->findAll( "texture" );

      for ( ConfBlock::childIter it = cip.first; it != cip.second; ++it ) {

        ConfBlock b = it->second;

        std::string file = b.getParam<std::string>( "file" );
        textures_.push_back( drvr_->getTexture( file.c_str() ) );
      }

      generateBlackBodyColorTable( 2. * maxTemp );

    } catch ( std::exception e ) {
      std::cerr << e.what() << std::endl;
      exit( -1 );
    } catch ( const char* e ) {
      std::cerr << e << std::endl;
      exit( -1 );
    }

    std::cout << "ParticleHandler setup finished!" << std::endl;
  }

  void run() {

    device_->setWindowCaption(L"LBM Reference");
    smgr_->addCameraSceneNode( smgr_->getRootSceneNode(),
                               core::vector3df( sizeX_/2, sizeY_/2, -sizeZ_ ),
                               core::vector3df( sizeX_/2, sizeY_/2, 0 ) );
//    smgr_->addCameraSceneNodeFPS();
//    device_->getCursorControl()->setVisible( false );

    int step = 0;
    while ( device_->run() ) {

      // emit particles
      emitParticles();
//      if ( povFileName_.length() > 0 ) writePovray( step );

      drvr_->beginScene(true, true, video::SColor(255,100,100,100));
      smgr_->drawAll();
      drvr_->endScene();

      solver_.runStep();
      updateParticles();

      core::stringw str = L"LBM fire simulation [";
      str += drvr_->getName();
      str += "L], FPS: ";
      str += (s32)drvr_->getFPS();
      str += L", Particles: ";
      str += numParticles_;
      str += L", Time step: ";
      str += step;
      str += L" / ";
      str += maxSteps_;
      device_->setWindowCaption( str.c_str() );

      if ( ++step >= maxSteps_ ) break;
    }
  }

protected:

  void updateParticles() {

    // Go over all particles of all emitters
    typename std::vector< Emitter<T> >::iterator ite, ite2;
    typename std::list< Particle<T> >::iterator itp, itp2;
    for ( ite = emitters_.begin(); ite != emitters_.end(); ++ite ) {
      for ( itp = (*ite).particles_.begin(); itp != (*ite).particles_.end(); ++itp) {
        assert( (*itp).sprite_ );

        // Remove particles that have left the domain or exceeded their lifetime
        while ( itp != (*ite).particles_.end() && (
                (*itp).getPos().X < 1 || (*itp).getPos().X > sizeX_ - 1 ||
                (*itp).getPos().Y < 1 || (*itp).getPos().Y > sizeY_ - 1 ||
                (*itp).getPos().Z < 1 || (*itp).getPos().Z > sizeZ_ -1  ||
                (*itp).lifetime_ < 1 || (*itp).temp_ < ambTemp_ ) ) {
          // DEBUG output
          std::cout << "Removing particle " << (*itp).sprite_->getID();
          std::cout << " at position <" << (*itp).getPos().X << "," << (*itp).getPos().Y << "," << (*itp).getPos().Z << ">";
          std::cout << " with lifetime " << (*itp).lifetime_ << ", Remaining: " << (*ite).particles_.size() << std::endl;
          (*itp).sprite_->remove();
          itp = (*ite).particles_.erase( itp );
        }

        if ( itp == (*ite).particles_.end() ) break;

        // Move particle according to fluid velocity at current position plus
        // the buoyancy force
        Vec3<T> tmp = solver_.getVelocity( (*itp).getPos().X, (*itp).getPos().Y, (*itp).getPos().Z );
        (*itp).updatePos( core::vector3df( tmp[0], tmp[1], tmp[2] )
                          + gravity_ * k_ * ( (*itp).temp_ - ambTemp_ ) );
        // DEBUG output
//        std::cout << "Velocity after update: <" << (*itp).pos_[0] << "," << (*itp).pos_[1] << "," << (*itp).pos_[2] << ">" << std::endl;

        if ( (*itp).type_ == FIRE ) {
          // Update temperature
          T tempExt = - gaussTable_[0] * (*itp).temp_;
          // Add up temperature contributions of all other particles weighted by distance
          for ( ite2 = emitters_.begin(); ite2 != emitters_.end(); ++ite2 ) {
            for ( itp2 = (*ite2).particles_.begin(); itp2 != (*ite2).particles_.end(); ++itp2) {
              tempExt += gaussTable_[ (int) (*itp).dist( *itp2 ) ] * (*itp2).temp_;
            }
          }
          (*itp).temp_ = alpha_ * (*itp).temp_ + beta_ * tempExt;
          // DEBUG output
//          std::cout << "Temperature of particle " << (*itp).sprite_->getID() << ": ";
//          std::cout << (*itp).temp_ << ", table entry " << (int) (((*itp).temp_ - smokeTemp_) / 50.) << std::endl;
          (*itp).sprite_->setColor( bbColorTable_[ (int) (((*itp).temp_ - smokeTemp_) / 50.) ] );
          // If temperature has fallen below threshold, convert to smoke particle
          if ( (*itp).temp_ < smokeTemp_ ) {
            (*itp).type_ = SMOKE;
            (*itp).sprite_->setColor( video::SColor( 50,180,180,180 ) );
          }
        }

        // Update lifetime and particle size
        if ( (*itp).lifetime_-- % 10 == 0 ) {
          T sz = 16.0 * (*itp).lifetime_ / ( (*ite).fuel_ * (*ite).lifetimeCoeff_ );
          (*itp).sprite_->setSize( core::dimension2df( sz, sz ) );
        }
      }
    }
  }

  void emitParticles() {
    // Go over all emitters
    typename std::vector< Emitter<T> >::iterator ite;
    for ( ite = emitters_.begin(); ite != emitters_.end(); ++ite ) {
//      (*ite).emit();
      if ( (*ite).fuel_ > std::rand() % (*ite).emitThreshold_ ) {
        (*ite).particles_.push_back( Particle<T>( smgr_->getRootSceneNode(),
                                     smgr_,
                                     numParticles_++,
                                     (*ite).pos_,
                                     textures_[ std::rand() % textures_.size() ],
                                     (*ite).temp_,
                                     (int) ((*ite).fuel_ * (*ite).lifetimeCoeff_ )) );
        (*ite).fuel_ *= (*ite).fuelConsumption_;
      }
    }
  }

  void writePovray( int step ) {

    // Open file for writing
    std::ostringstream oss;
    oss << povFileName_ << "." << step << ".pov";
    std::ofstream povFile( oss.str().c_str(), std::ios::out );

    povFile << "// -w640 -h480 +a0.3\n";
    povFile << "global_settings { ambient_light rgb<1,1,1> }\n";
    povFile << "camera {\n";
    povFile << "  location <" << sizeX_ / 2 << "," << sizeZ_ / 2 << ",-20>\n";
    povFile << "  look_at <" << sizeX_ / 2 << "," << sizeZ_ / 2 << ",0>\n";
    povFile << "}\n";
    povFile << "light_source { <" << sizeX_ / 2 << "," << sizeY_ / 2 << ",-20>, rgb<1,1,1> }\n";
    typename std::vector< Emitter<T> >::iterator ite;
    typename std::list< Particle<T> >::iterator itp;
    for ( ite = emitters_.begin(); ite != emitters_.end(); ++ite ) {
      for ( itp = (*ite).particles_.begin(); itp != (*ite).particles_.end(); ++itp) {
        Particle<T> p = *itp;
        povFile << "disc {\n";
        povFile << "  <" << p.pos_[0] << "," << p.pos_[2] << "," << p.pos_[1] << ">, ";
        povFile << "-z, " << p.lifetime_ / 100. << "\n";
        if ( p.type_ == FIRE ) povFile << "  texture { pigment { color rgb<1,0,0> } }\n";
        else povFile << "  texture { pigment { color rgb<.5,.5,.5> } }\n";
        povFile << "}\n";
      }
    }
  }

  void generateBlackBodyColorTable( T maxTemp ) {

    T cie_colour_match[81][3] = {
        {0.0014,0.0000,0.0065}, {0.0022,0.0001,0.0105}, {0.0042,0.0001,0.0201},
        {0.0076,0.0002,0.0362}, {0.0143,0.0004,0.0679}, {0.0232,0.0006,0.1102},
        {0.0435,0.0012,0.2074}, {0.0776,0.0022,0.3713}, {0.1344,0.0040,0.6456},
        {0.2148,0.0073,1.0391}, {0.2839,0.0116,1.3856}, {0.3285,0.0168,1.6230},
        {0.3483,0.0230,1.7471}, {0.3481,0.0298,1.7826}, {0.3362,0.0380,1.7721},
        {0.3187,0.0480,1.7441}, {0.2908,0.0600,1.6692}, {0.2511,0.0739,1.5281},
        {0.1954,0.0910,1.2876}, {0.1421,0.1126,1.0419}, {0.0956,0.1390,0.8130},
        {0.0580,0.1693,0.6162}, {0.0320,0.2080,0.4652}, {0.0147,0.2586,0.3533},
        {0.0049,0.3230,0.2720}, {0.0024,0.4073,0.2123}, {0.0093,0.5030,0.1582},
        {0.0291,0.6082,0.1117}, {0.0633,0.7100,0.0782}, {0.1096,0.7932,0.0573},
        {0.1655,0.8620,0.0422}, {0.2257,0.9149,0.0298}, {0.2904,0.9540,0.0203},
        {0.3597,0.9803,0.0134}, {0.4334,0.9950,0.0087}, {0.5121,1.0000,0.0057},
        {0.5945,0.9950,0.0039}, {0.6784,0.9786,0.0027}, {0.7621,0.9520,0.0021},
        {0.8425,0.9154,0.0018}, {0.9163,0.8700,0.0017}, {0.9786,0.8163,0.0014},
        {1.0263,0.7570,0.0011}, {1.0567,0.6949,0.0010}, {1.0622,0.6310,0.0008},
        {1.0456,0.5668,0.0006}, {1.0026,0.5030,0.0003}, {0.9384,0.4412,0.0002},
        {0.8544,0.3810,0.0002}, {0.7514,0.3210,0.0001}, {0.6424,0.2650,0.0000},
        {0.5419,0.2170,0.0000}, {0.4479,0.1750,0.0000}, {0.3608,0.1382,0.0000},
        {0.2835,0.1070,0.0000}, {0.2187,0.0816,0.0000}, {0.1649,0.0610,0.0000},
        {0.1212,0.0446,0.0000}, {0.0874,0.0320,0.0000}, {0.0636,0.0232,0.0000},
        {0.0468,0.0170,0.0000}, {0.0329,0.0119,0.0000}, {0.0227,0.0082,0.0000},
        {0.0158,0.0057,0.0000}, {0.0114,0.0041,0.0000}, {0.0081,0.0029,0.0000},
        {0.0058,0.0021,0.0000}, {0.0041,0.0015,0.0000}, {0.0029,0.0010,0.0000},
        {0.0020,0.0007,0.0000}, {0.0014,0.0005,0.0000}, {0.0010,0.0004,0.0000},
        {0.0007,0.0002,0.0000}, {0.0005,0.0002,0.0000}, {0.0003,0.0001,0.0000},
        {0.0002,0.0001,0.0000}, {0.0002,0.0001,0.0000}, {0.0001,0.0000,0.0000},
        {0.0001,0.0000,0.0000}, {0.0001,0.0000,0.0000}, {0.0000,0.0000,0.0000}
    };

    for ( T t = smokeTemp_; t < maxTemp; t += 50. ) {

      // Calculate x,y,z colors from solar spectrum

      int i;
      T lambda, x = 0, y = 0, z = 0, xyz;
      for (i = 0, lambda = 380; lambda < 780.1; i++, lambda += 5) {
          T Me;
          // Get black body radiation intensity for given temperature and
          // wavelength
          T wlm = lambda * 1e-9; // wavelength in meters
          Me = (3.74183e-16 * pow(wlm, -5.0)) / (exp(1.4388e-2 / (wlm * t)) - 1.0);
          x += Me * cie_colour_match[i][0];
          y += Me * cie_colour_match[i][1];
          z += Me * cie_colour_match[i][2];
      }
      xyz = (x + y + z);
      x /= xyz;
      y /= xyz;
      z /= xyz;

      // Calculate r,g,b colors from x,y,z colors

      T xr, yr, zr, xg, yg, zg, xb, yb, zb;
      T xw, yw, zw;
      T rx, ry, rz, gx, gy, gz, bx, by, bz;
      T rw, gw, bw;

      xr = 0.630;  yr = 0.340;  zr = 1 - (xr + yr);
      xg = 0.310;  yg = 0.595;  zg = 1 - (xg + yg);
      xb = 0.155;  yb = 0.070;  zb = 1 - (xb + yb);

      xw = 0.3127; yw = 0.3291; zw = 1 - (xw + yw);

      // xyz -> rgb matrix, before scaling to white.
      rx = (yg * zb) - (yb * zg);  ry = (xb * zg) - (xg * zb);  rz = (xg * yb) - (xb * yg);
      gx = (yb * zr) - (yr * zb);  gy = (xr * zb) - (xb * zr);  gz = (xb * yr) - (xr * yb);
      bx = (yr * zg) - (yg * zr);  by = (xg * zr) - (xr * zg);  bz = (xr * yg) - (xg * yr);

      // White scaling factors.
      // Dividing by yw scales the white luminance to unity, as conventional.
      rw = ((rx * xw) + (ry * yw) + (rz * zw)) / yw;
      gw = ((gx * xw) + (gy * yw) + (gz * zw)) / yw;
      bw = ((bx * xw) + (by * yw) + (bz * zw)) / yw;

      // xyz -> rgb matrix, correctly scaled to white.
      rx = rx / rw;  ry = ry / rw;  rz = rz / rw;
      gx = gx / gw;  gy = gy / gw;  gz = gz / gw;
      bx = bx / bw;  by = by / bw;  bz = bz / bw;

      // rgb of the desired point
      T r = (rx * x) + (ry * y) + (rz * z);
      T g = (gx * x) + (gy * y) + (gz * z);
      T b = (bx * x) + (by * y) + (bz * z);

      // Constrain to RGB gammut
      double w;
      // Amount of white needed is w = - min(0, *r, *g, *b)
      w = (0 < r) ? 0 : r;
      w = (w < g) ? w : g;
      w = (w < b) ? w : b;
      w = -w;
      // Add just enough white to make r, g, b all positive.
      if (w > 0) {
          r += w;  g += w; b += w;
      }

      // Normalize
      T max = (r > b) ? ( (r > g) ? r : g ) : ( (b > g) ? b : g );
      r *= 255. / max;
      g *= 255. / max;
      b *= 255. / max;

      bbColorTable_.push_back( video::SColor( 255, r, g, b ) );
      std::cout << "Temperature " << t << "K: RGB <" << r << ", " << g << ", " << b << ">" << std::endl;
    }
    std::cout << "Generated black body color table from " << smokeTemp_;
    std::cout << "K to " << maxTemp << "K (" << bbColorTable_.size() << " values)" << std::endl;
  }

  lbm::LBM<T> solver_;

  std::vector< Emitter<T> > emitters_;

  T alpha_;
  T beta_;
  T smokeTemp_;
  T ambTemp_;
  T k_;
  core::vector3df gravity_;
  std::vector<T> gaussTable_;
  std::vector< video::SColor > bbColorTable_;

  IrrlichtDevice* device_;
  scene::ISceneManager* smgr_;
  video::IVideoDriver* drvr_;
  std::vector< video::ITexture* > textures_;

  int sizeX_;
  int sizeY_;
  int sizeZ_;
  int maxSteps_;
  int numParticles_;

  std::string povFileName_;
};

}

#endif /* PARTICLEHANDLER_H_ */
