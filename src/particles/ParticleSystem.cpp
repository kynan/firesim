//! \file ParticleSystem.cpp
//! \date   Mar 12, 2009
//! \author Florian Rathgeber

#include "../confparser/ConfParser.h"
#include "ParticleSystem.h"
#include "../lbm/LBM_def.h"

namespace particles {

ParticleSystem::ParticleSystem ( std::string configFileName )
    : numParticles_( 0 ) {
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

    srand ( time(NULL) );

  } catch ( std::exception e ) {
    std::cerr << e.what() << std::endl;
    exit( -1 );
  } catch ( const char* e ) {
    std::cerr << e << std::endl;
    exit( -1 );
  }
}

void ParticleSystem::setup( ConfBlock& base ) {
  try {

    // Read parameters from config file

    std::cout << "Setting up ParticleSystem..." << std::endl;

    // Read domain specification
    ConfBlock* paramBlock = base.find( "domain" );
    if ( paramBlock == NULL ) throw "No domain size given.";
    sizeX_ = paramBlock->getParam<int>( "sizeX" );
    sizeY_ = paramBlock->getParam<int>( "sizeY" );
    sizeZ_ = paramBlock->getParam<int>( "sizeZ" );
    std::cout << "Read domain specification:" << std::endl;
    std::cout << "sizeX : " << sizeX_ << std::endl;
    std::cout << "sizeY : " << sizeY_ << std::endl;
    std::cout << "sizeZ : " << sizeZ_ << std::endl;

    // Read parameter specification
    paramBlock = base.find( "parameters" );
    if ( paramBlock == NULL ) throw "No parameter specification given.";
    alpha_ = paramBlock->getParam<float>( "alpha" );
    beta_ = paramBlock->getParam<float>( "beta" );
    k_ = paramBlock->getParam<float>( "k" );
    float g_x = paramBlock->getParam<float>( "g_x" );
    float g_y = paramBlock->getParam<float>( "g_y" );
    float g_z = paramBlock->getParam<float>( "g_z" );
    gravity_ = core::vector3df( g_x, g_y, g_z );
    smokeTemp_ = paramBlock->getParam<float>( "smokeTemp" );
    ambTemp_ = paramBlock->getParam<float>( "ambTemp" );
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
    float a = 1. / sqrt( 2. * M_PI );
    for ( int i = 0; i < maxElem; ++i ) {
      gaussTable_.push_back( a * exp(-  i * i / 2. ) );
    }

    // Set up povray output
    paramBlock = base.find( "povray" );
    if ( paramBlock != NULL ) {
      povFileName_ = paramBlock->getParam<std::string>( "povFileName" );
      std::cout << "Povray output specification:" << std::endl;
      std::cout << "Povray output file base name : " << povFileName_ << std::endl;
    } else {
      std::cout << "No povray block given in configuration file, no output will be created." << std::endl;
      povFileName_ = "";
    }

    // Set up bounding boxes of obstacles
    paramBlock = base.find( "obstacles" );
    if ( paramBlock == NULL ) {
      std::cout << "No obstacles defined." << std::endl;
    } else {

      std::cout << "Set up the obstacles..." << std::endl;

      ConfBlock::childIterPair bit = paramBlock->findAll( "cuboid_stationary" );
      for ( ConfBlock::childIter it = bit.first; it != bit.second; ++it ) {

        ConfBlock& b = it->second;

        int xStart = b.getParam<int>( "xStart" );
        int xEnd   = b.getParam<int>( "xEnd" ) + 1;
        int yStart = b.getParam<int>( "yStart" );
        int yEnd   = b.getParam<int>( "yEnd" ) + 1;
        int zStart = b.getParam<int>( "zStart" );
        int zEnd   = b.getParam<int>( "zEnd" ) + 1;
        std::cout << "Stationary cuboid ranging from <" << xStart << ",";
        std::cout << yStart << "," << zStart << "> to <" << xEnd << "," << yEnd;
        std::cout << "," << zEnd << ">" << std::endl;
        obstacles_.push_back(
            core::aabbox3df( xStart, yStart, zStart, xEnd, yEnd, zEnd ) );
      }
    }

    // Set up particle system
    paramBlock = base.find( "particles" );
    if ( paramBlock == NULL ) throw "No particle emitters specified.";
    std::cout << "Set up the particle system..." << std::endl;
    ConfBlock::childIterPair cip = paramBlock->findAll( "emitter" );
    sizeBase_ = paramBlock->getParam<float>( "sizeBase" );
    sizeVar_ = paramBlock->getParam<float>( "sizeVar" );
    numSprites_ = paramBlock->getParam<int>( "numSprites" );

    float maxTemp = 0;
    for ( ConfBlock::childIter it = cip.first; it != cip.second; ++it ) {

      ConfBlock b = it->second;

      float xStart = b.getParam<float>( "xStart" );
      float xEnd   = b.getParam<float>( "xEnd" );
      float yStart = b.getParam<float>( "yStart" );
      float yEnd   = b.getParam<float>( "yEnd" );
      float zStart = b.getParam<float>( "zStart" );
      float zEnd   = b.getParam<float>( "zEnd" );
      int xRange  = b.getParam<int>( "xRange" );
      int yRange  = b.getParam<int>( "yRange" );
      int zRange  = b.getParam<int>( "zRange" );
      float dx = (xRange <= 1) ? 0. : ( xEnd - xStart ) / ( xRange - 1 );
      float dy = (yRange <= 1) ? 0. : ( yEnd - yStart ) / ( yRange - 1 );
      float dz = (zRange <= 1) ? 0. : ( zEnd - zStart ) / ( zRange - 1 );
      float temp   = b.getParam<float>( "temp" );
      if ( temp > maxTemp ) maxTemp = temp;
      int fuel    = b.getParam<int>( "fuel" );
      int emitThreshold = b.getParam<int>( "emitThreshold" );
      float fuelConsumption = b.getParam<float>( "fuelConsumption" );
      float lifetimeCoeff = b.getParam<float>( "lifetimeCoeff" );
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
      float x, y, z;
      for ( z = zStart, i = 0; z <= zEnd && i < zRange; z += dz, ++i ) {
        for ( y = yStart, j = 0; y <= yEnd && j < yRange; y += dy, ++j ) {
          for ( x = xStart, k = 0; x <= xEnd && k < xRange; x += dx, ++k ) {
            std::cout << "Create emiter at <" << x << "," << y << "," << z << ">" << std::endl;
            emitters_.push_back( Emitter( core::vector3df( x, y, z ),
                                             temp,
                                             fuel,
                                             emitThreshold,
                                             fuelConsumption,
                                             lifetimeCoeff ) );
          }
        }
      }

    }

    // Set up irrlicht engine
    paramBlock = base.find( "irrlicht" );
    if ( paramBlock == NULL ) throw "No irrlicht configuration specified.";
    std::cout << "Set up irrlicht configuration..." << std::endl;
    cip = paramBlock->findAll( "texture" );

    for ( ConfBlock::childIter it = cip.first; it != cip.second; ++it ) {

      ConfBlock b = it->second;

      std::string file = b.getParam<std::string>( "file" );
      textures_.push_back( drvr_->getTexture( file.c_str() ) );
    }

    // Generate black body color table
    std::cout << "Generate black body color table..." << std::endl;
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

void ParticleSystem::run() {

  // Add camera to scene
  smgr_->addCameraSceneNode( smgr_->getRootSceneNode(),
                             core::vector3df( sizeX_/2, sizeY_/2, -sizeZ_ ),
                             core::vector3df( sizeX_/2, sizeY_/2, 0 ) );

  int step = 0;
  // Start the simulation loop
  while ( device_->run() ) {

    // emit particles
    emitParticles();
    if ( povFileName_.length() > 0 ) writePovray( step );

    // Draw all primitives
    drvr_->beginScene(true, true, video::SColor(255,255,255,255));
    smgr_->drawAll();
    drvr_->endScene();

    // Simulate one LBM step
    solver_.runStep();
    // Update the particles
    updateParticles();

    // Update the window caption
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

    // Check for end of simulation
    if ( ++step >= maxSteps_ ) break;
  }
}

inline void ParticleSystem::updateParticles() {

  // Go over all particles of all emitters
  std::vector< Emitter >::iterator ite, ite2;
  std::list< Particle >::iterator itp, itp2;
  for ( ite = emitters_.begin(); ite != emitters_.end(); ++ite ) {
    // Precompute modificator for particle size
    float szCoeff = 1. / ( (*ite).fuel_ * (*ite).lifetimeCoeff_ );
    for ( itp = (*ite).particles_.begin(); itp != (*ite).particles_.end(); ++itp) {

      // Remove particles that have left the domain or exceeded their lifetime
      while ( itp != (*ite).particles_.end() && ( (*itp).lifetime_ < 1 ||
              (*itp).getPos().X < 1 || (*itp).getPos().X > sizeX_ - 1 ||
              (*itp).getPos().Y < 1 || (*itp).getPos().Y > sizeY_ - 1 ||
              (*itp).getPos().Z < 1 || (*itp).getPos().Z > sizeZ_ -1  ||
              (*itp).temp_ < ambTemp_ ) ) {
        // DEBUG output
//          std::cout << "Removing particle " << (*itp).sprites_[0]->getID();
//          std::cout << " at position <" << (*itp).getPos().X << "," << (*itp).getPos().Y << "," << (*itp).getPos().Z << ">";
//          std::cout << " with lifetime " << (*itp).lifetime_ << ", Remaining: " << (*ite).particles_.size() << std::endl;
        (*itp).clear();
        itp = (*ite).particles_.erase( itp );
      }

      // If no particles are left anymore, go to next emitter
      if ( itp == (*ite).particles_.end() ) break;

      // Move particle according to fluid velocity at current position
      Vec3<float> vel = solver_.getVelocity( (*itp).pos_.X, (*itp).pos_.Y, (*itp).pos_.Z );
      core::vector3df v = core::vector3df( vel[0], vel[1], vel[2] );
      // Buoyancy force
      core::vector3df g =  gravity_ * k_ * ( (*itp).temp_ - ambTemp_ );
      // Check whether the buoyancy force would carry the particle into an
      // obstacle
      bool inObstacle = false;
      core::vector3df u = (*itp).pos_ + g + v;
      for ( int i = 0; i < obstacles_.size(); ++i )
        inObstacle &= obstacles_[i].isPointInside( u );
      // If not, apply it
      if ( !inObstacle ) (*itp).updatePos( v + g );
      else (*itp).updatePos( v );
      // DEBUG output
//        std::cout << "Velocity after update: <" << (*itp).pos_[0] << "," << (*itp).pos_[1] << "," << (*itp).pos_[2] << ">" << std::endl;

      if ( (*itp).type_ == FIRE ) {
        // Update temperature
        float tempExt = - gaussTable_[0] * (*itp).temp_;
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
        (*itp).setColor( bbColorTable_[ (int) (((*itp).temp_ - smokeTemp_) / 50.) ] );
        // If temperature has fallen below threshold, convert to smoke particle
        if ( (*itp).temp_ < smokeTemp_ ) {
          (*itp).type_ = SMOKE;
          (*itp).setColor( video::SColor( 25,220,220,220 ) );
        }
      }

      // Update lifetime and particle size
      if ( (*itp).lifetime_-- % 10 == 0 ) {
        float sz = sizeBase_ + sizeVar_ * (*itp).lifetime_ * szCoeff;
        (*itp).setSize( sz );
      }
    }
  }
}

inline void ParticleSystem::emitParticles() {
  // Go over all emitters
  std::vector< Emitter >::iterator ite;
  for ( ite = emitters_.begin(); ite != emitters_.end(); ++ite ) {
    // Check if emitter is to create new particle, depending on emit threshold
    if ( (*ite).fuel_ > std::rand() % (*ite).emitThreshold_ ) {
      (*ite).particles_.push_back( Particle( smgr_,
                                   numParticles_++,
                                   (*ite).pos_,
                                   textures_[ std::rand() % textures_.size() ],
                                   numSprites_,
                                   (*ite).temp_,
                                   bbColorTable_[ (int) (((*ite).temp_ - smokeTemp_) / 50.) ],
                                   (int) ((*ite).fuel_ * (*ite).lifetimeCoeff_ )) );
      // Reduce emitter's fuel
      (*ite).fuel_ *= (*ite).fuelConsumption_;
    }
  }
}

void ParticleSystem::writePovray( int step ) {

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
  std::vector< Emitter >::iterator ite;
  std::list< Particle >::iterator itp;
  for ( ite = emitters_.begin(); ite != emitters_.end(); ++ite ) {
    for ( itp = (*ite).particles_.begin(); itp != (*ite).particles_.end(); ++itp) {
      Particle p = *itp;
      povFile << "disc {\n";
      povFile << "  <" << p.pos_.X << "," << p.pos_.Y << "," << p.pos_.Z << ">, ";
      povFile << "-z, " << p.lifetime_ / 100. << "\n";
      if ( p.type_ == FIRE ) povFile << "  texture { pigment { color rgb<1,0,0> } }\n";
      else povFile << "  texture { pigment { color rgb<.5,.5,.5> } }\n";
      povFile << "}\n";
    }
  }
}

void ParticleSystem::generateBlackBodyColorTable( float maxTemp ) {

  float cie_colour_match[81][3] = {
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

  for ( float t = smokeTemp_; t < maxTemp; t += 50. ) {

    // Calculate x,y,z colors from solar spectrum

    int i;
    float lambda, x = 0, y = 0, z = 0, xyz;
    for (i = 0, lambda = 380; lambda < 780.1; i++, lambda += 5) {
        float Me;
        // Get black body radiation intensity for given temperature and
        // wavelength
        float wlm = lambda * 1e-9; // wavelength in meters
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

    float xr, yr, zr, xg, yg, zg, xb, yb, zb;
    float xw, yw, zw;
    float rx, ry, rz, gx, gy, gz, bx, by, bz;
    float rw, gw, bw;

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
    float r = (rx * x) + (ry * y) + (rz * z);
    float g = (gx * x) + (gy * y) + (gz * z);
    float b = (bx * x) + (by * y) + (bz * z);

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
    float max = (r > b) ? ( (r > g) ? r : g ) : ( (b > g) ? b : g );
    r *= 255. / max;
    g *= 255. / max;
    b *= 255. / max;

    bbColorTable_.push_back( video::SColor( 25, r, g, b ) );
    std::cout << "Temperature " << t << "K: RGB <" << r << ", " << g << ", " << b << ">" << std::endl;
  }
  std::cout << "Generated black body color table from " << smokeTemp_;
  std::cout << "K to " << maxTemp << "K (" << bbColorTable_.size() << " values)" << std::endl;
}

} // namespace particles
