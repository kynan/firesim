/*
 * LBM.h
 *
 *  Created on: Jan 2, 2009
 *      Author: Florian Rathgeber
 */

#ifndef LBM_H_
#define LBM_H_

#include "D3Q19.h"
#include "Vec.h"

namespace lbm {

class LBM {

public:

  // Constructors and destructors
  // ============================

  LBM();

  LBM( int sizeX,
       int sizeY,
       int sizeZ,
       std::vector< Vec3<int> > &boundaryCells,
       std::vector< Vec3<int> > &velocityCells,
       std::vector< Vec3<double> > &velocities );

  virtual ~LBM();

  void run( double omega, int maxSteps );

private:

  // Internal helper functions
  // =========================

  inline void collideStream( int x, int y, int z, double omega );

  inline void treatBoundary();

  inline void treatVelocities();

  // Data members
  // ============

  //! Distribution function fields
  dfField* grid0_;
  dfField* grid1_;

  //! velocity field
  vField u_;

  //! density field
  sField rho_;

  //! List with coordinates of all boundary cells
  std::vector< Vec3<int> > boundaryCells_;

  //! List with coordinates of all velocity cells
  std::vector< Vec3<int> > velocityCells_;

  //! List with velocities for all velocity cells
  std::vector< Vec3<double> > velocities_;
};

} // namespace lbm

#endif /* LBM_H_ */
