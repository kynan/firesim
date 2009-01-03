/*
 * LBM.h
 *
 *  Created on: Jan 2, 2009
 *      Author: Florian Rathgeber
 */

#ifndef LBM_H_
#define LBM_H_

#include "D3Q19.h"

namespace lbm {

class LBM {

public:

  LBM();
  LBM( int sizeX, int sizeY, int sizeZ );
  virtual ~LBM();

  void run( double omega, int maxSteps );

private:

  //! Distribution function fields
  dfField* grid0_;
  dfField* grid1_;

  //! velocity field
  vField u_;

  //! density field
  sField rho_;

  //! Fluid viscosity in range (0..2]
  //double omega_;
};

} // namespace lbm

#endif /* LBM_H_ */
