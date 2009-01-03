/*
 * D3Q19.h
 *
 *  Created on: Jan 2, 2009
 *      Author: fuan
 */

#ifndef D3Q19_H_
#define D3Q19_H_

#include "Grid.h"

namespace lbm {

const int dim = 19;

typedef Grid<double,19> dfField;
typedef Grid<double, 3> vField;
typedef Grid<double, 1> sField;

//! Directions of the 19 distribution values for each cell
enum Direction {
  C = 0,
  N = 1,
  E = 2,
  S = 3,
  W = 4,
  T = 5,
  B = 6,
  NE = 7,
  SE = 8,
  SW = 9,
  NW = 10,
  TN = 11,
  TE = 12,
  TS = 13,
  TW = 14,
  BN = 15,
  BE = 16,
  BS = 17,
  BW = 18
};

//! Inverse lattice directions for numerals 0 - 18
const int finv[] = {
    0,
    3,
    4,
    1,
    2,
    6,
    5,
    9,
    10,
    7,
    17,
    18,
    15,
    16,
    13,
    14,
    11,
    12
};

//! Corresponding Direction value for the numerals 0 - 18
const Direction fd[] = {
  C,
  N,
  E,
  S,
  W,
  T,
  B,
  NE,
  SE,
  SW,
  NW,
  TN,
  TE,
  TS,
  TW,
  BN,
  BE,
  BS,
  BW
};

//! Inverse Direction values for numerals 0 - 18
const Direction fdinv[] = {
    C,
    S,
    W,
    N,
    E,
    B,
    T,
    SW,
    NW,
    NE,
    BS,
    BW,
    BN,
    BE,
    TS,
    TW,
    TN,
    TE
};

//! Weights of the distribution values for collision
const double w[] = {
    1./3.,
    1./18.,
    1./18.,
    1./18.,
    1./18.,
    1./18.,
    1./18.,
    1./36.,
    1./36.,
    1./36.,
    1./36.,
    1./36.,
    1./36.,
    1./36.,
    1./36.,
    1./36.,
    1./36.,
    1./36.,
    1./36.
};

//! Lattice velocities
//                 C, N, E, S, W, T, B,NE,SE,SW,NW,TN,TE,TS,TW,BN,BE,BS,BW
const int ex[] = { 0, 0, 1, 0,-1, 0, 0, 1, 1,-1,-1, 0, 1, 0,-1, 0, 1, 0,-1 };
const int ey[] = { 0, 1, 0,-1, 0, 0, 0, 1,-1,-1, 1, 1, 0,-1, 0, 1, 0,-1, 0 };
const int ez[] = { 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1, 1, 1, 1,-1,-1,-1,-1 };

} // namespace lbm

#endif /* D3Q19_H_ */
