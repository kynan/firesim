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

const int Dim = 19;

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
const int finv[19] = {
    0,  // C   0 -> C
    3,  // N   1 -> S
    4,  // E   2 -> W
    1,  // S   3 -> N
    2,  // W   4 -> E
    6,  // T   5 -> B
    5,  // B   6 -> T
    9,  // NE  7 -> SW
    10, // SE  8 -> NW
    7,  // SW  9 -> NE
    8,  // NW 10 -> SE
    17, // TN 11 -> BS
    18, // TE 12 -> BW
    15, // TS 13 -> BN
    16, // TW 14 -> BE
    13, // BN 15 -> TS
    14, // BE 16 -> TW
    11, // BS 17 -> TN
    12  // BW 18 -> TE
};

//! Corresponding Direction value for the numerals 0 - 18
const Direction fd[19] = {
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
const Direction fdinv[19] = {
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
const double w[19] = {
    1. / 3.,
    1. / 18.,
    1. / 18.,
    1. / 18.,
    1. / 18.,
    1. / 18.,
    1. / 18.,
    1. / 36.,
    1. / 36.,
    1. / 36.,
    1. / 36.,
    1. / 36.,
    1. / 36.,
    1. / 36.,
    1. / 36.,
    1. / 36.,
    1. / 36.,
    1. / 36.,
    1. / 36.
};

//! Lattice velocities
//                 C, N, E, S, W, T, B,NE,SE,SW,NW,TN,TE,TS,TW,BN,BE,BS,BW
const int ex[19] = { 0, 0, 1, 0,-1, 0, 0, 1, 1,-1,-1, 0, 1, 0,-1, 0, 1, 0,-1 };
const int ey[19] = { 0, 1, 0,-1, 0, 0, 0, 1,-1,-1, 1, 1, 0,-1, 0, 1, 0,-1, 0 };
const int ez[19] = { 0, 0, 0, 0, 0, 1,-1, 0, 0, 0, 0, 1, 1, 1, 1,-1,-1,-1,-1 };

//! Products of lattice velocities
const int ep[6][19] = {
//    C, N, E, S, W, T, B,NE,SE,SW,NW,TN,TE,TS,TW,BN,BE,BS,BW
    { 0, 0, 1, 0, 1, 0, 0, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0, 1 }, // ex * ex
    { 0, 1, 0, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 1, 0, 1, 0, 1, 0 }, // ey * ey
    { 0, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1 }, // ez * ez
    { 0, 0, 0, 0, 0, 0, 0, 1,-1, 1,-1, 0, 0, 0, 0, 0, 0, 0, 0 }, // ex * ey
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,-1, 0,-1, 0, 1 }, // ex * ez
    { 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0,-1, 0,-1, 0, 1, 0 }  // ey * ez
};

} // namespace lbm

#endif /* D3Q19_H_ */
