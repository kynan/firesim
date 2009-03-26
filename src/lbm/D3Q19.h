//! \file D3Q19.h
//! Specification and definitions for the D3Q19 model of the Lattice Boltzmann
//! method
//! \date   Jan 2, 2009
//! \author Florian Rathgeber

#ifndef D3Q19_H_
#define D3Q19_H_

//! Common namespace for all LBM classes

namespace lbm {

//! Number of distribution functions for each cell

const int Dim = 19;

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

//! Inverse lattice directions corresponding to numeric directions 0 - 18

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

//! Corresponding Direction value for numeric directions 0 - 18

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

//! Inverse Direction values for numeric directions 0 - 18

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

//! Weights of the distribution values for the collision step

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

//! Lattice velocities in x-direction

//                   C, N, E, S, W, T, B,NE,SE,SW,NW,TN,TE,TS,TW,BN,BE,BS,BW
const int ex[19] = { 0, 0, 1, 0,-1, 0, 0, 1, 1,-1,-1, 0, 1, 0,-1, 0, 1, 0,-1 };

//! Lattice velocities in y-direction

//                   C, N, E, S, W, T, B,NE,SE,SW,NW,TN,TE,TS,TW,BN,BE,BS,BW
const int ey[19] = { 0, 1, 0,-1, 0, 0, 0, 1,-1,-1, 1, 1, 0,-1, 0, 1, 0,-1, 0 };

//! Lattice velocities in z-direction

//                   C, N, E, S, W, T, B,NE,SE,SW,NW,TN,TE,TS,TW,BN,BE,BS,BW
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

//! Normalized lattice velocities in x-direction

//                     C        , N        , E        , S        ,  W,  T,  B,
const Real exn[19] = { 0        , 0        , 1        , 0        , -1,  0,  0,
//                     NE       , SE       , SW       , NW       ,
                       1/sqrt(2), 1/sqrt(2),-1/sqrt(2),-1/sqrt(2),
//                     TN       , TE       , TS       , TW       ,
                       0        , 1/sqrt(2), 0        ,-1/sqrt(2),
//                     BN       , BE       , BS       , BW
                       0        , 1/sqrt(2), 0        ,-1/sqrt(2) };

//! Normalized lattice velocities in y-direction

//                     C, N, E, S, W, T, B, NE       , SE       , SW       , NW       , TN       , TE       , TS       , TW       , BN       , BE       , BS       , BW
const Real eyn[19] = { 0, 1, 0,-1, 0, 0, 0, 1/sqrt(2),-1/sqrt(2),-1/sqrt(2), 1/sqrt(2), 1/sqrt(2), 0        ,-1/sqrt(2), 0        , 1/sqrt(2), 0        ,-1/sqrt(2), 0 };

//! Normalized lattice velocities in z-direction

//                     C, N, E, S, W, T, B, NE       , SE       , SW       , NW       , TN       , TE       , TS       , TW       , BN       , BE       , BS       , BW
const Real ezn[19] = { 0, 0, 0, 0, 0, 1,-1, 0        , 0        , 0        , 0        , 1/sqrt(2), 1/sqrt(2), 1/sqrt(2), 1/sqrt(2),-1/sqrt(2),-1/sqrt(2),-1/sqrt(2),-1/sqrt(2) };

//! Lengths of the lattice links
//                     C, N, E, S, W, T, B, NE    , SE    , SW    , NW    , TN    , TE    , TS    , TW    , BN    , BE    , BS    , BW
const Real le[19] = { 1.,1.,1.,1.,1.,1.,1.,sqrt(2),sqrt(2),sqrt(2),sqrt(2),sqrt(2),sqrt(2),sqrt(2),sqrt(2),sqrt(2),sqrt(2),sqrt(2),sqrt(2) };

} // namespace lbm

#endif /* D3Q19_H_ */
