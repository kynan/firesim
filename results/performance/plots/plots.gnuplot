set terminal postscript eps enhanced color size 3.2,3.2

set output 'lbm_mlup.eps'
load "lbm_mlup.plot"
reset

#pause -1

set output 'lbm_time.eps'
load "lbm_time.plot"
reset

#pause -1

set output 'spheres_time.eps'
load "spheres_time.plot"
reset

#pause -1

set output 'spheres_mlup.eps'
load "spheres_mlup.plot"
reset

#pause -1

#set output 'fmlup.eps'
#load "fmlup.plot"
#reset

#pause -1

#set output 'lbm_share.histogram.eps'
#load "lbm_share.histogram.plot"
#reset

#pause -1

#set output 'lbm_share.filledcurve.nosmago.eps'
#load "lbm_share.filledcurve.plot"
#reset

#pause -1

#load "particle_time.plot"
#reset

#pause -1

#load "particle_share.plot"
#reset

#pause -1

load "irrlicht_time.plot"
reset

#pause -1

load "irrlicht_share.plot"
reset

#pause -1

set terminal postscript eps enhanced color size 6,3.6

set output 'update_render.eps'
load "update_render.plot"
reset

#pause -1

set output 'lighting_time.eps'
load "lighting_time.plot"
reset
