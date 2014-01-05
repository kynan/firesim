Firesim
=======

Lattice Boltzmann Method (LBM) fluid solver.

Uses the BGK collision model and Smagorinsky turbulence correction.

LBM is used to drive a particle engine for the simulation of fire.

[Irrlicht] is used for the real-time 3D OpenGL visualization.

Dependencies
------------

* The [CMake] build system
* The [Boost] C++ libraries, in particular libbost-regex
* The [Irrlicht] 3D engine
* [Doxygen] to generate the documentation (optional)

Installation
------------

Run the setup script to build a `Release` or `Debug` build:
```
./setup.sh <Build type>
```

Documentation
-------------

* [Documentation generated with Doxygen](http://kynan.github.io/firesim)
* My [BSc thesis] this code was developed for

License
-------

MIT

[Boost]: http://boost.org
[CMake]: http://cmake.org
[Doxygen]: http://doxygen.org
[Irrlicht]: http://irrlicht.sourceforge.net
[BSc thesis]: http://www10.informatik.uni-erlangen.de/Publications/Theses/2009/Rathgeber_BA09.pdf
