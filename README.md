DAlitz Fitter aNd Event generator (DAFNE)
=========================================

Dafne is a C++ framework for fitting and generation of time-(in)dependent 3-body decays.

It has external dependencies on Hydra, ROOT and TCLAP. Hydra and ROOT have their own
dependencies. Please refer to the istructions provided in the
[Hydra GitHub repository](https://github.com/MultithreadCorner/Hydra)
and on the [ROOT website](https://root.cern.ch) for more information.

Dafne requires Hydra release 3.2.3 (or higher). Any recent version of ROOT should work with
Dafne, we recommend to use the latest production version available. If TCLAP is not
installed in your system, download it from [sourceforge.net](http://tclap.sourceforge.net).

While ROOT needs to be compiled, Hydra and TCLAP are header-only libraries, so downloading
them is enough for Dafne to work.

You can compile all Dafne's examples by executing the following commands:
			
```
mkdir build
cd build
cmake -D HYDRA_INCLUDE_DIR=path_to_hydra_dir -D TCLAP_INCLUDE_DIR=path_to_tclap_dir ../
make
```

Some example programs may have a ``.cu`` and a ``.cpp`` source file in the ``src`` directory.
These file have the same include statement of a ``.inl`` file containing the actual source
code implementing the program. Cmake determines which backends are available and generates
a Makefile to compile the code accordingly. The compiled executables will be in the
``build/bin`` directory. For each example you will find different executables
compiled against the available backends. The name convention is `example_{cpp,tbb,omp,cu}`
for standard (single-threaded) C++, TBB (multi-threaded CPU), OpenMP (multi-threaded CPU),
CUDA (Nvidia GPU) executables, respectively.			               
