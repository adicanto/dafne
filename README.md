DAlitz Fitter aNd Event generator (DAFNE)
=========================================

Dafne is a C++ framework for fitting and generation of time-(in)dependent 3-body decays.

It has external dependencies on [Hydra](https://github.com/MultithreadCorner/Hydra), [ROOT](https://root.cern.ch), [Eigen3](https://gitlab.com/libeigen/eigen/) and  [TCLAP](http://tclap.sourceforge.net). Hydra and ROOT have their own
dependencies. Please refer to the istructions provided in the 
[Hydra GitHub repository](https://github.com/MultithreadCorner/Hydra)
and on the [ROOT website](https://root.cern.ch) for more information.  

Dafne requires Hydra release 4.0 (or higher). Any recent version of ROOT should work with
Dafne, we recommend to use the latest production version available.  

Dafne is a header-only library. Not installation procedure is required.
You can compile all Dafne's examples by executing the following commands:
			

```
mkdir build
cd build
cmake -D CMAKE_C_COMPILER=`which gcc` -D CMAKE_CXX_COMPILER=`which g++` -D HYDRA_INCLUDE_DIR=<path to hydra> -D TCLAP_INCLUDE_PATH=<path to tclap>/include -D EIGEN3_INCLUDE_DIR=<path to eigen3>/include/eigen3 -D TBB_INCLUDE_DIRS=<path to tbb>/include/  -DTBB_tbb_LIBRARY_RELEASE=<path to tbb>/lib64/libtbb.so ../
make -j 10 # if you have 10 cores available for compilation
```
