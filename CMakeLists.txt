cmake_minimum_required(VERSION 3.2)

project(dafne)

if(NOT UNIX)
	message(FATAL_ERROR "This is an unsupported system.")
endif()

set(Project_CMAKE_DIR "${PROJECT_SOURCE_DIR}/cmake")
set(CMAKE_MODULE_PATH "${Project_CMAKE_DIR}" ${CMAKE_MODULE_PATH})
set(CMAKE_VERBOSE_MAKEFILE  ON)

set(EXECUTABLE_OUTPUT_PATH ${CMAKE_BINARY_DIR}/bin)

#check if compiler is C++17 compliant
include(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("--std=c++17" COMPILER_SUPPORTS_CXX17)
if(NOT COMPILER_SUPPORTS_CXX17)
	message(FATAL "The compiler ${CMAKE_CXX_COMPILER} does not support C++17. Please use a different C++ compiler.")
endif()

#compiler flags
set(CMAKE_CXX_FLAGS "--std=c++17 -W -Wno-unused-parameter -march=native -fPIC -O3" CACHE STRING "compile flags" FORCE)

#-----------------------
# get Hydra
find_package(Hydra REQUIRED)
if(HYDRA_FOUND)
	include_directories(${HYDRA_INCLUDE_DIR})
endif(HYDRA_FOUND)

#-----------------------
# get TCLAP
find_package(TCLAP REQUIRED)
if(TCLAP_FOUND)
	include_directories(${TCLAP_INCLUDE_DIR})
endif(TCLAP_FOUND)

#-----------------------
# get ROOT
find_package(ROOT COMPONENTS Minuit2)
if(ROOT_FOUND)
	include_directories(${ROOT_INCLUDE_DIR})
	link_directories(${ROOT_LIBRARY_DIR})
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -D_ROOT_AVAILABLE_")
	if(${ROOT_Minuit2_LIBRARY} MATCHES "libMinuit2.so")
		set(Minuit2_FOUND ON)
	else(Minuit2_FOUND OFF)
	endif()
endif(ROOT_FOUND)

#-----------------------
# get Eigen for K-Matrix examples
find_package( Eigen3 3.3.7 REQUIRED  )
include_directories( EIGEN3_INCLUDE_DIR  )
if(EIGEN3_FOUND)
	include_directories(${EIGEN3_INCLUDE_DIR})
endif(EIGEN3_FOUND)

# CUDA is not supported for the moment
#
# #-----------------------
# #get CUDA
# find_package(CUDA 8.0)
# if(CUDA_FOUND)
# 	link_directories( ${CUDA_TOOLKIT_ROOT_DIR}/lib64/)
# 
# 	#set cuda flags
# 	set(CUDA_NVCC_FLAGS ${CUDA_NVCC_FLAGS}  -ftemplate-backtrace-limit=0; --cudart; static;  -O4 ; --expt-relaxed-constexpr; -fmad=true; --expt-extended-lambda;--relocatable-device-code=false; -Xptxas -dlcm=ca;-Xptxas --opt-level=4 )
# 
# 	set(CUDA_PROPAGATE_HOST_FLAGS ON)
# 	set(CUDA_SEPARABLE_COMPILATION OFF)
# 	set(CUDA_VERBOSE_BUILD OFF)
# 
# 	if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER 5.4)
# 		list(APPEND CUDA_NVCC_FLAGS " -Xcompiler -D__CORRECT_ISO_CPP11_MATH_H_PROTO ")
# 	endif()
# 
# 	# Detect CUDA architecture and get best NVCC flags
# 	include(${CMAKE_CURRENT_SOURCE_DIR}/cmake/FindCudaArch.cmake)
# 	select_nvcc_arch_flags(NVCC_FLAGS_EXTRA)
# 	list(APPEND CUDA_NVCC_FLAGS ${NVCC_FLAGS_EXTRA})
# endif(CUDA_FOUND)

#-----------------------
#get TBB
find_package(TBB )
if(TBB_FOUND)
	include_directories(${TBB_INCLUDE_DIRS})
	link_directories(${TBB_LIBRARY})
endif(TBB_FOUND)

#-----------------------
#get OpenMP
find_package(OpenMP)
if(OPENMP_CXX_FOUND OR OPENMP_FOUND)
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif(OPENMP_CXX_FOUND OR OPENMP_FOUND)

#-----------------------
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

#always on
set(BUILD_CPP_TARGETS  ON)

#crucial for clang build, where openmp support is not embedded in the compiler
if( OPENMP_CXX_FOUND OR OPENMP_FOUND)
	set(BUILD_OMP_TARGETS  ON)
else()
	set(BUILD_OMP_TARGETS  OFF)
endif()

# CUDA is not supported for the moment
#
# #not every computer has a GPU and/or cuda installed or 
# if( CUDA_FOUND )
# 	set(BUILD_CUDA_TARGETS  ON)
# else()
# 	set(BUILD_CUDA_TARGETS  OFF)
# endif()

# to enable tbb builds
if(TBB_FOUND)
	set(BUILD_TBB_TARGETS  ON)
else()
	set(BUILD_TBB_TARGETS OFF)
endif(TBB_FOUND)

# messages
message(STATUS "Hydra include path: ${HYDRA_INCLUDE_DIR}")
message(STATUS "TCLAP include path: ${TCLAP_INCLUDE_DIR}")
message(STATUS "ROOT include path: ${ROOT_INCLUDE_DIR}")
message(STATUS "ROOT library path: ${ROOT_LIBRARY_DIR}" )
message(STATUS "ROOT libraries: ${ROOT_LIBRARIES}")
message(STATUS "ROOT::Minuit2 library: ${ROOT_Minuit2_LIBRARY}")
message(STATUS "Build CUDA/NVCC-based targets: ${BUILD_CUDA_TARGETS}")
message(STATUS "CUDA include: ${CUDA_INCLUDE_DIRS}" )
message(STATUS "CUDA libraries: ${CUDA_LIBRARIES}"  )
message(STATUS "NVCC flags: ${CUDA_NVCC_FLAGS}" )
message(STATUS "Build TBB targets: ${BUILD_TBB_TARGETS}")
message(STATUS "TBB include: ${TBB_INCLUDE_DIRS}" )
message(STATUS "TBB libraries: ${TBB_LIBRARIES}" )
message(STATUS "Build OpenMP targets: ${BUILD_OMP_TARGETS}")
message(STATUS "OpenMP flags: ${OpenMP_CXX_FLAGS}" )
message(STATUS "OpenMP libraries:  ${OpenMP_CXX_LIBRARIES}" )
message(STATUS "Host compiler flags: ${CMAKE_CXX_FLAGS}" )
message(STATUS "Install directory: ${CMAKE_INSTALL_PREFIX}"  )
message(STATUS "Dafne source dir: ${PROJECT_SOURCE_DIR}" )
message(STATUS "Dafne build dir: ${PROJECT_BINARY_DIR}" )
message(STATUS "Dafne binary dir: ${EXECUTABLE_OUTPUT_PATH}" )

#+++++++++++++++++++++++++++
#   targets                +
#+++++++++++++++++++++++++++
include(${Project_CMAKE_DIR}/AddAnalysisTarget.cmake)

ADD_ANALYSIS(dalitz-fit)
ADD_ANALYSIS(time-dependent-dalitz-fit)
ADD_ANALYSIS(generate-Kpipi0)
ADD_ANALYSIS(generate-KSpipi-time-dependent)
ADD_ANALYSIS(generate-KSpipi-time-dependent-time-resolution)
ADD_ANALYSIS(generate-KSpipi-time-dependent-time-resolution-with-background)
ADD_ANALYSIS(generate-KSpipi-time-independent)
ADD_ANALYSIS(generate-KSpipi-time-independent-with-background)
ADD_ANALYSIS(fit-KSpipi-time-independent)
ADD_ANALYSIS(fit-KSpipi-time-independent-with-background)
ADD_ANALYSIS(fit-KSpipi-time-dependent)
ADD_ANALYSIS(fit-KSpipi-time-dependent-time-resolution)
