#	Required version of CMake
cmake_minimum_required(VERSION 3.0)

#	Project title
project(MPM)

include_directories(../ ${PROJECT_SOURCE_DIR}/)

#	C++


cmake_host_system_information(RESULT _host_name QUERY HOSTNAME)
message("-- _host_name variable is set to: " ${_host_name})
if(${_host_name} STREQUAL "mine-win10")

set(CMAKE_CXX_FLAGS "/std:c++14 /EHsc")
include_directories(
	"C:/SDKs/cuda/v9.2/Toolkit/include"
	"C:/SDKs/cuda/v9.2/Samples/common/inc")
link_directories("C:/SDKs/cuda/v9.2/Toolkit/lib/x64")

include_directories(C:/SDKs/Partio)   
link_directories(C:/SDKs/Partio)

else()

set(CMAKE_C_COMPILER /usr/bin/gcc-7)
set(CMAKE_CXX_COMPILER /usr/bin/g++-7)
set(CMAKE_CXX_FLAGS "-std=c++14 -O3 -c")

include_directories(D:/Libs/CUDAToolkit/v8.0/include)
link_directories(D:/Libs/CUDAToolkit/v8.0/lib/x64)

include_directories(D:/Libs/partio)   
link_directories(D:/Libs/partio)    
endif()

add_library(mpm_cpp
	../Simulation/Material/Material.cpp
	../Simulation/Particle/ParticleDomain.cpp
	../Simulation/Grid/GridDomain.cpp
	../Simulation/DomainTransform/DomainTransformer.cpp
	../Simulation/MPM/Simulator.cpp
	../Simulation/MPM/SimulatorBuilder.cpp
	../SPGrid_Utilities.cpp
	../SPGrid_Geometry.cpp
	../MnBase/Math/Probability/Probability.cpp
)

target_link_libraries(mpm_cpp
                      cudart
                      cudadevrt
                      partio
                     )
