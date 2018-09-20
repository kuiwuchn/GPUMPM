#include "Benchmarks.cuh"
#include <System/CudaDevice/CudaKernelLauncher.cu>
#include <System/CudaDevice/CudaDeviceUtils.cuh>

namespace mn {

    Benchmarks::Benchmarks(MPMSimulator& simulator) 
        : _kSimulator(simulator) {}

    Benchmarks::~Benchmarks() {
    }

    void Benchmarks::run() {
        SimulatorBuilder builder;
        builder.build(_kSimulator, GEOMETRY_TYPE);
        _kSimulator.simulateToFrame(60);
    }

}
