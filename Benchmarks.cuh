#ifndef __BENCHMARKS_CUH_
#define __BENCHMARKS_CUH_
#include <array>
#include <cuda_runtime.h>
#include <Setting.h>
#include <Simulation/MPM/SimulatorBuilder.h>

namespace mn {

    class Benchmarks {
    public:
        //Benchmarks(Particles& p, SPGrid& g, DomainTransformer<TEST_STRUCT<T>>& t, uint64_t* _masks, T* _sigma, T* _massRef, uint64_t* _hmasks);
        Benchmarks(MPMSimulator& simulator);
        ~Benchmarks();
        void run();

    public:
        MPMSimulator&   _kSimulator;
    };

}

#endif
