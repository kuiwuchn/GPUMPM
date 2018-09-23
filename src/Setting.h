#ifndef __SETTING_H_
#define __SETTING_H_

using T = float;
constexpr const int N = 128;
constexpr const T dx = 1. / N;
constexpr const T one_over_dx = N;
constexpr const T D_inverse = 4.f * one_over_dx * one_over_dx;
constexpr int Dim = 3;

// 0: explicit 1: implicit
#define MPM_SIM_TYPE 0

// 0: flip 1: apic 2: mls
#define TRANSFER_SCHEME 1

// 0: 7M cube 1: two dragons collide
#define GEOMETRY_TYPE 1

// the amount of grid being used
#define MEMORY_SCALE 0.6//0.1

#define SAVE_DATA

template<typename T>
struct TEST_STRUCT {
    unsigned flags;
    T ch0;
    T ch1;
    T ch2;
    T ch3;
    T ch4;
    T ch5;
    T ch6;
    T ch7;
    T ch8;
    T ch9;
    T ch10;
    T ch11;
    T ch12;
    T ch13;
    T ch14;
};

enum {
    MPM_DIRICHLET = 0x00000001u
};

#endif
