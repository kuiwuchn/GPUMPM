#ifndef __SPGrid_Utilities_h__
#define __SPGrid_Utilities_h__

#define __attribute__(x)
#include <sstream>

#ifdef _WIN32
#include <algorithm>
#endif

//#define HASWELL

namespace SPGrid{

template<unsigned d> struct BitLength;
template<> struct BitLength<0> {enum {value=0};};
template<unsigned d> struct BitLength {enum {value=1+BitLength<(d>>1)>::value};};
template<unsigned d> struct NextLogTwo {enum {value=BitLength<d-1>::value};};

template<class T,class T_FIELD> size_t
OffsetOfMember(T_FIELD T::*field)
{return (size_t)((char*)&(((T*)0)->*field)-(char*)0);}

template<class T,class T1,class T2> struct EnableIfSame {};
template<class T,class T1> struct EnableIfSame<T,T1,T1> {typedef T type;};

inline unsigned Next_Power_Of_Two(unsigned i)
{i--;i|=i>>1;i|=i>>2;i|=i>>4;i|=i>>8;i|=i>>16;return i+1;}

template<class T> inline T max(const T& x1, const T& x2)
{return std::max(x1,x2);}

template<class T> inline T max(const T& x1, const T& x2,const T& x3)
{return std::max(x1,std::max(x2,x3));}

template<class T> inline T max(const T& x1, const T& x2,const T& x3, const T& x4)
{return std::max(std::max(x1,x2),std::max(x3,x4));}

template<class T> std::string Value_To_String(const T& value)
{std::ostringstream output;output<<value;return output.str();}

#ifdef HASWELL
unsigned long Bit_Spread(const unsigned long data,const unsigned long mask);
unsigned long Bit_Spread(const unsigned int data,const unsigned long mask);
unsigned long Bit_Spread(const int data,const unsigned long mask);
#else
template<uint64_t mask, class T_INT >
uint64_t Bit_Spread(const T_INT data)
{
    uint64_t uldata=data;
    uint64_t result=0;

    if(0x0000000000000001UL & mask) result |= uldata & 0x0000000000000001UL; else uldata <<= 1;
    if(0x0000000000000002UL & mask) result |= uldata & 0x0000000000000002UL; else uldata <<= 1;
    if(0x0000000000000004UL & mask) result |= uldata & 0x0000000000000004UL; else uldata <<= 1;
    if(0x0000000000000008UL & mask) result |= uldata & 0x0000000000000008UL; else uldata <<= 1;
    if(0x0000000000000010UL & mask) result |= uldata & 0x0000000000000010UL; else uldata <<= 1;
    if(0x0000000000000020UL & mask) result |= uldata & 0x0000000000000020UL; else uldata <<= 1;
    if(0x0000000000000040UL & mask) result |= uldata & 0x0000000000000040UL; else uldata <<= 1;
    if(0x0000000000000080UL & mask) result |= uldata & 0x0000000000000080UL; else uldata <<= 1;
    if(0x0000000000000100UL & mask) result |= uldata & 0x0000000000000100UL; else uldata <<= 1;
    if(0x0000000000000200UL & mask) result |= uldata & 0x0000000000000200UL; else uldata <<= 1;
    if(0x0000000000000400UL & mask) result |= uldata & 0x0000000000000400UL; else uldata <<= 1;
    if(0x0000000000000800UL & mask) result |= uldata & 0x0000000000000800UL; else uldata <<= 1;
    if(0x0000000000001000UL & mask) result |= uldata & 0x0000000000001000UL; else uldata <<= 1;
    if(0x0000000000002000UL & mask) result |= uldata & 0x0000000000002000UL; else uldata <<= 1;
    if(0x0000000000004000UL & mask) result |= uldata & 0x0000000000004000UL; else uldata <<= 1;
    if(0x0000000000008000UL & mask) result |= uldata & 0x0000000000008000UL; else uldata <<= 1;
    if(0x0000000000010000UL & mask) result |= uldata & 0x0000000000010000UL; else uldata <<= 1;
    if(0x0000000000020000UL & mask) result |= uldata & 0x0000000000020000UL; else uldata <<= 1;
    if(0x0000000000040000UL & mask) result |= uldata & 0x0000000000040000UL; else uldata <<= 1;
    if(0x0000000000080000UL & mask) result |= uldata & 0x0000000000080000UL; else uldata <<= 1;
    if(0x0000000000100000UL & mask) result |= uldata & 0x0000000000100000UL; else uldata <<= 1;
    if(0x0000000000200000UL & mask) result |= uldata & 0x0000000000200000UL; else uldata <<= 1;
    if(0x0000000000400000UL & mask) result |= uldata & 0x0000000000400000UL; else uldata <<= 1;
    if(0x0000000000800000UL & mask) result |= uldata & 0x0000000000800000UL; else uldata <<= 1;
    if(0x0000000001000000UL & mask) result |= uldata & 0x0000000001000000UL; else uldata <<= 1;
    if(0x0000000002000000UL & mask) result |= uldata & 0x0000000002000000UL; else uldata <<= 1;
    if(0x0000000004000000UL & mask) result |= uldata & 0x0000000004000000UL; else uldata <<= 1;
    if(0x0000000008000000UL & mask) result |= uldata & 0x0000000008000000UL; else uldata <<= 1;
    if(0x0000000010000000UL & mask) result |= uldata & 0x0000000010000000UL; else uldata <<= 1;
    if(0x0000000020000000UL & mask) result |= uldata & 0x0000000020000000UL; else uldata <<= 1;
    if(0x0000000040000000UL & mask) result |= uldata & 0x0000000040000000UL; else uldata <<= 1;
    if(0x0000000080000000UL & mask) result |= uldata & 0x0000000080000000UL; else uldata <<= 1;
    if(0x0000000100000000UL & mask) result |= uldata & 0x0000000100000000UL; else uldata <<= 1;
    if(0x0000000200000000UL & mask) result |= uldata & 0x0000000200000000UL; else uldata <<= 1;
    if(0x0000000400000000UL & mask) result |= uldata & 0x0000000400000000UL; else uldata <<= 1;
    if(0x0000000800000000UL & mask) result |= uldata & 0x0000000800000000UL; else uldata <<= 1;
    if(0x0000001000000000UL & mask) result |= uldata & 0x0000001000000000UL; else uldata <<= 1;
    if(0x0000002000000000UL & mask) result |= uldata & 0x0000002000000000UL; else uldata <<= 1;
    if(0x0000004000000000UL & mask) result |= uldata & 0x0000004000000000UL; else uldata <<= 1;
    if(0x0000008000000000UL & mask) result |= uldata & 0x0000008000000000UL; else uldata <<= 1;
    if(0x0000010000000000UL & mask) result |= uldata & 0x0000010000000000UL; else uldata <<= 1;
    if(0x0000020000000000UL & mask) result |= uldata & 0x0000020000000000UL; else uldata <<= 1;
    if(0x0000040000000000UL & mask) result |= uldata & 0x0000040000000000UL; else uldata <<= 1;
    if(0x0000080000000000UL & mask) result |= uldata & 0x0000080000000000UL; else uldata <<= 1;
    if(0x0000100000000000UL & mask) result |= uldata & 0x0000100000000000UL; else uldata <<= 1;
    if(0x0000200000000000UL & mask) result |= uldata & 0x0000200000000000UL; else uldata <<= 1;
    if(0x0000400000000000UL & mask) result |= uldata & 0x0000400000000000UL; else uldata <<= 1;
    if(0x0000800000000000UL & mask) result |= uldata & 0x0000800000000000UL; else uldata <<= 1;
    if(0x0001000000000000UL & mask) result |= uldata & 0x0001000000000000UL; else uldata <<= 1;
    if(0x0002000000000000UL & mask) result |= uldata & 0x0002000000000000UL; else uldata <<= 1;
    if(0x0004000000000000UL & mask) result |= uldata & 0x0004000000000000UL; else uldata <<= 1;
    if(0x0008000000000000UL & mask) result |= uldata & 0x0008000000000000UL; else uldata <<= 1;
    if(0x0010000000000000UL & mask) result |= uldata & 0x0010000000000000UL; else uldata <<= 1;
    if(0x0020000000000000UL & mask) result |= uldata & 0x0020000000000000UL; else uldata <<= 1;
    if(0x0040000000000000UL & mask) result |= uldata & 0x0040000000000000UL; else uldata <<= 1;
    if(0x0080000000000000UL & mask) result |= uldata & 0x0080000000000000UL; else uldata <<= 1;
    if(0x0100000000000000UL & mask) result |= uldata & 0x0100000000000000UL; else uldata <<= 1;
    if(0x0200000000000000UL & mask) result |= uldata & 0x0200000000000000UL; else uldata <<= 1;
    if(0x0400000000000000UL & mask) result |= uldata & 0x0400000000000000UL; else uldata <<= 1;
    if(0x0800000000000000UL & mask) result |= uldata & 0x0800000000000000UL; else uldata <<= 1;
    if(0x1000000000000000UL & mask) result |= uldata & 0x1000000000000000UL; else uldata <<= 1;
    if(0x2000000000000000UL & mask) result |= uldata & 0x2000000000000000UL; else uldata <<= 1;
    if(0x4000000000000000UL & mask) result |= uldata & 0x4000000000000000UL; else uldata <<= 1;
    if(0x8000000000000000UL & mask) result |= uldata & 0x8000000000000000UL; else uldata <<= 1;

    return result;
}
#endif

int Bit_Pack(const uint64_t data, const uint64_t mask);

template<uint64_t data,uint64_t mask,uint64_t bit=1UL> struct BitSpread;
template<uint64_t data,uint64_t mask> struct BitSpread<data,mask,0> { static const uint64_t value=0; };
template<uint64_t data,uint64_t mask,uint64_t bit> struct BitSpread
{ static const uint64_t value = (bit & mask ) ? BitSpread<data,mask,bit<<1>::value | (data & bit) : BitSpread<data<<1,mask,bit<<1>::value; };

#define FATAL_ERROR(...) \
    Fatal_Error((const char*)__FUNCTION__,__FILE__,__LINE__,##__VA_ARGS__)

void Fatal_Error(const char* function,const char* file,unsigned int line) __attribute__ ((noreturn)) ;
void Fatal_Error(const char* function,const char* file,unsigned int line,const char* message) __attribute__ ((noreturn)) ;
void Fatal_Error(const char* function,const char* file,unsigned int line,const std::string& message) __attribute__ ((noreturn)) ;

#define SPGRID_JOIN( X, Y ) SPGRID_DO_JOIN( X, Y )
#define SPGRID_DO_JOIN( X, Y ) SPGRID_DO_JOIN2(X,Y)
#define SPGRID_DO_JOIN2( X, Y ) X##Y

template <bool x> struct STATIC_ASSERTION_FAILURE;
template <> struct STATIC_ASSERTION_FAILURE<true> { enum { value = 1 }; };
template<int x> struct static_assert_test{};
#define SPGRID_STATIC_ASSERT( B ) \
   typedef static_assert_test<\
      sizeof(STATIC_ASSERTION_FAILURE< (bool)( B ) >)>\
         SPGRID_JOIN(boost_static_assert_typedef_, __LINE__)
#define Static_Assert(...) SPGRID_STATIC_ASSERT((__VA_ARGS__))

}
#endif
