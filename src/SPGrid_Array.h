//#####################################################################
// Class SPGrid_Array
//#####################################################################
#ifndef __SPGrid_Array_h__
#define __SPGrid_Array_h__

// #define SPGRID_CHECK_BOUNDS

#include <SPGrid_Geometry.h>

namespace SPGrid{

template<class T,class T_MASK>
class SPGrid_Array
{
public:
    enum {dim=T_MASK::dim};
    typedef T DATA;
    typedef T_MASK MASK;

private:
    void* const data_ptr;

public:
    const SPGrid_Geometry<dim>& geometry;

    SPGrid_Array(void* const data_ptr_input, const SPGrid_Geometry<dim>& geometry_input)
        :data_ptr(data_ptr_input),geometry(geometry_input)
    {}

    inline T& operator()(const std_array<unsigned int,3>& coord)
    {
        Static_Assert(dim==3);
#ifdef SPGRID_CHECK_BOUNDS
        geometry.Check_Bounds(coord.data[0],coord.data[1],coord.data[2]);
#endif        
        return *reinterpret_cast<T*>(reinterpret_cast<uint64_t>(data_ptr)+T_MASK::Linear_Offset(coord));
    }
    
    inline T& operator()(const std_array<unsigned int,2>& coord)
    {
        Static_Assert(dim==2);
#ifdef SPGRID_CHECK_BOUNDS
        geometry.Check_Bounds(coord.data[0],coord.data[1]);
#endif        
        return *reinterpret_cast<T*>(reinterpret_cast<uint64_t>(data_ptr)+T_MASK::Linear_Offset(coord));
    }

    inline T& operator()(const std_array<int,3>& coord)
    {
        Static_Assert(dim==3);
#ifdef SPGRID_CHECK_BOUNDS
        geometry.Check_Bounds(coord.data[0],coord.data[1],coord.data[2]);
#endif        
        return *reinterpret_cast<T*>(reinterpret_cast<uint64_t>(data_ptr)+T_MASK::Linear_Offset(coord));
    }
    
    inline T& operator()(const std_array<int,2>& coord)
    {
        Static_Assert(dim==2);
#ifdef SPGRID_CHECK_BOUNDS
        geometry.Check_Bounds(coord.data[0],coord.data[1]);
#endif        
        return *reinterpret_cast<T*>(reinterpret_cast<uint64_t>(data_ptr)+T_MASK::Linear_Offset(coord));
    }
    
    inline T& operator()(const unsigned int i,const unsigned int j,const unsigned int k)
    {
        Static_Assert(dim==3);
#ifdef SPGRID_CHECK_BOUNDS
        geometry.Check_Bounds(i,j,k);
#endif        
        return *reinterpret_cast<T*>(reinterpret_cast<uint64_t>(data_ptr)+T_MASK::Linear_Offset(i,j,k));
    }
    
    inline T& operator()(const uint64_t offset)
    {
        return *reinterpret_cast<T*>(reinterpret_cast<uint64_t>(data_ptr)+offset);
    }

// *** New templatized versions
    template<int di, int dj>
    inline T& operator()(const uint64_t offset)
    {
        Static_Assert(dim==2);
        return *reinterpret_cast<T*>(reinterpret_cast<uint64_t>(data_ptr)+T_MASK::Packed_Offset<di,dj>(offset));
    }
    template<int di, int dj, int dk>
    inline T& operator()(const uint64_t offset)
    {
        Static_Assert(dim==3);
        return *reinterpret_cast<T*>(reinterpret_cast<uint64_t>(data_ptr)+T_MASK::Packed_Offset<di,dj,dk>(offset));
    }

    inline T& operator()(const unsigned int i,const unsigned int j)
    {
        Static_Assert(dim==2);
#ifdef SPGRID_CHECK_BOUNDS
        geometry.Check_Bounds(i,j);
#endif        
        return *reinterpret_cast<T*>(reinterpret_cast<uint64_t>(data_ptr)+T_MASK::Linear_Offset(i,j));
    }

// *** Templatized Get's
    template<int di, int dj>
    inline T& Get(uint64_t offset)
    {
        Static_Assert(dim==2);
        return *reinterpret_cast<T*>(reinterpret_cast<uint64_t>(data_ptr)+T_MASK::Packed_Offset<di,dj>(offset));
    }

    template<int di, int dj, int dk>
    inline T& Get(uint64_t offset)
    {
        Static_Assert(dim==3);
        return *reinterpret_cast<T*>(reinterpret_cast<uint64_t>(data_ptr)+T_MASK::Packed_Offset<di,dj,dk>(offset));
    }

    // Debug_Get() functions operate like the operator parenthesis, but also check if the memory address is resident
    T& Debug_Get(const unsigned int i,const unsigned int j,const unsigned int k)
    {
        Static_Assert(dim==3);
#ifdef SPGRID_CHECK_BOUNDS
        geometry.Check_Bounds(i,j,k);
#endif        
        T* addr=reinterpret_cast<T*>(reinterpret_cast<uint64_t>(data_ptr)+T_MASK::Linear_Offset(i,j,k));
        if(!Check_Address_Resident(addr)) FATAL_ERROR("In Check_Address_Resident() : Input address "+Value_To_String(addr)+" is not resident in physical memory");
        return *addr;
    }

    T& Debug_Get(const unsigned int i,const unsigned int j)
    {
        Static_Assert(dim==2);
#ifdef SPGRID_CHECK_BOUNDS
        geometry.Check_Bounds(i,j);
#endif        
        T* addr=reinterpret_cast<T*>(reinterpret_cast<uint64_t>(data_ptr)+T_MASK::Linear_Offset(i,j));
        if(!Check_Address_Resident(addr)) FATAL_ERROR("In Check_Address_Resident() : Input address "+Value_To_String(addr)+" is not resident in physical memory");
        return *addr;
    }

    inline const void* Get_Data_Ptr() const { return data_ptr; }
    inline void* Get_Data_Ptr() { return data_ptr; }

//#####################################################################
};
}

#endif
