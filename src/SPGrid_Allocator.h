//#####################################################################
#ifndef __SPGrid_Allocator_h__
#define __SPGrid_Allocator_h__

#include <SPGrid_Mask.h>
#include <SPGrid_Array.h>

namespace SPGrid{

//#####################################################################
// Class SPGrid_Allocator
//#####################################################################
template<class T,int dim>
class SPGrid_Allocator: public SPGrid_Geometry<dim>, public SPGrid_Mask_base<NextLogTwo<sizeof(T)>::value,dim>
{
    typedef SPGrid_Geometry<dim> T_Geometry_base;
    typedef SPGrid_Mask_base<NextLogTwo<sizeof(T)>::value,dim> T_Mask_base;
    using T_Mask_base::block_bits;using T_Mask_base::block_xbits;using T_Mask_base::block_ybits;using T_Mask_base::elements_per_block;
    using T_Geometry_base::Padded_Volume;

private:
    void* data_ptr;

public:
    template<class T_FIELD=T> struct Array
    {
      typedef SPGrid_Mask<NextLogTwo<sizeof(T)>::value,NextLogTwo<sizeof(T_FIELD)>::value,dim> mask;
      typedef SPGrid_Array<T_FIELD,mask> type;
    };

    template<class T1,class T_FIELD> typename EnableIfSame<typename Array<T_FIELD>::type,T1,T>::type
    Get_Array(T_FIELD T1::* field)
    {
        typedef typename Array<T_FIELD>::type Array_type;
        size_t offset=OffsetOfMember(field)<<block_bits;
        void* offset_ptr=reinterpret_cast<void*>(reinterpret_cast<uint64_t>(data_ptr)+offset);
        return Array_type(offset_ptr,*this);
    }

    template<class T1,class T_FIELD> typename EnableIfSame<typename Array<const T_FIELD>::type,T1,T>::type
    Get_Array(T_FIELD T1::* field) const
    {
        typedef typename Array<const T_FIELD>::type Array_type;
        size_t offset=OffsetOfMember(field)<<block_bits;
        void* offset_ptr=reinterpret_cast<void*>(reinterpret_cast<uint64_t>(data_ptr)+offset);
        return Array_type(offset_ptr,*this);
    }

    template<class T1,class T_FIELD> typename EnableIfSame<typename Array<const T_FIELD>::type,T1,T>::type
    Get_Const_Array(T_FIELD T1::* field) const
    {
        typedef typename Array<const T_FIELD>::type Array_type;
        size_t offset=OffsetOfMember(field)<<block_bits;
        void* offset_ptr=reinterpret_cast<void*>(reinterpret_cast<uint64_t>(data_ptr)+offset);
        return Array_type(offset_ptr,*this);
    }

    typename Array<>::type Get_Array()
    {
        typedef typename Array<>::type Array_type;
        return Array_type(data_ptr,*this);
    }

    typename Array<const T>::type Get_Array() const
    {
        typedef typename Array<const T>::type Array_type;
        return Array_type(data_ptr,*this);
    }

    typename Array<const T>::type Get_Const_Array() const
    {
        typedef typename Array<const T>::type Array_type;
        return Array_type(data_ptr,*this);
    }

protected:    
    inline void* Get_Data_Ptr() {return data_ptr;}

//#####################################################################
};
}
#endif
