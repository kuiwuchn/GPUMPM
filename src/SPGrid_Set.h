//#####################################################################
#ifndef __SPGrid_Set_h__
#define __SPGrid_Set_h__

#include <stdio.h>
#include <string.h>
#include <sys/mman.h>
#include <vector>
#include <pthread.h>

#include <std_array.h>

namespace SPGrid{

//#####################################################################
// Class SPGrid_Set
//#####################################################################
template<class T_ARRAY>
class SPGrid_Set
{
    enum {dim=T_ARRAY::dim};
    typedef typename T_ARRAY::DATA T;
    pthread_mutex_t pm_mutex;

public:
    typedef typename T_ARRAY::MASK T_MASK;
    T_ARRAY array;
    uint64_t* page_mask_array;
    uint64_t array_length;
    uint64_t max_linear_offset; // TODO: Change semantics to make this the first offset that is *invalid*
    std::vector<uint64_t> block_offsets;
    bool dirty; // Indicates that block offsets are inconsistent with the bitmap (perhaps for a good reason, if only one of them is used)

    // Single size argument constructor    
    SPGrid_Set(T_ARRAY array_input)
        :array(array_input)
    {
        uint64_t number_of_elements = array.geometry.Padded_Volume();
        uint64_t number_of_pages = number_of_elements>>T_MASK::block_bits;
        array_length = (number_of_pages+0x3fUL)>>6;
        page_mask_array=new uint64_t[array_length]; // TODO: Why "new" and not Raw_Allocate() ?
        memset(reinterpret_cast<void*>(page_mask_array),0,array_length*sizeof(uint64_t)); // TODO: Is this really needed?
        max_linear_offset = array.geometry.Padded_Volume()*T_ARRAY::MASK::Bytes_Per_Element(); // TODO: Check this is correct
        pthread_mutex_init(&pm_mutex,0); // TODO: Check to see where these mutexes are really used
        dirty = false; // They start as consistent -- both empty
    }

    ~SPGrid_Set()
    {
        delete page_mask_array;
    }
    
    inline bool CheckBounds(uint64_t linear_offset)
    {
        return (linear_offset < max_linear_offset);
    }

    void MarkPageActive(uint64_t linear_offset)
    {
        if( linear_offset < max_linear_offset )
        {
            uint64_t page_mask = 1UL << (linear_offset>>12 & 0x3f);
            page_mask_array[linear_offset>>18] |= page_mask;
            dirty = true;
        } 
        else 
            FATAL_ERROR("Linear offset "+Value_To_String(linear_offset)+" is out of range (upper limit = "+Value_To_String(max_linear_offset)+")");
    }

    bool IsPageActive(const uint64_t linear_offset) const
    {
        if( linear_offset < max_linear_offset )
        {
            uint64_t page_mask = 1UL << (linear_offset>>12 & 0x3f);
            return page_mask_array[linear_offset>>18] & page_mask;
        }else
            return false;
    }

    bool IsPageActive(const std_array<unsigned int,dim> coord) const
    {
        uint64_t linear_offset = T_MASK::Linear_Offset(coord);
        if( linear_offset < max_linear_offset )
        {
            uint64_t page_mask = 1UL << (linear_offset>>12 & 0x3f);
            return page_mask_array[linear_offset>>18] & page_mask;
        }else
            return false;
    }

    void Mask(const std_array<unsigned int,dim> coord,const T mask)
    {
        T& data=array(coord);
        uint64_t linear_offset=reinterpret_cast<uint64_t>(&data)-reinterpret_cast<uint64_t>(array.Get_Data_Ptr());
        uint64_t page_mask = 1UL << (linear_offset>>12 & 0x3f);
        assert((linear_offset>>18) < array_length);
        page_mask_array[linear_offset>>18] |= page_mask; dirty = true;
        data |= mask;

        // TODO: Examine using the following syntax instead
        // uint64_t linear_offset=T_ARRAY::MASK::Linear_Offset(coord);
        // MarkPageActive(linear_offset);
        // array(linear_offset) |= mask;
    }

    void Mask(uint64_t linear_offset, const T mask)
    {
        if( linear_offset < max_linear_offset)
        {
            T* data=reinterpret_cast<T*>(reinterpret_cast<uint64_t>(array.Get_Data_Ptr()) + linear_offset);
            uint64_t page_mask = 1UL << (linear_offset>>12 & 0x3f);
            page_mask_array[linear_offset>>18] |= page_mask; dirty = true;
            *data |= mask;
        } 
        else 
        {
            std::cout<<"Linear offset out of range:"<<linear_offset<<std::endl;
            std::cout<<max_linear_offset<<std::endl;
            exit(0);
        }
    }

    void Unmask(uint64_t linear_offset)
    {
        if( linear_offset < max_linear_offset)
        {
            uint64_t page_mask = ~(1UL << (linear_offset>>12 & 0x3f));
            page_mask_array[linear_offset>>18] &= page_mask; dirty = true;
        }
        else 
        {
            std::cout<<"Linear offset out of range:"<<linear_offset<<std::endl;
            std::cout<<max_linear_offset<<std::endl;
            exit(0);
        }
    }

    bool MaskMatch(uint64_t linear_offset, const T mask)
    {
        if( linear_offset < max_linear_offset)
        {
            T* data=reinterpret_cast<T*>(reinterpret_cast<uint64_t>(array.Get_Data_Ptr()) + linear_offset);
            return (*data & mask);
        } 
        return false; // TODO: Do we need a FATAL ERROR?
    }

    bool Is_Set(uint64_t linear_offset,const T mask) const
    {
        if(linear_offset < max_linear_offset)
        {
            if(page_mask_array[linear_offset>>18] & (1UL<<(linear_offset>>12 & 0x3f)))
            {
                T* data=reinterpret_cast<T*>(reinterpret_cast<uint64_t>(array.Get_Data_Ptr()) + linear_offset);
                return (*data & mask);
            }
        }
        return false;        
    }
    
    bool Is_Set(std_array<int,dim> coord,const T mask) const
    {
        uint64_t linear_offset = T_MASK::Linear_Offset(coord);
        if(linear_offset < max_linear_offset)
        {
            if(page_mask_array[linear_offset>>18] & (1UL<<(linear_offset>>12 & 0x3f)))
            {
                T* data=reinterpret_cast<T*>(reinterpret_cast<uint64_t>(array.Get_Data_Ptr()) + linear_offset);
                return (*data & mask);
            }
        }
        return false;        
    }

    std::pair<const uint64_t*,unsigned> Get_Blocks() const
    {
        if(block_offsets.size())
            return std::pair<const uint64_t*,unsigned>(&block_offsets[0],block_offsets.size());
        else
            return std::pair<const uint64_t*,unsigned>((const uint64_t*)0,0);
    }

    void Refresh_Block_Offsets()
    {
        if(dirty)
            block_offsets = GenerateBlockOffsets();
        dirty = false;
    }

    std::vector<uint64_t> GenerateBlockOffsets()
    {
        std::vector<uint64_t> block_offsets;
        for (uint64_t i = 0; i<array_length; i++)
        {
            if(page_mask_array[i])
            {
                for (uint64_t pos=0; pos<64; pos++)
                {
                    if(page_mask_array[i] & (1UL<<pos))
                        block_offsets.push_back((i<<18)|(pos<<12));
                }
            }
        }
        return block_offsets;
    }

    void Clear_Bitmap()
    {
        memset(reinterpret_cast<void*>(page_mask_array),0,array_length*sizeof(uint64_t));
        dirty=true;
    }

    void Clear_Blocks()
    {
        std::vector<uint64_t>().swap(block_offsets);
        dirty=true;
    }

    void Clear()
    {Clear_Bitmap();Clear_Blocks();dirty=false;}

    void FillBitmapWithBlockOffsets(std::vector<uint64_t> block_offsets)
    {
        dirty = true;
        for (std::vector<uint64_t>::iterator it = block_offsets.begin(); it != block_offsets.end(); ++it)
        {
            uint64_t cur_offset = *it;
            page_mask_array[cur_offset>>18] |= ( 1UL << (cur_offset>>12 & 0x3f) );
        }
    }
};
}
#endif
