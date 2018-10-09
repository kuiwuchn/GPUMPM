//#####################################################################
#ifndef __SPGrid_Mask_h__
#define __SPGrid_Mask_h__

#include <std_array.h>

namespace SPGrid{

template<int log2_struct,int D> class SPGrid_Mask_base;
template<int log2_struct,int log2_field,int D> class SPGrid_Mask;

//#####################################################################
// Class SPGrid_Mask_base (3D)
//#####################################################################
//#define BLOCKED_GRID_SIZE 1024

template<int log2_struct>
class SPGrid_Mask_base<log2_struct,3>
{
protected:

    enum : uint64_t {
        data_bits=log2_struct,            // Bits needed for indexing individual bytes within type T
        block_bits=12-data_bits,                   // Bits needed for indexing data elements within a block
        block_zbits=block_bits/3+(block_bits%3>0), // Bits needed for the z-coordinate of a data elements within a block
        block_ybits=block_bits/3+(block_bits%3>1), // Bits needed for the y-coordinate of a data elements within a block
        block_xbits=block_bits/3                   // Bits needed for the x-coordinate of a data elements within a block
    };

#ifndef BLOCKED_GRID_SIZE
    enum : uint64_t { // Bit masks for the upper 52 bits of memory addresses (page indices)
        page_zmask=(uint64_t) (0x9249249249249249UL<<(3-block_bits%3))&0xfffffffffffff000UL,
        page_ymask=(uint64_t) (0x2492492492492492UL<<(3-block_bits%3))&0xfffffffffffff000UL,
        page_xmask=(uint64_t) (0x4924924924924924UL<<(3-block_bits%3))&0xfffffffffffff000UL
    };
#else
    // WARNING - using this scheme will break all upsampling and downsampling!!
    // If blocked scheme we need to know the size of the total grid
    enum : uint64_t { // Bit masks for the upper 52 bits of memory addresses (page indices)
        page_zbits=NextLogTwo<BLOCKED_GRID_SIZE>::value - block_zbits,
        page_ybits=NextLogTwo<BLOCKED_GRID_SIZE>::value - block_ybits,
        page_xbits=NextLogTwo<BLOCKED_GRID_SIZE>::value - block_xbits,
        page_zmask= ((1UL<<page_zbits) - 1) << (12),
        page_ymask= ((1UL<<page_ybits) - 1) << (12+page_zbits),
        page_xmask= ((1UL<<page_xbits) - 1) << (12+page_zbits+page_ybits)
    };
#endif

public:
    enum {elements_per_block=1u<<block_bits};
};

//#####################################################################
// Class SPGrid_Mask (3D)
//#####################################################################

template<int log2_struct,int log2_field>
class SPGrid_Mask<log2_struct,log2_field,3>: public SPGrid_Mask_base<log2_struct,3>
{
public:
    enum {dim=3};
    enum {field_size=1<<log2_field};
    typedef SPGrid_Mask_base<log2_struct,dim> T_Mask_base;
    using T_Mask_base::data_bits;using T_Mask_base::block_bits;
    using T_Mask_base::block_xbits;using T_Mask_base::block_ybits;using T_Mask_base::block_zbits;
    using T_Mask_base::elements_per_block;

private:
    using T_Mask_base::page_xmask;using T_Mask_base::page_ymask;using T_Mask_base::page_zmask;

    enum : uint64_t { // Bit masks for the lower 12 bits of memory addresses (element indices within a page)
        element_zmask=(uint64_t) ((1<<block_zbits)-1)<<log2_field,
        element_ymask=(uint64_t) ((1<<block_ybits)-1)<<(log2_field+block_zbits),
        element_xmask=(uint64_t) ((1<<block_xbits)-1)<<(log2_field+block_zbits+block_ybits)
    };
    
public:
    enum  : uint64_t { 
        // Same as the corresponding element bit masks, but with the most significant bit the respective coordinate zeroed out
        element_z_lsbits=(element_zmask>>1)&element_zmask,
        element_y_lsbits=(element_ymask>>1)&element_ymask,
        element_x_lsbits=(element_xmask>>1)&element_xmask,
        
        // Same as the corresponding element bit masks, but with the least significant bit the respective coordinate zeroed out
        element_z_msbits=(element_zmask<<1)&element_zmask,
        element_y_msbits=(element_ymask<<1)&element_ymask,
        element_x_msbits=(element_xmask<<1)&element_xmask,

        // Just the most significant bit of the element bit mask for the respective coordinate
        element_z_msbit=element_zmask^element_z_lsbits,
        element_y_msbit=element_ymask^element_y_lsbits,
        element_x_msbit=element_xmask^element_x_lsbits,
        
        downsample_lower_mask = element_z_lsbits | element_y_lsbits | element_x_lsbits,
        upsample_lower_mask   = element_z_msbits | element_y_msbits | element_x_msbits,
        
        // "Left over bits" - lob=0 means page address starts with z-bit, lob=1 is x-bit, lob=2 is y-bit
        lob = data_bits%3,

        xloc = lob==0 ? 14 : ( lob==1 ? 12:13),
        yloc = lob==0 ? 13 : ( lob==1 ? 14:12),
        zloc = lob==0 ? 12 : ( lob==1 ? 13:14),

        u_zbit_shift = zloc - (log2_field+block_zbits-1),
        u_ybit_shift = yloc - (log2_field+block_zbits+block_ybits-1),
        u_xbit_shift = xloc - (log2_field+block_bits-1),

        bit12_mask = lob==0 ? element_z_msbit : ( lob==1 ? element_x_msbit : element_y_msbit ),
        bit13_mask = lob==0 ? element_y_msbit : ( lob==1 ? element_z_msbit : element_x_msbit ),
        bit14_mask = lob==0 ? element_x_msbit : ( lob==1 ? element_y_msbit : element_z_msbit )

    };
    
    enum  : uint64_t { // Bit masks for aggregate addresses
		zmask = (uint64_t)page_zmask | (uint64_t)element_zmask,
		ymask = (uint64_t)page_ymask | (uint64_t)element_ymask,
		xmask = (uint64_t)page_xmask | (uint64_t)element_xmask
    };
    enum  : uint64_t { 
        MXADD_Zmask=~zmask, 
        MXADD_Ymask=~ymask, 
        MXADD_Xmask=~xmask, 
        MXADD_Wmask=xmask|ymask|zmask
    };
    enum  : uint64_t {
        ODD_BITS=BitSpread<1,xmask>::value | BitSpread<1,ymask>::value | BitSpread<1,zmask>::value
    };


public:

    static unsigned int Bytes_Per_Element()
    {return 1u<<data_bits;}

    static unsigned int Elements_Per_Block()
    {return elements_per_block;}

    template<int i, int j, int k> struct LinearOffset
    {
      static const uint64_t value = BitSpread<i,xmask>::value | BitSpread<j,ymask>::value | BitSpread<k,zmask>::value;
    };

    inline static uint64_t Linear_Offset(const int i, const int j, const int k)
    {
#ifdef HASWELL
        return Bit_Spread(i,xmask)|Bit_Spread(j,ymask)|Bit_Spread(k,zmask);
#else
        return Bit_Spread<xmask>(i)|Bit_Spread<ymask>(j)|Bit_Spread<zmask>(k);
#endif
    }

    inline static uint64_t Linear_Offset(const std_array<int,3>& coord)
    {
#ifdef HASWELL
        return Bit_Spread(coord.data[0],xmask)|Bit_Spread(coord.data[1],ymask)|Bit_Spread(coord.data[2],zmask);
#else
        return Bit_Spread<xmask>(coord.data[0])|Bit_Spread<ymask>(coord.data[1])|Bit_Spread<zmask>(coord.data[2]);
#endif
    }
    inline static uint64_t Linear_Offset(const std_array<unsigned int,3>& coord)
    {
#ifdef HASWELL
        return Bit_Spread(coord.data[0],xmask)|Bit_Spread(coord.data[1],ymask)|Bit_Spread(coord.data[2],zmask);
#else
        return Bit_Spread<xmask>(coord.data[0])|Bit_Spread<ymask>(coord.data[1])|Bit_Spread<zmask>(coord.data[2]);
#endif
    }

    template<class T_mask_other>
    inline static uint64_t Translate_Linear_Offset(const uint64_t linear_offset_other)
    {
        if(T_mask_other::data_bits>data_bits){
            enum { page_spread_mask = LinearOffset<0x1fffffu<<T_mask_other::block_xbits,
                                                   0x1fffffu<<T_mask_other::block_ybits,
                                                   0x3fffffu<<T_mask_other::block_zbits>::value,
                   element_spread_mask = ((1<<log2_field)-1) |
                     ((1<<T_mask_other::block_zbits)-1)<<log2_field |
                     ((1<<T_mask_other::block_ybits)-1)<<(log2_field+block_zbits) |
                     ((1<<T_mask_other::block_xbits)-1)<<(log2_field+block_zbits+block_ybits) };
#ifdef HASWELL
            return Bit_Spread(linear_offset_other>>12,page_spread_mask) |
                Bit_Spread(linear_offset_other&0xfff,element_spread_mask);
#else
            return Bit_Spread<page_spread_mask>(linear_offset_other>>12) |
                   Bit_Spread<element_spread_mask>(linear_offset_other&0xfff);
#endif
        }
        else { int i,j,k; T_mask_other::LinearToCoord(linear_offset_other,&i,&j,&k); return Linear_Offset(i,j,k); }
    }

    inline static void LinearToCoord(uint64_t linear_offset, int* i, int* j, int* k)
    {
        *i = Bit_Pack(linear_offset,xmask);
        *j = Bit_Pack(linear_offset,ymask);
        *k = Bit_Pack(linear_offset,zmask);
    }
    inline static std_array<int,3> LinearToCoord(uint64_t linear_offset)
    {
        std_array<int,3> coord;
        coord.data[0] = Bit_Pack(linear_offset,xmask);
        coord.data[1] = Bit_Pack(linear_offset,ymask);
        coord.data[2] = Bit_Pack(linear_offset,zmask);
        return coord;
    }

    // 1-based offset calculation
    inline static uint64_t DownsampleOffset(uint64_t linear_offset)
    {
        return Packed_RightShift(Packed_Offset<1,1,1>(linear_offset));
    }


    inline static uint64_t Packed_RightShift(uint64_t linear_offset)
    {
        static uint64_t my_array[8] = {
            0,
            bit12_mask,
            bit13_mask,
            bit12_mask|bit13_mask,
            bit14_mask,
            bit12_mask|bit14_mask,
            bit13_mask|bit14_mask,
            bit12_mask|bit13_mask|bit14_mask};

        uint64_t upper = (linear_offset >> 3) & 0xfffffffffffff000UL;
        uint64_t lower = (linear_offset >> 1) & downsample_lower_mask;
        uint64_t result = upper | lower | my_array[(linear_offset>>12) & 0x7UL];
        return result;
    }
    
    // 1-based offset calculation
    inline static uint64_t UpsampleOffset(uint64_t linear_offset)
    {
        return Packed_Offset<-1,-1,-1>(Packed_LeftShift(linear_offset));
    }

    inline static uint64_t Packed_LeftShift(uint64_t linear_offset)
    {
        uint64_t upper = (linear_offset << 3) & 0xfffffffffffff000UL;
        uint64_t lower = (linear_offset << 1) & upsample_lower_mask;
					  
        uint64_t x_bit_place = (linear_offset & element_x_msbit)<<u_xbit_shift;
        uint64_t y_bit_place = (linear_offset & element_y_msbit)<<u_ybit_shift;
        uint64_t z_bit_place = (linear_offset & element_z_msbit)<<u_zbit_shift;
					  
        uint64_t result = upper | lower | x_bit_place | y_bit_place | z_bit_place;
        return result;
    }

    inline static uint64_t Packed_Add(const uint64_t i,const uint64_t j)
    {
        uint64_t x_result=( (i | MXADD_Xmask) + (j & ~MXADD_Xmask) ) & ~MXADD_Xmask;
        uint64_t y_result=( (i | MXADD_Ymask) + (j & ~MXADD_Ymask) ) & ~MXADD_Ymask;
        uint64_t z_result=( (i | MXADD_Zmask) + (j & ~MXADD_Zmask) ) & ~MXADD_Zmask;
        uint64_t w_result=( (i | MXADD_Wmask) + (j & ~MXADD_Wmask) ) & ~MXADD_Wmask;
        uint64_t result=x_result | y_result | z_result | w_result;
        return result;
    }

    template<uint64_t j>
    inline static uint64_t Packed_Add_Imm(const uint64_t i)
    {
        uint64_t x_result=( (i | MXADD_Xmask) + (j & ~MXADD_Xmask) ) & ~MXADD_Xmask;
        uint64_t y_result=( (i | MXADD_Ymask) + (j & ~MXADD_Ymask) ) & ~MXADD_Ymask;
        uint64_t z_result=( (i | MXADD_Zmask) + (j & ~MXADD_Zmask) ) & ~MXADD_Zmask;
        uint64_t w_result=( (i | MXADD_Wmask) + (j & ~MXADD_Wmask) ) & ~MXADD_Wmask;
        uint64_t result=x_result | y_result | z_result | w_result;
        return result;
    }
    
    template<int di, int dj, int dk>
    inline static uint64_t Packed_Offset(const uint64_t pI)
    {
        static const uint64_t pJ = LinearOffset<di,dj,dk>::value;

        uint64_t x_result=( (pI | MXADD_Xmask) + (pJ & ~MXADD_Xmask) ) & ~MXADD_Xmask;
        uint64_t y_result=( (pI | MXADD_Ymask) + (pJ & ~MXADD_Ymask) ) & ~MXADD_Ymask;
        uint64_t z_result=( (pI | MXADD_Zmask) + (pJ & ~MXADD_Zmask) ) & ~MXADD_Zmask;
        
        uint64_t result=x_result | y_result | z_result;
        return result;
    }

    template<int di>
    inline static uint64_t Packed_OffsetXdim(const uint64_t pI)
    {
        static const uint64_t pJ = LinearOffset<di,0,0>::value;

		uint64_t x_result = ((pI | MXADD_Xmask) + (pJ & ~MXADD_Xmask)) & ~MXADD_Xmask;
		uint64_t not_x_result = pI & MXADD_Xmask;
        
        uint64_t result=x_result | not_x_result;
        return result;
    }

    template<int dj>
    inline static uint64_t Packed_OffsetYdim(const uint64_t pI)
    {
        static const uint64_t pJ = LinearOffset<0,dj,0>::value;

        uint64_t y_result=( (pI | MXADD_Ymask) + (pJ & ~MXADD_Ymask) ) & ~MXADD_Ymask;
        uint64_t not_y_result= pI & MXADD_Ymask;
        			
        uint64_t result=y_result | not_y_result;
        return result;
    }
    
    template<int dk>
    inline static uint64_t Packed_OffsetZdim(const uint64_t pI)
    {
        static const uint64_t pJ = LinearOffset<0,0,dk>::value;

		uint64_t z_result = ((pI | MXADD_Zmask) + (pJ & ~MXADD_Zmask)) & ~MXADD_Zmask;
		uint64_t not_z_result = pI & MXADD_Zmask;

		uint64_t result = z_result | not_z_result;
        return result;
    }

};
}
#endif
