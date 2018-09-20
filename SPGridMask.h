#ifndef __SPGRID_MASK_H_
#define __SPGRID_MASK_H_

#include <array>
#include "SPGridUtility.h"

namespace SPGrid {

	template<typename Channels, int Dim> struct SPGridMask;
	template<typename Channels, typename Field, int Dim> struct SPGridFieldMask;

	template<typename Channels> 
	struct SPGridMask<Channels, 3> {
		enum { Dim = 3 };
		enum Block {
			channelSize = sizeof(Channels),
			channelBitCount = bitCount(channelSize),

			blockBitCount = 12 - channelBitCount,
			blockZBitCount = blockBitCount / 3 + (blockBitCount % 3 > 0),
			blockYBitCount = blockBitCount / 3 + (blockBitCount % 3 > 1),
			blockXBitCount = blockBitCount / 3,

			blockSize = 1u << blockBitCount
		};

		/// upper 52 bits of memory addresses (page indices)
		enum BlockMask : uint64_t {
			blockZMask = (0x9249249249249249ULL << (3 - blockBitCount % 3)) & 0xfffffffffffff000ULL,
			blockYMask = (0x2492492492492492ULL << (3 - blockBitCount % 3)) & 0xfffffffffffff000ULL,
			blockXMask = (0x4924924924924924ULL << (3 - blockBitCount % 3)) & 0xfffffffffffff000ULL
		};
	};

	template<typename Channels, typename Field>
	struct SPGridFieldMask<Channels, Field, 3> : SPGridMask<Channels, 3> {
		using GridMask = SPGridMask<Channels, 3>;
		enum {
			fieldSize = sizeof(Field),
			fieldBitCount = bitCount(fieldSize)
		};

		/// lower 12 bits [X..Y..Z..Field..]
		enum CellMask : uint64_t {
			cellZMask = ((1ULL << GridMask::blockZBitCount) - 1) << fieldBitCount,
			cellYMask = ((1ULL << GridMask::blockYBitCount) - 1) << (fieldBitCount + GridMask::blockZBitCount),
			cellXMask = ((1ULL << GridMask::blockXBitCount) - 1) << (fieldBitCount + GridMask::blockZBitCount + GridMask::blockYBitCount)
		};

		enum OffsetMask : uint64_t {
			/// Same as the corresponding element bit masks, but with the most significant bit the respective coordinate zeroed out
			cellZLSBits = (cellZMask >> 1) & cellZMask,
			cellYLSBits = (cellYMask >> 1) & cellYMask,
			cellXLSBits = (cellXMask >> 1) & cellXMask,

			/// Same as the corresponding element bit masks, but with the least significant bit the respective coordinate zeroed out
			cellZMSBits = (cellZMask << 1) & cellZMask,
			cellYMSBits = (cellYMask << 1) & cellYMask,
			cellXMSBits = (cellXMask << 1) & cellXMask,

			/// Just the most significant bit of the element bit mask for the respective coordinate
			cellZMSBit = cellZMask ^ cellZLSBits,
			cellYMSBit = cellYMask ^ cellYLSBits,
			cellXMSBit = cellXMask ^ cellXLSBits,

			downsampleLowerMask = cellZLSBits | cellYLSBits  | cellXLSBits,
			upsampleLowerMask = cellZMSBits | cellYMSBits  | cellXMSBits,

			/// "Left over bits" - lob=0 means page address starts with z-bit, lob=1 is x-bit, lob=2 is y-bit
			lob = (66 - GridMask::blockBitCount) % 3,

			xloc = lob == 0 ? 14 : (lob == 1 ? 12 : 13),
			yloc = lob == 0 ? 13 : (lob == 1 ? 14 : 12),
			zloc = lob == 0 ? 12 : (lob == 1 ? 13 : 14),

			uZBitShift = zloc - (fieldBitCount + GridMask::blockZBitCount - 1),
			uYBitShift = yloc - (fieldBitCount + GridMask::blockZBitCount + GridMask::blockYBitCount - 1),
			uXBitShift = xloc - (fieldBitCount + GridMask::blockBitCount - 1),

			bit12Mask = lob == 0 ? cellZMSBit : (lob == 1 ? cellXMSBit : cellYMSBit),
			bit13Mask = lob == 0 ? cellYMSBit : (lob == 1 ? cellZMSBit : cellXMSBit),
			bit14Mask = lob == 0 ? cellXMSBit : (lob == 1 ? cellYMSBit : cellZMSBit)
		};

		/// final mask for each dimension
		enum Mask : uint64_t { 
			ZMask = GridMask::blockZMask | static_cast<uint64_t>(cellZMask),
			YMask = GridMask::blockYMask | static_cast<uint64_t>(cellYMask),
			XMask = GridMask::blockXMask | static_cast<uint64_t>(cellXMask)
		};
		enum Compute : uint64_t {
			MXADDZmask = ~ZMask,
			MXADDYmask = ~YMask,
			MXADDXmask = ~XMask,
			MXADDWmask = XMask | YMask | ZMask
		};
		//enum {
		//	OddBits = BitSpread<1, xmask>::value | BitSpread<1, ymask>::value | BitSpread<1, zmask>::value
		//};

		static inline uint32_t bytesPerElement() noexcept { return 1u << GridMask::channelBitCount; }
		static inline uint32_t elementsPerBlock() noexcept { return GridMask::blockSize; }

		static constexpr uint64_t Linear_Offset(const int i, const int j, const int k) noexcept {
			return Bit_Spread(i, XMask) | Bit_Spread(j, YMask) | Bit_Spread(k, ZMask);
		}
		static inline uint64_t linear_offset(const int i, const int j, const int k) noexcept {
			return bit_spread(i, XMask) | bit_spread(j, YMask) | bit_spread(k, ZMask);
		}
		static inline uint64_t linear_offset(const std::array<int, 3>& coord) noexcept {
			return bit_spread(coord[0], XMask) | bit_spread(coord[1], YMask) | bit_spread(coord[2], ZMask);
		}

		static constexpr std::array<int, 3> Offset_To_Coord(const uint64_t offset) {
			std::array<int, 3>	coord;
			coord[0] = Bit_Pack(offset, XMask);
			coord[1] = Bit_Pack(offset, YMask);
			coord[2] = Bit_Pack(offset, ZMask);
			return coord;
		}
		static constexpr std::array<int, 3> offset_to_coord(const uint64_t offset) {
			std::array<int, 3>	coord;
			coord[0] = bit_pack(offset, XMask);
			coord[1] = bit_pack(offset, YMask);
			coord[2] = bit_pack(offset, ZMask);
			return coord;
		}
		static inline void offset_to_coord(const uint64_t offset, int* i, int* j, int* k) {
			*i = bit_pack(offset, XMask);
			*j = bit_pack(offset, YMask);
			*k = bit_pack(offset, ZMask);
		}
		static inline uint64_t packed_right_shift(uint64_t linearOffset) noexcept {
			static uint64_t my_array[8] = {
				0,
				bit12Mask,
				bit13Mask,
				bit12Mask | bit13Mask,
				bit14Mask,
				bit12Mask | bit14Mask,
				bit13Mask | bit14Mask,
				bit12Mask | bit13Mask | bit14Mask };
			uint64_t upper = (linearOffset >> 3) & 0xfffffffffffff000UL;
			uint64_t lower = (linearOffset >> 1) & downsampleLowerMask;
			uint64_t result = upper | lower | my_array[(linearOffset >> 12) & 0x7UL];
			return result;
		}
		static constexpr uint64_t packed_left_shift(uint64_t linearOffset) noexcept {
			uint64_t upper = (linearOffset << 2) & 0xfffffffffffff000ULL;
			uint64_t lower = (linearOffset << 1) & upsampleLowerMask;

			uint64_t xBitPlace = (linearOffset & cellXMSBit) << uXBitShift;
			uint64_t yBitPlace = (linearOffset & cellYMSBit) << uYBitShift;

			uint64_t result = upper | lower | xBitPlace | yBitPlace;
			return result;
		}
		/// fill the gap with 1s (by '|' the negated mask) between bits of the same dimension
		static constexpr uint64_t packed_add(const uint64_t i, const uint64_t j) noexcept {
			uint64_t x_result = ((i | MXADDXmask) + (j & XMask)) & XMask;
			uint64_t y_result = ((i | MXADDYmask) + (j & YMask)) & YMask;
			uint64_t z_result = ((i | MXADDZmask) + (j & ZMask)) & ZMask;
			uint64_t w_result = ((i | MXADDWmask) + (j & ~MXADDWmask)) & ~MXADDWmask;
			uint64_t result = x_result | y_result | z_result | w_result;
			return result;
		}
		template<int di, int dj, int dk>
		static uint64_t Packed_Offset(const uint64_t pI) noexcept {
			constexpr uint64_t pJ = Linear_Offset(di, dj, dk);
			uint64_t x_result = ((pI | MXADDXmask) + (pJ & XMask)) & XMask;
			uint64_t y_result = ((pI | MXADDYmask) + (pJ & YMask)) & YMask;
			uint64_t z_result = ((pI | MXADDZmask) + (pJ & ZMask)) & ZMask;
			return x_result | y_result | z_result;
		}
		template<int di>
		static constexpr uint64_t Packed_Offset_X(const uint64_t pI) noexcept {
			constexpr uint64_t pJ = Linear_Offset(di, 0, 0);
			uint64_t x_result = ((pI | MXADDXmask) + (pJ & XMask)) & XMask;
			uint64_t not_x_result = pI & MXADDXmask;
			return x_result | not_x_result;
		}
		template<int dj>
		static constexpr uint64_t Packed_Offset_Y(const uint64_t pI) noexcept {
			constexpr uint64_t pJ = Linear_Offset(0, dj, 0);
			uint64_t y_result = ((pI | MXADDYmask) + (pJ & YMask)) & YMask;
			uint64_t not_y_result = pI & MXADDYmask;
			return y_result | not_y_result;
		}
		template<int dk>
		static constexpr uint64_t Packed_Offset_Z(const uint64_t pI) noexcept {
			constexpr uint64_t pJ = Linear_Offset(0, 0, dk);
			uint64_t z_result = ((pI | MXADDZmask) + (pJ & ZMask)) & ZMask;
			uint64_t not_z_result = pI & MXADDZmask;
			return z_result | not_z_result;
		}
		// 1-based offset calculation
		static constexpr uint64_t down_sample_offset(const uint64_t linearOffset) noexcept {
			return packed_right_shift(Packed_Offset<1, 1, 1>(linearOffset));
		}
		// 1-based offset calculation
		static constexpr uint64_t up_sample_offset(const uint64_t linearOffset) noexcept {
			return Packed_Offset<-1, -1, -1>(packed_left_shift(linearOffset));
		}

	};

}

#endif