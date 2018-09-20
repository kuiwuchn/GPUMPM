#ifndef __SPGRID_UTILITY_H_
#define __SPGRID_UTILITY_H_

#include <stdint.h>
#include <algorithm>

namespace SPGrid {

	//template<uint32_t N> constexpr uint32_t bitLength() { return bitLength < (N >> 1) >() + 1; }
	//template<> constexpr uint32_t bitLength<0>() { return 0; }
	/**
	 *	\fn uint32_t bitLength(uint32_t N)
	 *	\brief compute the count of significant digits of a number
	 *	\param N the number
	 */
	constexpr uint32_t bitLength(uint32_t N) noexcept {
		if (N) return bitLength(N >> 1) + 1;
		else return 0;
	}
	/**
	 *	\fn uint32_t bitCount(uint32_t N)
	 *	\brief compute the count of digits required to express integers in [0, N)
	 *	\param N the maximum of the range
	 */
	constexpr uint32_t bitCount(uint32_t N) noexcept { return bitLength(N - 1); }
	/**
	 *	\fn uint32_t nextPowerOfTwo(uint32_t i)
	 *	\brief compute the next power of two bigger than the number i
	 *	\param i the number
	 */
	constexpr uint32_t nextPowerOfTwo(uint32_t i) noexcept {
		i--; i |= i >> 1; i |= i >> 2; i |= i >> 4; i |= i >> 8; i |= i >> 16; return i + 1;
	}
	/**
	 *	\fn size_t OffsetOfMember(T_FIELD T::*field)
	 *	\brief compute the offset of a data member in its class
	 *	\param field the data member of the class
	 */
	template<typename T, typename T_FIELD> 
	size_t OffsetOfMember(T_FIELD T::*field) {
		return (size_t)((char*)&(((T*)0)->*field) - (char*)0);
	}
	/**
	 *	\fn T maximum(const T& x1, const T& x2)
	 */
	template<class T> constexpr T maximum(const T& x1, const T& x2) {
		return std::max(x1, x2);
	}
	/**
	 *	\fn T maximum(const T& x1, const T& x2, cosnt T& x3)
	 */
	template<class T> constexpr T maximum(const T& x1, const T& x2, const T& x3) {
		return std::max(x1, std::max(x2, x3));
	}
	/**
	 *	\fn T maximum(const T& x1, const T& x2, cosnt T& x3, cosnt T& x4)
	 */
	template<class T> constexpr T maximum(const T& x1, const T& x2, const T& x3, const T& x4) {
		return std::max(std::max(x1, x2), std::max(x3, x4));
	}
	template<class T> constexpr T minimum(const T& x1, const T& x2) {
		return std::min(x1, x2);
	}
	template<class T> constexpr T minimum(const T& x1, const T& x2, const T& x3) {
		return std::min(x1, std::min(x2, x3));
	}
	template<class T> constexpr T minimum(const T& x1, const T& x2, const T& x3, const T& x4) {
		return std::min(std::min(x1, x2), std::min(x3, x4));
	}

	constexpr uint32_t reverse_bits(uint32_t n) noexcept {
		n = (n >> 1) & 0x55555555 | (n << 1) & 0xaaaaaaaa;
		n = (n >> 2) & 0x33333333 | (n << 2) & 0xcccccccc;
		n = (n >> 4) & 0x0f0f0f0f | (n << 4) & 0xf0f0f0f0;
		n = (n >> 8) & 0x00ff00ff | (n << 8) & 0xff00ff00;
		n = (n >> 16) & 0x0000ffff | (n << 16) & 0xffff0000;
		return n;
	}
	constexpr uint64_t reverse_bits(uint64_t n) {
		uint64_t result = reverse_bits(static_cast<uint32_t>(n & 0xffffffff));
		return result << 32 | reverse_bits(static_cast<uint32_t>(n >> 32));
	}
	constexpr uint64_t clz64(uint64_t num) noexcept {
		unsigned char cnt = 0;
		uint64_t offset = 0x8000000000000000;
		while (!(num & offset) && offset) { cnt++; offset >>= 1; }
		return cnt;
	}
	constexpr uint64_t Bit_Spread(uint32_t bits, uint64_t mask) noexcept {
		// gcc: __builtin_ctz
		// msvd: __tzcnt64
		uint64_t rmask = reverse_bits(mask);	//_byteswap_uint64 doesn't apply
		uint64_t result = 0;
		unsigned char lz = 0;
		while (rmask) {
			lz = clz64(rmask) + 1;
			result = result << lz | (bits & 1);
			bits >>= 1, rmask <<= lz;
		}
		result = reverse_bits(result) >> clz64(mask);
		return result;
	}
	constexpr uint32_t Bit_Pack(uint64_t bits, uint64_t mask) noexcept {
		uint64_t ulresult = 0;
		uint64_t uldata = bits;
		int count = 0;
		ulresult = 0;

		uint64_t rmask = reverse_bits(mask);
		unsigned char lz = 0;

		while (rmask) {
			lz = clz64(rmask);
			uldata >>= lz;
			ulresult = ulresult << 1 | (uldata & 1);
			uldata >>= 1;
			rmask <<= lz + 1;
			count++;
		}
		ulresult <<= 64 - count; // 64 means 64 bits ... maybe not use a constant 64 ...
		ulresult = reverse_bits(ulresult);
		return static_cast<uint32_t>(ulresult);
	}
	/*
	constexpr uint32_t Index_To_Morton_2D_32(uint32_t x, uint32_t y) noexcept {
		;
	}
	constexpr void Morton_To_Index_2D_32(uint64_t m, uint32_t& x, uint32_t& y) noexcept {
		;
	}
	*/

	uint64_t bit_spread(uint32_t bits, uint64_t mask) noexcept;
	uint32_t bit_pack(uint64_t bits, uint64_t mask) noexcept;
	uint64_t index_to_morton_2d_64(uint32_t x, uint32_t y) noexcept;
	void morton_to_index_2d_64(uint64_t m, uint32_t& x, uint32_t& y) noexcept;
	uint32_t index_to_morton_2d_32(uint32_t x, uint32_t y) noexcept;
	void morton_to_index_2d_32(uint32_t m, uint32_t& x, uint32_t& y) noexcept;

	void checkTotalMemory();
	void checkProcessMemory();
	bool checkAddressResident(void * const addr);
}

#endif