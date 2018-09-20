#ifndef __SPGRID_ARRAY_H_
#define __SPGRID_ARRAY_H_

#include <limits>
#include <algorithm>
#include "SPGridAllocator.h"
#include "SPGridMask.h"
#include "SPGridBlock.h"

namespace SPGrid {

	// inline constexpr bool g_bCheckBound = true;
//#define SPGRID_CHECK_BOUNDS

	template<typename T, typename MASK>
	class SPGridArray {
	public:
		using Data = T;
		using DataType = std::remove_const_t<T>;
		using Mask = MASK;
		enum { Dim = MASK::Dim };

		SPGridArray(void* const ptr, const SPGridBlock<Dim>& domain) : _gridPtr(ptr), _domain(domain) {}

		inline T& operator()(const std::array<int, 3>& coord) {
			//static_assert(Dim == 3);
#ifdef SPGRID_CHECK_BOUNDS
			_domain.checkBounds(coord[0], coord[1], coord[2]);
#endif
			return *reinterpret_cast<T*>(reinterpret_cast<uint64_t>(_gridPtr) + Mask::linear_offset(coord));
		}
		inline T& operator()(const int i, const int j, const int k) {
			//static_assert(Dim == 3);
#ifdef SPGRID_CHECK_BOUNDS
			_domain.checkBounds(i, j, k);
#endif
			return *reinterpret_cast<T*>(reinterpret_cast<uint64_t>(_gridPtr) + Mask::linear_offset(i, j, k));
		}
		inline T& operator()(const uint64_t offset) {
			return *reinterpret_cast<T*>(reinterpret_cast<uint64_t>(_gridPtr) + offset);
		}
		inline DataType get(const std::array<int, 3>& coord) const {
			//static_assert(Dim == 3);
#ifdef SPGRID_CHECK_BOUNDS
			_domain.checkBounds(coord[0], coord[1], coord[2]);
#endif
			return *reinterpret_cast<DataType*>(reinterpret_cast<uint64_t>(_gridPtr) + Mask::linear_offset(coord));
		}
		inline DataType get(const int i, const int j, const int k) const {
			//static_assert(Dim == 3);
#ifdef SPGRID_CHECK_BOUNDS
			_domain.checkBounds(i, j, k);
#endif
			return *reinterpret_cast<DataType*>(reinterpret_cast<uint64_t>(_gridPtr) + Mask::linear_offset(i, j, k));
		}
		inline DataType get(const uint64_t offset) const {
			return *reinterpret_cast<DataType*>(reinterpret_cast<uint64_t>(_gridPtr) + offset);
		}

		inline bool probeValue(const std::array<int, 3>& coord, DataType& val) const {
			val = *reinterpret_cast<DataType*>(reinterpret_cast<uint64_t>(_gridPtr) + Mask::linear_offset(coord[0], coord[1], coord[2]));
			return std::abs(val) > std::numeric_limits<T>::epsilon();
		}
		inline bool probeValue(const int i, const int j, const int k, DataType& val) const {
			val = *reinterpret_cast<DataType*>(reinterpret_cast<uint64_t>(_gridPtr) + Mask::linear_offset(i, j, k));
			return std::abs(val) > std::numeric_limits<T>::epsilon();
		}
		inline bool probeValue(const uint64_t offset, DataType& val) const {
			val = *reinterpret_cast<DataType*>(reinterpret_cast<uint64_t>(_gridPtr) + offset);
			return std::abs(val) > std::numeric_limits<T>::epsilon();
		}

		template<int di, int dj, int dk>
		inline T& operator()(const uint64_t offset) {
			//static_assert(Dim == 3);
#ifdef SPGRID_CHECK_BOUNDS
			_domain.checkBounds(i, j, k);
#endif
			auto ptr = reinterpret_cast<T*>(reinterpret_cast<uint64_t>(_gridPtr) + Mask::Packed_Offset<di, dj, dk>(offset));
			return *ptr;
		}
		template<int di, int dj, int dk>
		inline T& get(const uint64_t offset) {
			//static_assert(Dim == 3);
#ifdef SPGRID_CHECK_BOUNDS
			_domain.checkBounds(i, j, k);
#endif
			auto ptr = reinterpret_cast<T*>(reinterpret_cast<uint64_t>(_gridPtr) + Mask::Packed_Offset<di, dj, dk>(offset));
			return *ptr;
		}

		inline void* getPtr() noexcept { return _gridPtr; }
		inline const auto& getDomainRef() const { return _domain; }
		inline auto& refDomain() { return _domain; }

	private:
		const SPGridBlock<Dim>& _domain;
		void* const _gridPtr;
	};

}

#endif