#ifndef __ALLOC_META_H_
#define __ALLOC_META_H_

#include "Meta.h"

namespace mn {

	template<typename... Args>
	auto make_vector(Args&&... args) {
		using Item = std::common_type_t<Args...>;
		std::vector<Item>	result(sizeof...(Args));
		// works as a building block
		forArgs(
			[&result](auto&& x) {result.emplace_back(std::forward<decltype(x)>(x)); },
			std::forward<Args>(args)...
		);
		return result;
	}

	template<typename Tuple>
	auto cuda_allocs(unsigned int size) {
		std::array<void*, std::tuple_size<std::decay_t<Tuple>>::value> addrs;
		forIndexAlloc<Tuple>(
			[&addrs, size](size_t&& i, size_t&& typeSize) { cudaMalloc(&addrs[i], typeSize * size); },
			//[&addrs, size](size_t&& i, size_t&& typeSize) { printf("alloc %llu-th attrib of %llu * %u\n", i, typeSize, size); },
			std::make_index_sequence<std::tuple_size<std::decay_t<Tuple>>::value>{}
		);
		return addrs;
	}
	
	template<typename Tuple>
	auto cuda_allocs(const unsigned int* sizes) {
		std::array<void*, std::tuple_size<std::decay_t<Tuple>>::value> addrs;
		forIndexAlloc<Tuple>(
			[&addrs, &sizes](size_t&& i, size_t&& typeSize) { cudaMalloc(&addrs[i], typeSize * sizes[i]); },
			std::make_index_sequence<std::tuple_size<std::decay_t<Tuple>>::value>{}
		);
		return addrs;
	}

	template<int num>
	void cuda_free(std::array<void*, num> &_attribs) {
		for (auto& attrib : _attribs)
			cudaFree(attrib);
	}
}

#endif