#ifndef __META_H_
#define __META_H_

#include <utility>
#include <tuple>
#include <array>
#include <vector>
#include <cuda_runtime.h>

namespace mn {

	template <class... ValueTypes>
	struct type_list {
		template <std::size_t N>
		using type = typename std::tuple_element<N, std::tuple<ValueTypes...>>::type;
	};
	template <typename Tuple>
	struct tuple_type_list {
		template <std::size_t N>
		using type = typename std::tuple_element<N, Tuple>::type;
	};

	template<unsigned int i, typename... ValueTypes>
	using AttribType = typename type_list<ValueTypes...>::template type<i>;
	template<unsigned int i, typename Tuple>
	using TupleAttribType = typename tuple_type_list<Tuple>::template type<i>;

	template<typename Func, typename... Args>
	void forArgs(Func&& f, Args&&... args) {
		return (void)std::initializer_list<int> {
			(f(std::forward<Args>(args)), 0)...
		};
	}
	template<typename Tuple, typename AllocFunc, size_t... I>
	void forIndexAlloc(AllocFunc&& f, std::index_sequence<I...>) {
		return (void)std::initializer_list<int> {
			(f(std::forward<size_t>(I), sizeof(TupleAttribType<I, Tuple>)), 0)...
		};
	}
	template<typename F, typename Tuple, size_t... I>
	decltype(auto) apply_impl(F&& f, Tuple t, std::index_sequence<I...>) {
		return std::forward<F>(f)(std::get<I>(std::forward<Tuple>(t))...);
	}
	template<typename F, typename Tuple>
	decltype(auto) apply(F&& f, Tuple t) {
		using Indices = std::make_index_sequence<std::tuple_size<std::decay_t<Tuple>>::value>;
		return apply_impl(std::forward<F>(f), std::forward<Tuple>(t), Indices{});
	}
	template<typename Func, typename Tpl>
	void forTuple(Func&& f, Tpl&& tpl) {
		apply(
			[&f](auto&&... xs) {forArgs(f, std::forward<decltype(xs)>(xs)...); },
			std::forward<Tpl>(tpl)
		);
	}

}

#endif
