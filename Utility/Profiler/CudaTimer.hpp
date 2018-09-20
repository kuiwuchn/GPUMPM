/**	\file	CudaTimer.hpp
 *	\brief	CUDA GPU timer based on cudaEvent
 */

#ifndef __CUDA_TIMER_HPP_
#define __CUDA_TIMER_HPP_

#include <fstream>
#include <ostream>
#include <vector>

#include <cuda_runtime.h>

namespace mn {

	class CudaTimer {
		cudaEvent_t	_start, _end;
		std::vector<std::pair<std::string, float>>	_records;
		int lastSectionLen{ 0 };
	public:
		CudaTimer() {
			cudaEventCreate(&_start);
			cudaEventCreate(&_end);
		}
		~CudaTimer() {
			cudaEventDestroy(_start);
			cudaEventDestroy(_end);
		}

		void	tick(cudaStream_t streamid = 0) { cudaEventRecord(_start, streamid); }
		void	tock(cudaStream_t streamid = 0) { cudaEventRecord(_end, streamid); }
		float	elapsed() const {
			float ms;
			cudaEventSynchronize(_end);
			cudaEventElapsedTime(&ms, _start, _end);
			return ms;
		}
		void	clear() { _records.clear(); }
		void	record(std::string tag) {
			tock();
			_records.emplace_back(tag, elapsed());
		}

		void	recordSection(std::string msg) {
			float sec = 0;
			auto it = _records.end();
			for (--it; it != _records.begin() && !it->first.empty() && it->first[0] != '@' && it->first[0] != '#'; --it)
				if (it->first[0] != '#' && it->first[0] != '$')
					sec += it->second;
			if (it == _records.begin())
				sec += it->second;
			_records.emplace_back("@" + std::string{ msg }, sec);
			std::cout << "@" + std::string{ msg } + ": " << sec << "\n\n";
			_records.emplace_back("", 0);
		}

		void	blankLine() {
			float sum = 0;
			auto it = _records.end();
			for (--it; it != _records.begin() && it->first != "###"; --it) 
				if (it->first[0] != '@' && it->first[0] != '#' && it->first[0] != '$')
					sum += it->second;
			if (it == _records.begin())
				sum += it->second;
			_records.emplace_back("#timestep_time", sum);
			std::cout << "#timestep_time: " << sum << "\n\n";
			_records.emplace_back("###", 0);
		}

		void	recordFrame() {
			float sum = 0;
			auto it = _records.end();
			for (--it; it != _records.begin() && it->first != "$$$"; --it) 
				if (it->first == "#timestep_time")
					sum += it->second;
			_records.emplace_back("$frame_time", sum);
			std::cout << "$frame_time: " << sum << "\n\n";
			_records.emplace_back("$$$", 0);
		}

		void	log(std::ofstream& fs, bool bclear = true) {
			for (auto& pair : _records)
				if (!pair.first.empty())
					fs << pair.first.c_str() << ": " << pair.second << '\n';
				else
					fs << '\n';
			if (bclear)
				clear();
		}

		template <typename T, typename Traits>
		friend std::basic_ostream<T, Traits>& operator<<(std::basic_ostream<T, Traits>& out, const CudaTimer& timer) {
			return out << timer.elapsed();
		}
	};

}

#endif