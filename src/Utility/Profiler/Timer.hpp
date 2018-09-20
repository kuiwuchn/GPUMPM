/**	\file	Timer.hpp
 *	\brief	CPU timer based on std::chrono and std::basic_ostream
 *	\ref	https://www.daniweb.com/programming/software-development/code/445750/simple-timer-using-c-11-s-chrono-library
 */

#ifndef __TIMER_HPP_
#define __TIMER_HPP_

#include <fstream>
#include <ostream>
#include <vector>
#include <chrono>

namespace mn {

	class Timer {
		typedef std::chrono::high_resolution_clock	HRC;
		typedef std::chrono::milliseconds			MS;	///< default timer unit
		typedef std::chrono::microseconds			US;	///< default timer unit

		HRC::time_point	_start, _end;
		//std::vector<std::pair<std::string, MS>>	_records;
		std::vector<std::pair<std::string, US>>	_records;

		//MS		elapsed() const { return std::chrono::duration_cast<MS>(_end - _start); }
		US		elapsed() const { return std::chrono::duration_cast<US>(_end - _start); }
	public:
		Timer() {}
		~Timer() {}

		void	tick() { _start = HRC::now(); }
		void	tock() { _end = HRC::now(); }
		void	clear() { _records.clear(); }

		void	record(std::string tag) {
			tock();
			_records.emplace_back(tag, elapsed());
		}

		void	blankLine() { _records.emplace_back("", US()); }
		void	log(std::ofstream& fs, bool bclear = true) {
			for (auto& pair : _records)
				if (!pair.first.empty())
					//fs << pair.first.c_str() << ": " << pair.second.count() << std::endl;
					fs << pair.first.c_str() << ": " << pair.second.count() / 1000.f << std::endl;
				else
					fs << '\n';
			if (bclear)
				clear();
		}

		template <typename T, typename Traits>
		friend std::basic_ostream<T, Traits>& operator<<(std::basic_ostream<T, Traits>& out, const Timer& timer) {
			return out << timer.elapsed().count() / 1000.f;
		}
	};

}

#endif
