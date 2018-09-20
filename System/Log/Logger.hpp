#ifndef __LOGGER_HPP_
#define __LOGGER_HPP_

#include <string>
#include <iostream>
#include <MnBase/Singleton.h>
#include <Utility/Profiler/CudaTimer.hpp>
#include <Utility/Profiler/Timer.hpp>

namespace mn {

	enum class TimerType { CPU, GPU };

	class Logger : public ManagedSingleton<Logger> {
	public:
		Logger()	{}
		~Logger()	{}

		template<TimerType>
		static void tick();
		template<TimerType>
		static void tock(std::string message);
		template<TimerType>
		static void blankLine();
		template<TimerType>
		static void recordSection(std::string msg);
		template<TimerType>
		static void recordFrame();

		static void message(std::string filename);

		//static void record(std::string filename);

	private:
		CudaTimer	_kGPUTimer;
		Timer		_kCPUTimer;
		std::vector<std::string>	_kInfos;
	};

	template <>
	inline void Logger::tick<TimerType::GPU>() { getInstance()->_kGPUTimer.tick(); }
	template <>
	inline void Logger::tick<TimerType::CPU>() { getInstance()->_kCPUTimer.tick(); }

	template <>
	inline void Logger::tock<TimerType::GPU>(std::string message) { 
		getInstance()->_kGPUTimer.record(message);
		std::cout << message << ": " << getInstance()->_kGPUTimer << std::endl; 
	}
	template <>
	inline void Logger::tock<TimerType::CPU>(std::string message) { 
		getInstance()->_kCPUTimer.record(message);
		std::cout << message << ": " << getInstance()->_kCPUTimer << std::endl; 
	}
	template<>
	inline void Logger::recordSection<TimerType::GPU>(std::string msg) { getInstance()->_kGPUTimer.recordSection(msg); }
	template<>
	inline void Logger::recordFrame<TimerType::GPU>() { getInstance()->_kGPUTimer.recordFrame(); }
	template<>
	inline void Logger::blankLine<TimerType::GPU>() { getInstance()->_kGPUTimer.blankLine(); }
	template<>
	inline void Logger::blankLine<TimerType::CPU>() { getInstance()->_kCPUTimer.blankLine(); }

	/*
	inline void Logger::record(std::string filename) {
		using namespace std::experimental::filesystem;
		path outputTarget(filename.data());
		if (outputTarget.empty()) return;
		if (!exists(outputTarget.parent_path()))
			create_directory(outputTarget.parent_path());

		std::ofstream ofs;
		while (exists(outputTarget)) {
			outputTarget = outputTarget.parent_path().string() + "\\" + outputTarget.stem().string() + "_" + outputTarget.extension().string();
		}
		ofs.open(outputTarget.string());
		if (ofs.is_open()) {
			getInstance()->_kCPUTimer.log(ofs, true);
			ofs << '\n';

			for (auto& str : getInstance()->_kInfos)
				ofs << str << '\n';
			ofs << '\n';
			getInstance()->_kInfos.clear();

			getInstance()->_kGPUTimer.log(ofs, true);
			ofs << '\n';

			ofs.close();
		}
	}
	*/

	inline void Logger::message(std::string filename) {
		getInstance()->_kInfos.emplace_back(filename);
	}
}

#endif
