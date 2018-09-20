#ifndef __SINGLETON_H_
#define __SINGLETON_H_

#include <assert.h>
// http ://www.boost.org/doc/libs/1_39_0/boost/pool/detail/singleton.hpp

namespace mn {

	/*
	*	@note	Singleton
	*/
	// T must be: no-throw default constructible and no-throw destructible
	template <typename T>
	struct Singleton {
	private:
		struct SingletonCreator {
			// This constructor does nothing more than ensure that instance()
			//  is called before main() begins, thus creating the static
			//  T object before multithreading race issues can come up.
			SingletonCreator() { Singleton<T>::instance(); }
			inline void doNothing() const {}
		};
		static SingletonCreator g_Creator;

		// private on purpose
		Singleton();

	public:
		// If, at any point (in user code), singleton_default<T>::instance()
		//  is called, then the following function is instantiated.
		static T & instance() {
			// This is the object that we return a reference to.
			// It is guaranteed to be created before main() begins because of
			//  the next line.
			static T _instance;

			// The following line does nothing else than force the instantiation
			//  of singleton_default<T>::create_object, whose constructor is
			//  called before main() begins.
			g_Creator.do_nothing();

			return _instance;
		}
	};
	template <typename T>
	typename Singleton<T>::SingletonCreator
		Singleton<T>::g_Creator;

	/*
	*	\class	ManagedSingleton
	*	\note	Used in global systems
	*/
	template <typename T>
	struct ManagedSingleton {

		static void startup() {
			assert(_pInstance == nullptr);
			_pInstance = new T();
		}
		static void shutdown() {
			assert(_pInstance != nullptr);
			delete _pInstance;
		}

	protected:
		static T* _pInstance;

	public:
		/// 
		T* operator ->() { return _pInstance; }

		static T* getInstance() {
			assert(_pInstance != nullptr);
			return _pInstance;
		}
		static const T* getConstInstance() {
			assert(_pInstance != nullptr);
			return _pInstance;
		}
	};
	template <typename T>
	T*	ManagedSingleton<T>::_pInstance = nullptr;
}

#endif