#ifndef ATOMICPOINTER_H
#define ATOMICPOINTER_H

#include <atomic>

class AtomicPointer {
public:
	AtomicPointer() {}
	explicit AtomicPointer(void *p) : rep_(p) {}
	inline void* Acquire_Load() const {
		return rep_.load(std::memory_order_acquire); 
	}
	inline void Release_Store(void* v) {
		rep_.store(v, std::memory_order_release);
	}
	inline void* NoBarrier_Load() const {
		return rep_.load(std::memory_order_relaxed);
	}
	inline void NoBarrier_Store(void* v) { 
		rep_.store(v, std::memory_order_relaxed); 
	}
private:
	std::atomic<void*> rep_;
};

#endif
