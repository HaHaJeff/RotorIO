#ifndef IO_H
#define IO_H

#include "Mutex.h"
#include "AtomicPointer.h"

#include <stdarg.h>
#include <stdint.h>
#include <string>
#include <sys/types.h>
#include <unistd.h>
#include <sys/stat.h>
#include <fcntl.h>


class Limiter {
public:
	Limiter(intptr_t n) {
		SetAllowed(n);
	}

	bool Acquire() {
		if (GetAllowed() <= 0) {
			return false;
		}

		MutexLock l(&mu_);
		intptr_t x = GetAllowed();
		if (x <= 0) {
			return false;
		} else {
			SetAllowed(x - 1);
			return true;
		}
	}

	void Release() {
		MutexLock l(&mu_);
		SetAllowed(GetAllowed() + 1);
	}

private:
	void SetAllowed(intptr_t v) {
		allowed_.Release_Store(reinterpret_cast<void*>(v));
	}
	intptr_t GetAllowed() const {
		return reinterpret_cast<intptr_t>(allowed_.Acquire_Load());
	}
	Mutex mu_;
	AtomicPointer allowed_;
};

class SequentialFile {
public:
	SequentialFile() {}
	virtual ~SequentialFile();

	// Read n bytes from file.
	// REQUIRES: External synchronizaion
	virtual bool Read(size_t n, void* out) = 0;
	virtual bool Skip(uint64_t n) = 0;
private:
	//No copying allowed
	SequentialFile(const SequentialFile&);
	void operator=(const SequentialFile&);
};

class RandomAccessFile {
public:
	RandomAccessFile() {}
	virtual ~RandomAccessFile();

	//Safe for concurrent use by multipe threads
	virtual bool Read(uint64_t offset, size_t n, void* out) = 0;
private:
	//No copying allowed
	RandomAccessFile(const RandomAccessFile&);
	void operator=(const RandomAccessFile&);
};

class WritableFile {
public:
	WritableFile() {}
	virtual ~WritableFile();

	virtual bool Append(const void* buf) = 0;
	virtual bool Close() = 0;
	virtual bool Flush() = 0;
	virtual bool Sync()  = 0;

private:
	//No copying allowed
	WritableFile(const WritableFile&);
	void operator=(const WritableFile&);
};

class Logger {
public:
	Logger() {}
	virtual ~Logger();
	virtual void Logv(const char* format, va_list ap) = 0;

private:
	//No copying allowed
	Logger(const Logger&);
	void operator=(const Logger&);
};

class FileLock {
public:
	FileLock() {}
	virtual ~FileLock();

private:
	FileLock(const FileLock&);
	void operator=(const FileLock&);
};

#endif
