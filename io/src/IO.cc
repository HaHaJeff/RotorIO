
#include "IO.h"
/*
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
	Mutex mu_;
	AtomicPoiter allowed_;

	intptr_t GetAllowed() const {
		return reinterpret_cast<intptr_t>(allowed_.Acquire_Load());
	}

	void SetAllowed(intptr_t v) {
		allowed_.Release_Store(reinterpret_cast<void*>(v));
	}
};
*/
class PosixSequentialFile : public SequentialFile {
public:
	PosixSequentialFile(const std::string& fname, int fd) : filename_(fname), fd_(fd) {}
	virtual ~PosixSequentialFile() { close(fd_); }

	virtual bool Read(size_t n, void* out) {
		bool ret = false;
		while (true) {
			ssize_t r = read(fd_, out, n);
			if (r < 0) {
				if (errno == EINTR) {
					continue;
				}

				break;
			}
			ret = true;
			break;
		}
		return ret;
	}

	virtual bool Skip(uint64_t n) {
		if (lseek(fd_, n, SEEK_CUR) == static_cast<off_t>(-1)) {
			return false;
		}
		return true;
	}
private:
	std::string filename_;
	int fd_;

};

class PosixRandomAccessFile: public RandomAccessFile {
public:
	PosixRandomAccessFile(const std::string& fname, int fd, Limiter* limiter)
		: filename_(fname), fd_(fd), limiter_(limiter) {
		temporary_fd_ = !limiter->Acquire();
		if (temporary_fd_)	{
			close(fd_);
			fd_ = -1;
		}
	}

	virtual ~PosixRandomAccessFile() {
		if (!temporary_fd_) {
			close(fd_);
			limiter_->Release();
		}
	}

	virtual bool Read(uint64_t offset, size_t n, void* out) {
		int fd = fd_;
		if (temporary_fd_) {
			fd = open(filename_.c_str(), O_RDONLY);
			if (fd < 0) {
				return false;
			}
		}
		ssize_t r = pread(fd, out, n, static_cast<off_t>(offset));

		if (r < 0) {
			return false;
		}

		if (temporary_fd_) {
			close(fd);
		}

		return true;
	}
private:
    std::string filename_;
	int fd_;
	bool temporary_fd_;
	Limiter* limiter_;
};
