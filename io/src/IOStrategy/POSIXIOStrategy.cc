#include "IOStrategy/POSIXIOStrategy.h"

#include <cstdint>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <cassert>

//<-------------------------POSIXIOStrategy------------------------------->
POSIXIOStrategy::POSIXIOStrategy() : Strategy(), fd_(-1), filename_("") {

}

//@override Write()
ssize_t POSIXIOStrategy::Write(const Data_3D& data) {
	return 0;
}

//@override Write()
ssize_t POSIXIOStrategy::Write(const Data_3D& data, bool formated) {
	return 0;
}

//@override Read()
ssize_t POSIXIOStrategy::Read(Data_3D& data) {
	return 0;
}

//@override Read()
ssize_t POSIXIOStrategy::Read(Data_3D& data, bool formated) {
	return 0;
}

//@override Lseek()
off_t POSIXIOStrategy::Lseek(off_t off) {
	return 0;
}

//@override Open()
int POSIXIOStrategy::Open(const std::string& filename) {
	return 0;
}

//@override Close()
int POSIXIOStrategy::Close() {
	return 0;
}
//<-------------------------POSIXIOStrategy------------------------------->


//<-------------------------POSIXIOStrategyOneFilePerProcessAllWrite-------------------------------->
POSIXIOStrategyOneFilePerProcessAllWrite::POSIXIOStrategyOneFilePerProcessAllWrite() : POSIXIOStrategy() {
}

//@override Write()
ssize_t POSIXIOStrategyOneFilePerProcessAllWrite::Write(const Data_3D& data) {
	double *ptr = data.pData_;

	int left = (data.nx_ - 1) * (data.ny_ - 1) * (data.nz_ - 1) * sizeof(double);
	int done = 0;
	int offset = 0;

	while (left > 0 && (done = write(fd_, ptr + offset, left - done)) != -1) {
		offset += done/sizeof(double);
		left -= done;
	}
	return offset;
}

ssize_t POSIXIOStrategyOneFilePerProcessAllWrite::Write(const Data_3D& data, bool formated) {
	if(formated == false) return Write(data);
	
	double* ptr = data.pData_;

	FILE* fPtr = fdopen(fd_, "w");

	assert(ptr != NULL);

	int left = data.GetCount();
	int done = 0;

	int offset = 0;

	
	//TODO: pricision need to set correclly;
	while (left > 0 && (fprintf(fPtr, "%.15lf %.15lf %.15lf", ptr[done], ptr[done+1], ptr[done+2]))){
		done += 3;
		offset += done*sizeof(double);
		left -= done;
	}

	if (NULL != fPtr){
		fclose(fPtr);
	}

	return offset;
}
//@override Read()
ssize_t POSIXIOStrategyOneFilePerProcessAllWrite::Read(Data_3D& data) {
	if (fd_ == -1) {
		Open(filename_);
	}

	double *ptr = data.pData_;
	int left = (data.nx_ - 1) * (data.ny_ - 1) * (data.nz_ - 1) * sizeof(double);
	int done = 0;
	int offset = 0;

	while (left > 0 && (done = read(fd_, ptr + offset, left - done)) != -1) {
		offset += done/sizeof(double);
		left -= done;
	}

	return offset;
}

//@override Read()
ssize_t POSIXIOStrategyOneFilePerProcessAllWrite::Read(Data_3D& data, bool formated) {
	return 0;
}

//@override Lseek()
off_t POSIXIOStrategyOneFilePerProcessAllWrite::Lseek(off_t off) {

	return 0;
}

//@override Open()
int POSIXIOStrategyOneFilePerProcessAllWrite::Open(const std::string& filename) {
	filename_ = filename;
	fd_ = open(filename_.c_str(), O_CREAT | O_WRONLY, 0664);
	return fd_;
}

//@override Close()
int POSIXIOStrategyOneFilePerProcessAllWrite::Close() {
	int ret = 0;
	filename_.clear();
	
	if (fd_ != -1) {
		ret = close(fd_);
	}

	return ret;
}
//<-------------------------POSIXIOStrategyOneFilePerProcessAllWrite-------------------------------->



//<-------------------------POSIXIOStrategySingleSharedFileOneWrites-------------------------------->
POSIXIOStrategySingleSharedFileOneWrites::POSIXIOStrategySingleSharedFileOneWrites() : POSIXIOStrategy() {

}

//@override Write()
ssize_t POSIXIOStrategySingleSharedFileOneWrites::Write(const Data_3D& data) {

	return 0;
}

//@override Write()
ssize_t POSIXIOStrategySingleSharedFileOneWrites::Write(const Data_3D& data, bool formated) {

	return 0;
}

//@override Read()
ssize_t POSIXIOStrategySingleSharedFileOneWrites::Read(Data_3D& data) {

	return 0;
}

//@override Read()
ssize_t POSIXIOStrategySingleSharedFileOneWrites::Read(Data_3D& data, bool formated) {

	return 0;
}

//@override Lseek()
off_t POSIXIOStrategySingleSharedFileOneWrites::Lseek(off_t off) {

	return 0;
}

//@override Open()
int POSIXIOStrategySingleSharedFileOneWrites::Open(const std::string& filename) {

	return 0;
}

//@override Close()
int POSIXIOStrategySingleSharedFileOneWrites::Close() {

	return 0;
}
//<-------------------------POSIXIOStrategySingleSharedFileOneWrites-------------------------------->




//<-------------------------POSIXIOStrategySingleSharedFileAllWrite-------------------------------->
POSIXIOStrategySingleSharedFileAllWrite::POSIXIOStrategySingleSharedFileAllWrite() : POSIXIOStrategy() {

}

//@override Write()
ssize_t POSIXIOStrategySingleSharedFileAllWrite::Write(const Data_3D& data) {

	return 0;
}

//@override Write()
ssize_t POSIXIOStrategySingleSharedFileAllWrite::Write(const Data_3D& data, bool formated) {

	return 0;
}

//@override Read()
ssize_t POSIXIOStrategySingleSharedFileAllWrite::Read(Data_3D& data) {

	return 0;
}

//@override Read()
ssize_t POSIXIOStrategySingleSharedFileAllWrite::Read(Data_3D& data, bool formated) {

	return 0;
}

//@override Lseek()
off_t POSIXIOStrategySingleSharedFileAllWrite::Lseek(off_t off) {

	return 0;
}

//@override Open()
int POSIXIOStrategySingleSharedFileAllWrite::Open(const std::string& filename) {

	return 0;
}

//@override Close()
int POSIXIOStrategySingleSharedFileAllWrite::Close() {

	return 0;
}
//<-------------------------POSIXIOStrategySingleSharedFileAllWrite-------------------------------->



//<-------------------------POSIXIOStrategySingleSharedFileSubsetWrite-------------------------------->
POSIXIOStrategySingleSharedFileSubsetWrite::POSIXIOStrategySingleSharedFileSubsetWrite() : POSIXIOStrategy() {

}

//@override Write()
ssize_t POSIXIOStrategySingleSharedFileSubsetWrite::Write(const Data_3D& data) {

	return 0;
}

//@override Write()
ssize_t POSIXIOStrategySingleSharedFileSubsetWrite::Write(const Data_3D& data, bool formated) {

	return 0;
}

//@override Read()
ssize_t POSIXIOStrategySingleSharedFileSubsetWrite::Read(Data_3D& data) {

	return 0;
}

//@override Read()
ssize_t POSIXIOStrategySingleSharedFileSubsetWrite::Read(Data_3D& data, bool formated) {

	return 0;
}

//@override Lseek()
off_t POSIXIOStrategySingleSharedFileSubsetWrite::Lseek(off_t off) {

	return 0;
}

//@override Open()
int POSIXIOStrategySingleSharedFileSubsetWrite::Open(const std::string& filename) {

	return 0;
}

//@override Close()
int POSIXIOStrategySingleSharedFileSubsetWrite::Close() {

	return 0;
}
//<-------------------------POSIXIOStrategySingleSharedFileSubsetWrite-------------------------------->




