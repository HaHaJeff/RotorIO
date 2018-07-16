#include "IOStrategy/MPIIOStrategy.h"
#include <stdio.h>
#include <errno.h>
#include <string.h>

//<-------------------------MPIIOStrategy------------------------------->
MPIIOStrategy::MPIIOStrategy(const MPI_Comm& comm, int rank) : Strategy(), fh_(0), filename_(""), rank_(rank) {
	MPI_Comm_dup(comm, &comm_);
	MPI_Info_create(&info_);
}

//@override Write()
ssize_t MPIIOStrategy::Write(const Data_3D& data) {
	return 0;
}

//@override Write()
ssize_t MPIIOStrategy::Write(const Data_3D& data, bool formated) {
	return 0;
}

//@override Read()
ssize_t MPIIOStrategy::Read(Data_3D& data) {
	return 0;
}

//@override Read()
ssize_t MPIIOStrategy::Read(Data_3D& data, bool formated) {
	return 0;
}

//@override Lseek()
off_t MPIIOStrategy::Lseek(off_t off) {
	return 0;
}

//@override Open()
int MPIIOStrategy::Open(const std::string& filename) {
	return 0;
}

//@override Close()
int MPIIOStrategy::Close() {
	return 0;
}

int MPIIOStrategy::Split(int groupsize) {
	MPI_Comm_split(comm_, rank_/groupsize, rank_%groupsize, &comm_);
	MPI_Comm_rank(comm_, &rank_);
	return rank_;
}

int MPIIOStrategy::SetView(int offset, MPI_Datatype etype, MPI_Datatype ftype) {
	return MPI_File_set_view(fh_, (MPI_Offset)offset, etype, ftype, "native", info_);
}

int MPIIOStrategy::SetInfo(const std::string& key, const std::string& value) {
	MPI_Info_set(info_, key.c_str(), value.c_str());
}
//<-------------------------MPIIOStrategy------------------------------->



//<-------------------------OneFilePerProcessAllWrite-------------------------------->
MPIIOStrategyA::MPIIOStrategyA(const MPI_Comm& comm, int rank) : MPIIOStrategy(comm, rank) {
}

//@override Write()
ssize_t MPIIOStrategyA::Write(const Data_3D& data) {
	int count = data.GetCount();
	double *pData = data.pData_;
	ssize_t ret = MPI_File_write_all(fh_, pData, count, MPI_DOUBLE, NULL);
	if (ret != MPI_SUCCESS) fprintf(stderr, "write filename: %s\n, strerror: %s\n, ret: %d\n", filename_.c_str(), strerror(errno), ret);
	return ret;
}

//@override Write()
ssize_t MPIIOStrategyA::Write(const Data_3D& data, bool formated) {
	return 0;
}


//@override Read()
ssize_t MPIIOStrategyA::Read(Data_3D& data) {
	int count = data.GetCount();
	double *pData = data.pData_;
	ssize_t ret = MPI_File_read_all(fh_, pData, count, MPI_DOUBLE, NULL);
	if (ret != MPI_SUCCESS) fprintf(stderr, "read filename: %s\n, strerror: %s\n, ret: %d\n", filename_.c_str(), strerror(errno), ret);
	return ret;
}

//@override Write()
ssize_t MPIIOStrategyA::Read(Data_3D& data, bool formated) {
	return 0;
}


//@override Lseek()
off_t MPIIOStrategyA::Lseek(off_t off) {

	return 0;
}

//@override Open()
int MPIIOStrategyA::Open(const std::string& filename) {
	filename_ = filename;
	int ret = MPI_File_open(comm_, filename_.c_str(), MPI_MODE_CREATE | MPI_MODE_RDWR | MPI_MODE_UNIQUE_OPEN, MPI_INFO_NULL, &fh_);
	if (ret != MPI_SUCCESS) fprintf(stderr, "open filename: %s\n, strerror: %s\n, ret: %d\n", filename_.c_str(), strerror(errno), ret);
	return ret;
}

//@override Close()
int MPIIOStrategyA::Close() {
	filename_.clear();
	int ret = MPI_File_close(&fh_);
	if (ret != MPI_SUCCESS) fprintf(stderr, "open filename: %s\n, strerror: %s\n, ret: %d\n", filename_.c_str(), strerror(errno), ret);
	fh_ = 0;
	return ret;
}
//<-------------------------OneFilePerProcessAllWrite-------------------------------->


//<-------------------------SingleSharedFileOneWrites-------------------------------->
MPIIOStrategyB::MPIIOStrategyB(const MPI_Comm& comm, int rank) : MPIIOStrategy(comm, rank) {

}

//@override Write()
ssize_t MPIIOStrategyB::Write(const Data_3D& data) {
	return 0;
}

//@override Write()
ssize_t MPIIOStrategyB::Write(const Data_3D& data, bool formated) {

	return 0;
}

//@override Read()
ssize_t MPIIOStrategyB::Read(Data_3D& data) {
	return 0;
}

//@override Read()
ssize_t MPIIOStrategyB::Read(Data_3D& data, bool formated) {

	return 0;
}

//@override Lseek()
off_t MPIIOStrategyB::Lseek(off_t off) {

	return 0;
}

//@override Open()
int MPIIOStrategyB::Open(const std::string& filename) {
	return 0;
}

//@override Close()
int MPIIOStrategyB::Close() {
	return 0;
}
//<-------------------------SingleSharedFileOneWrites-------------------------------->


//<-------------------------SingleSharedFileAllWrite-------------------------------->
MPIIOStrategyC::MPIIOStrategyC(const MPI_Comm& comm, int rank) : MPIIOStrategy(comm, rank) {

}

//@override Write()
ssize_t MPIIOStrategyC::Write(const Data_3D& data) {

	int count = data.GetCount();
	double *pData = data.pData_;


	return 0;

}

//@override Write()
ssize_t MPIIOStrategyC::Write(const Data_3D& data, bool formated) {

	return 0;
}

//@override Read()
ssize_t MPIIOStrategyC::Read(Data_3D& data) {
	int count = data.GetCount();
	double *pData = data.pData_;

  return 0;
}

//@override Read()
ssize_t MPIIOStrategyC::Read(Data_3D& data, bool formated) {

	return 0;
}

//@override Lseek()
off_t MPIIOStrategyC::Lseek(off_t off) {

	return 0;
}

//@override Open()
int MPIIOStrategyC::Open(const std::string& filename) {

	filename_ = filename;
	int ret = MPI_File_open(comm_, filename_.c_str(), MPI_MODE_CREATE | MPI_MODE_RDWR | MPI_MODE_UNIQUE_OPEN, MPI_INFO_NULL, &fh_);
	if (ret != MPI_SUCCESS) fprintf(stderr, "open filename: %s\n, strerror: %s\n, ret: %d\n", filename_.c_str(), strerror(errno), ret);
	return ret;
}

//@override Close()
int MPIIOStrategyC::Close() {
	filename_.clear();
	int ret = MPI_File_close(&fh_);
	if (ret != MPI_SUCCESS) fprintf(stderr, "open filename: %s\n, strerror: %s\n, ret: %d\n", filename_.c_str(), strerror(errno), ret);
	fh_ = 0;
	return ret;
}
//<-------------------------SingleSharedFileAllWrite-------------------------------->


//<-------------------------SingleSharedFileSubsetWrite-------------------------------->


MPIIOStrategyD::MPIIOStrategyD(const MPI_Comm& comm, int rank) : MPIIOStrategy(comm, rank) {
}

//@override Write()
ssize_t MPIIOStrategyD::Write(const Data_3D& data) {
	return 0;
}

//@override Write()
ssize_t MPIIOStrategyD::Write(const Data_3D& data, bool formated) {
	return 0;
}

//@override Read()
ssize_t MPIIOStrategyD::Read(Data_3D& data) {
	return 0;
}

//@override Read()
ssize_t MPIIOStrategyD::Read(Data_3D& data, bool formated) {
	return 0;
}

//@override Lseek()
off_t MPIIOStrategyD::Lseek(off_t off) {
	return 0;
}

//@override Open()
int MPIIOStrategyD::Open(const std::string& filename) {
	return 0;
}

//@override Close()
int MPIIOStrategyD::Close() {
	return 0;
}

//<-------------------------SingleSharedFileSubsetWrite-------------------------------->

