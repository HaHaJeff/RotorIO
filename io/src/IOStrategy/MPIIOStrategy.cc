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



//<-------------------------MPIIOStrategyOneFilePerProcessAllWrite-------------------------------->
MPIIOStrategyOneFilePerProcessAllWrite::MPIIOStrategyOneFilePerProcessAllWrite(const MPI_Comm& comm, int rank) : MPIIOStrategy(comm, rank) {
}

//@override Write()
ssize_t MPIIOStrategyOneFilePerProcessAllWrite::Write(const Data_3D& data) {
	int count = data.GetCount();
	double *pData = data.pData_;
	ssize_t ret = MPI_File_write_all(fh_, pData, count, MPI_DOUBLE, NULL);
	if (ret != MPI_SUCCESS) fprintf(stderr, "write filename: %s\n, strerror: %s\n, ret: %d\n", filename_.c_str(), strerror(errno), ret);
	return ret;
}

//@override Write()
ssize_t MPIIOStrategyOneFilePerProcessAllWrite::Write(const Data_3D& data, bool formated) {
	return 0;
}


//@override Read()
ssize_t MPIIOStrategyOneFilePerProcessAllWrite::Read(Data_3D& data) {
	int count = data.GetCount();
	double *pData = data.pData_;
	ssize_t ret = MPI_File_read_all(fh_, pData, count, MPI_DOUBLE, NULL);
	if (ret != MPI_SUCCESS) fprintf(stderr, "read filename: %s\n, strerror: %s\n, ret: %d\n", filename_.c_str(), strerror(errno), ret);
	return ret;
}

//@override Write()
ssize_t MPIIOStrategyOneFilePerProcessAllWrite::Read(Data_3D& data, bool formated) {
	return 0;
}


//@override Lseek()
off_t MPIIOStrategyOneFilePerProcessAllWrite::Lseek(off_t off) {

	return 0;
}

//@override Open()
int MPIIOStrategyOneFilePerProcessAllWrite::Open(const std::string& filename) {
	filename_ = filename;
	char buf[32] = {};
	int len = 32;
	int ret = MPI_File_open(comm_, filename_.c_str(), MPI_MODE_CREATE | MPI_MODE_RDWR | MPI_MODE_UNIQUE_OPEN, MPI_INFO_NULL, &fh_);
	if (ret != MPI_SUCCESS) fprintf(stderr, "open filename: %s\n, strerror: %s\n, ret: %d\n", filename_.c_str(), strerror(errno), ret);
	return ret;
}

//@override Close()
int MPIIOStrategyOneFilePerProcessAllWrite::Close() {
	filename_.clear();
	int ret = MPI_File_close(&fh_);
	if (ret != MPI_SUCCESS) fprintf(stderr, "open filename: %s\n, strerror: %s\n, ret: %d\n", filename_.c_str(), strerror(errno), ret);
	fh_ = 0;
	return ret;
}
//<-------------------------MPIIOStrategyOneFilePerProcessAllWrite-------------------------------->


//<-------------------------MPIIOStrategySingleSharedFileOneWrites-------------------------------->
MPIIOStrategySingleSharedFileOneWrites::MPIIOStrategySingleSharedFileOneWrites(const MPI_Comm& comm, int rank) : MPIIOStrategy(comm, rank) {

}

//@override Write()
ssize_t MPIIOStrategySingleSharedFileOneWrites::Write(const Data_3D& data) {
	return 0;
}

//@override Write()
ssize_t MPIIOStrategySingleSharedFileOneWrites::Write(const Data_3D& data, bool formated) {

	return 0;
}

//@override Read()
ssize_t MPIIOStrategySingleSharedFileOneWrites::Read(Data_3D& data) {
	return 0;
}

//@override Read()
ssize_t MPIIOStrategySingleSharedFileOneWrites::Read(Data_3D& data, bool formated) {

	return 0;
}

//@override Lseek()
off_t MPIIOStrategySingleSharedFileOneWrites::Lseek(off_t off) {

	return 0;
}

//@override Open()
int MPIIOStrategySingleSharedFileOneWrites::Open(const std::string& filename) {
	return 0;
}

//@override Close()
int MPIIOStrategySingleSharedFileOneWrites::Close() {
	return 0;
}
//<-------------------------MPIIOStrategySingleSharedFileOneWrites-------------------------------->


//<-------------------------MPIIOStrategySingleSharedFileAllWrite-------------------------------->
MPIIOStrategySingleSharedFileAllWrite::MPIIOStrategySingleSharedFileAllWrite(const MPI_Comm& comm, int rank) : MPIIOStrategy(comm, rank) {
	
}

//@override Write()
ssize_t MPIIOStrategySingleSharedFileAllWrite::Write(const Data_3D& data) {
				
	int count = data.GetCount();
	double *pData = data.pData_;
	ssize_t ret = MPI_File_write_all(fh_, pData, count, MPI_DOUBLE, NULL);
	if (ret != MPI_SUCCESS) fprintf(stderr, "write filename: %s\n, strerror: %s\n, ret: %d\n", filename_.c_str(), strerror(errno), ret);
	return ret;

}

//@override Write()
ssize_t MPIIOStrategySingleSharedFileAllWrite::Write(const Data_3D& data, bool formated) {

	return 0;
}

//@override Read()
ssize_t MPIIOStrategySingleSharedFileAllWrite::Read(Data_3D& data) {
	int count = data.GetCount();
	double *pData = data.pData_;
	ssize_t ret = MPI_File_read_all(fh_, pData, count, MPI_DOUBLE, NULL);
	if (ret != MPI_SUCCESS) fprintf(stderr, "read filename: %s\n, strerror: %s\n, ret: %d\n", filename_.c_str(), strerror(errno), ret);
	return ret;

}

//@override Read()
ssize_t MPIIOStrategySingleSharedFileAllWrite::Read(Data_3D& data, bool formated) {

	return 0;
}

//@override Lseek()
off_t MPIIOStrategySingleSharedFileAllWrite::Lseek(off_t off) {

	return 0;
}

//@override Open()
int MPIIOStrategySingleSharedFileAllWrite::Open(const std::string& filename) {

	filename_ = filename;
	char buf[32] = {};
	int len = 32;
	int ret = MPI_File_open(comm_, filename_.c_str(), MPI_MODE_CREATE | MPI_MODE_RDWR | MPI_MODE_UNIQUE_OPEN, MPI_INFO_NULL, &fh_);
	if (ret != MPI_SUCCESS) fprintf(stderr, "open filename: %s\n, strerror: %s\n, ret: %d\n", filename_.c_str(), strerror(errno), ret);
	return ret;
}

//@override Close()
int MPIIOStrategySingleSharedFileAllWrite::Close() {

	filename_.clear();
	int ret = MPI_File_close(&fh_);
	if (ret != MPI_SUCCESS) fprintf(stderr, "open filename: %s\n, strerror: %s\n, ret: %d\n", filename_.c_str(), strerror(errno), ret);
	fh_ = 0;
	return ret;
}
//<-------------------------MPIIOStrategySingleSharedFileAllWrite-------------------------------->


//<-------------------------MPIIOStrategySingleSharedFileSubsetWrite-------------------------------->


MPIIOStrategySingleSharedFileSubsetWrite::MPIIOStrategySingleSharedFileSubsetWrite(const MPI_Comm& comm, int rank) : MPIIOStrategy(comm, rank) {

}
//@override Write()
ssize_t MPIIOStrategySingleSharedFileSubsetWrite::Write(const Data_3D& data) {
	return 0;
}

//@override Write()
ssize_t MPIIOStrategySingleSharedFileSubsetWrite::Write(const Data_3D& data, bool formated) {
	return 0;
}

//@override Read()
ssize_t MPIIOStrategySingleSharedFileSubsetWrite::Read(Data_3D& data) {
	return 0;
}

//@override Read()
ssize_t MPIIOStrategySingleSharedFileSubsetWrite::Read(Data_3D& data, bool formated) {
	return 0;
}

//@override Lseek()
off_t MPIIOStrategySingleSharedFileSubsetWrite::Lseek(off_t off) {
	return 0;
}

//@override Open()
int MPIIOStrategySingleSharedFileSubsetWrite::Open(const std::string& filename) {
	return 0;
}


//@override Close()
int MPIIOStrategySingleSharedFileSubsetWrite::Close() {
	return 0;
}

//<-------------------------MPIIOStrategySingleSharedFileSubsetWrite-------------------------------->

