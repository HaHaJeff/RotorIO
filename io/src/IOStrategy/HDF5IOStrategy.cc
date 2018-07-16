#include "HDF5IOStrategy.h"

//<-------------------------HDF5IOStrategy------------------------------->
HDF5IOStrategy::HDF5IOStrategy(const MPI_Comm& comm, int rank) : Strategy(), fileid_(0), filename_(""), rank_(rank) {
	MPI_Comm_dup(comm, &comm_);
}

//@override Write()
ssize_t HDF5IOStrategy::Write(const Data_3D& data) {
	return 0;
}

//@override Write()
ssize_t HDF5IOStrategy::Write(const Data_3D& data, bool formated) {
	return 0;
}

//@override Read()
ssize_t HDF5IOStrategy::Read(Data_3D& data) {
	return 0;
}

//@override Read()
ssize_t HDF5IOStrategy::Read(Data_3D& data, bool formated) {
	return 0;
}

//@override Lseek()
off_t HDF5IOStrategy::Lseek(off_t off) {
	return 0;
}

//@override Open()
int HDF5IOStrategy::Open(const std::string& filename) {
	return 0;
}

//@override Close()
int HDF5IOStrategy::Close() {
	return 0;
}
//<-------------------------HDF5IOStrategy------------------------------->


//<-------------------------OneFilePerProcessAllWrite-------------------------------->
HDF5IOStrategyA::HDF5IOStrategyA(const MPI_Comm& comm, int rank) : HDF5IOStrategy(comm, rank) {
}

//@override Write()
ssize_t HDF5IOStrategyA::Write(const Data_3D& data) {
	int count = data.GetCount();
	double *pData = data.pData_;

  std::string name = "/"+data.name_;
  hsize_t dims[3]{data.nx_, data.ny_, data.nz_};
  hid_t dataspace_id = H5Screate_simple(3, dims, NULL);
  hid_t dataset_id = H5Dcreate(fileid_, name.c_str(), H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
  H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pData);

  return 0;
}

//@override Write()
ssize_t HDF5IOStrategyA::Write(const Data_3D& data, bool formated) {
	int count = data.GetCount();
	double *pData = data.pData_;

	return 0;
}


//@override Read()
ssize_t HDF5IOStrategyA::Read(Data_3D& data) {
	int count = data.GetCount();
	double *pData = data.pData_;
	return 0;
}

//@override Write()
ssize_t HDF5IOStrategyA::Read(Data_3D& data, bool formated) {
	return 0;
}


//@override Lseek()
off_t HDF5IOStrategyA::Lseek(off_t off) {

	return 0;
}

//@override Open()
int HDF5IOStrategyA::Open(const std::string& filename) {
	filename_ = filename;
  fileid_   = H5Fcreate(filename_.c_str(), H5F_ACC_EXCL, H5P_DEFAULT, H5P_DEFAULT);
	return fileid_;
}

//@override Close()
int HDF5IOStrategyA::Close() {
	filename_.clear();
  H5Fclose(fileid_);
  fileid_ = 0;
	return fileid_;
}
//<-------------------------OneFilePerProcessAllWrite-------------------------------->


//<-------------------------SingleSharedFileOneWrites-------------------------------->
HDF5IOStrategyB::HDF5IOStrategyB(const MPI_Comm& comm, int rank) : HDF5IOStrategy(comm, rank) {
}

//@override Write()
ssize_t HDF5IOStrategyB::Write(const Data_3D& data) {
	return 0;
}

//@override Write()
ssize_t HDF5IOStrategyB::Write(const Data_3D& data, bool formated) {

	return 0;
}

//@override Read()
ssize_t HDF5IOStrategyB::Read(Data_3D& data) {
	return 0;
}

//@override Read()
ssize_t HDF5IOStrategyB::Read(Data_3D& data, bool formated) {

	return 0;
}

//@override Lseek()
off_t HDF5IOStrategyB::Lseek(off_t off) {

	return 0;
}

//@override Open()
int HDF5IOStrategyB::Open(const std::string& filename) {
	return 0;
}

//@override Close()
int HDF5IOStrategyB::Close() {
	return 0;
}
//<-------------------------SingleSharedFileOneWrites-------------------------------->

//<-------------------------SingleSharedFileAllWrite-------------------------------->
HDF5IOStrategyC::HDF5IOStrategyC(const MPI_Comm& comm, int rank) : HDF5IOStrategy(comm, rank) {

}

//@override Write()
ssize_t HDF5IOStrategyC::Write(const Data_3D& data) {
	int count = data.GetCount();
	double *pData = data.pData_;
	return 0;

}

//@override Write()
ssize_t HDF5IOStrategyC::Write(const Data_3D& data, bool formated) {

	return 0;
}

//@override Read()
ssize_t HDF5IOStrategyC::Read(Data_3D& data) {
	int count = data.GetCount();
	double *pData = data.pData_;
	return 0;

}

//@override Read()
ssize_t HDF5IOStrategyC::Read(Data_3D& data, bool formated) {

	return 0;
}

//@override Lseek()
off_t HDF5IOStrategyC::Lseek(off_t off) {

	return 0;
}

//@override Open()
int HDF5IOStrategyC::Open(const std::string& filename) {

	filename_ = filename;
	return 0;
}

//@override Close()
int HDF5IOStrategyC::Close() {

	filename_.clear();
	return 0;
}
//<-------------------------SingleSharedFileAllWrite-------------------------------->


//<-------------------------SingleSharedFileSubsetWrite-------------------------------->
HDF5IOStrategyD::HDF5IOStrategyD(const MPI_Comm& comm, int rank) : HDF5IOStrategy(comm, rank) {

}
//@override Write()
ssize_t HDF5IOStrategyD::Write(const Data_3D& data) {
	return 0;
}

//@override Write()
ssize_t HDF5IOStrategyD::Write(const Data_3D& data, bool formated) {
	return 0;
}

//@override Read()
ssize_t HDF5IOStrategyD::Read(Data_3D& data) {
	return 0;
}

//@override Read()
ssize_t HDF5IOStrategyD::Read(Data_3D& data, bool formated) {
	return 0;
}

//@override Lseek()
off_t HDF5IOStrategyD::Lseek(off_t off) {
	return 0;
}

//@override Open()
int HDF5IOStrategyD::Open(const std::string& filename) {
	return 0;
}

//@override Close()
int HDF5IOStrategyD::Close() {
	return 0;
}

//<-------------------------SingleSharedFileSubsetWrite-------------------------------->

