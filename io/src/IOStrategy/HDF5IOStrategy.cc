#include "HDF5IOStrategy.h" 
#include <iostream>
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
/*
   void HDF5IOStrategy::SetDimStride(int nx, int ny, int nz) {

   }

   void HDF5IOStrategy::SetDataspace(const hsize_t dimsf[3], const hsize_t chunk_dims[3]) {

   }
   void HDF5IOStrategy::SetDatasetid(const hsize_t chunk_dims[3], const std::string& dataname) {

   }
   void HDF5IOStrategy::SetChunk(int blockid, const hsize_t chunk_dims[3], const hsize_t chunk_count[3]) {

   }
   */

void HDF5IOStrategy::SetDataspace(const hsize_t dimsf[3], const hsize_t chunk_dims[3]) {
	dataspace_ = H5Screate_simple(3, dimsf, NULL);
	chunkspace_ = H5Screate_simple(3, chunk_dims, NULL);
}

void HDF5IOStrategy::SetDatasetid(const hsize_t chunk_dims[3], const std::string& dataname) {
	hid_t plist_id = H5Pcreate(H5P_DATASET_CREATE);
	H5Pset_chunk(plist_id, 3, chunk_dims);
	datasetid_ = H5Dcreate(fileid_, dataname.c_str(), H5T_NATIVE_DOUBLE, dataspace_, H5P_DEFAULT, plist_id, H5P_DEFAULT);
	H5Pclose(plist_id);
}

/*chunk 33 65 193, x = 2, y = 2, z = 3*/
void HDF5IOStrategy::SetChunk(int blockid, const hsize_t chunk_dims[3], const hsize_t chunk_count[3]) {
	hsize_t count[3]{1, 1, 1}, stride[3]{1, 1, 1}, block[3]{chunk_dims[0], chunk_dims[1], chunk_dims[2]};
	int nx = chunk_count[0], ny = chunk_count[1], nz = chunk_count[2];
	//FIXME: on z direction, topo is blockid/(nx*ny)
	hsize_t offset[3]{((blockid/ny)%nx)*chunk_dims[0], (blockid%ny)*chunk_dims[1], (blockid/(nx*ny))*chunk_dims[2]};
	dataspace_ = H5Dget_space(datasetid_);
	H5Sselect_hyperslab(dataspace_, H5S_SELECT_SET, offset, stride, count, block);
}

void HDF5IOStrategy::SetView(int block_id, const hsize_t file_dims[3], const hsize_t chunk_dims[3], const std::string& name) {

	hsize_t chunk_count[3] = {file_dims[0]/chunk_dims[0], file_dims[1]/chunk_dims[1], file_dims[2]/chunk_dims[2]};
	SetDataspace(file_dims, chunk_dims);
	SetDatasetid(chunk_dims, name);
	SetChunk(block_id, chunk_dims, chunk_count);
}

//<-------------------------HDF5IOStrategy------------------------------->


//<-------------------------OneFilePerProcessAllWrite-------------------------------->
HDF5IOStrategyA::HDF5IOStrategyA(const MPI_Comm& comm, int rank) : HDF5IOStrategy(comm, rank) {
}

//@override Write()
ssize_t HDF5IOStrategyA::Write(const Data_3D& data) {
	double *pData = data.pData_;

	std::string name = "/"+data.name_;
	/*
	   hsize_t dims[3]{data.nx_, data.ny_, data.nz_};
	   hid_t dataspace_id = H5Screate_simple(3, dims, NULL);
	   hid_t dataset_id = H5Dcreate(fileid_, name.c_str(), H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
	   */
	H5Dwrite(datasetid_, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, pData);
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
	H5Sclose(dataspace_);
	H5Sclose(chunkspace_);
	H5Dclose(datasetid_);
	H5Fclose(fileid_);
	fileid_ = 0;
	return fileid_;
}
//<-------------------------OneFilePerProcessAllWrite-------------------------------->


//<-------------------------SingleSharedFileOneWrites-------------------------------->
HDF5IOStrategyB::HDF5IOStrategyB(const MPI_Comm& comm, int rank) : HDF5IOStrategy(comm, rank) {
}

//@override Write()
//FIXME: need to change topo to nx ny nz and modify offset
ssize_t HDF5IOStrategyB::Write(const Data_3D& data) {
	double *pData = data.pData_;
	/*
	   hsize_t chunk_dims[3]{data.nx_, data.ny_, data.nz_};
	   hsize_t chunk_count[3]{nx_, ny_, nz_};
	   hsize_t dimsf[3]{data.nx_*nx_, data.ny_*ny_, data.nz_*nz_};

	   int blockid = data.blockid_;
	   std::string dataname = "/"+data.name_;

	   SetDataspace(dimsf, chunk_dims);
	   SetDatasetid(chunk_dims, dataname);
	   SetChunk(blockid, chunk_dims, chunk_count);
	   */
	hsize_t plist_id = H5Pcreate(H5P_DATASET_XFER);
	H5Pset_dxpl_mpio(plist_id, H5FD_MPIO_COLLECTIVE);
	H5Dwrite(datasetid_, H5T_NATIVE_DOUBLE, chunkspace_, dataspace_, plist_id, pData);
	H5Pclose(plist_id);

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

	filename_ = filename;

	MPI_Info info = MPI_INFO_NULL;

	hid_t plist_id = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plist_id, comm_, info);

	fileid_ = H5Fcreate(filename_.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plist_id);

	H5Pclose(plist_id);
	return fileid_;

}

//@override Close()
int HDF5IOStrategyB::Close() {

	filename_.clear();
	H5Sclose(dataspace_);
	H5Sclose(chunkspace_);
	H5Dclose(datasetid_);
	H5Fclose(fileid_);
	fileid_ = 0;
	return fileid_;

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

void HDF5IOStrategyC::DatasetSetChunks(int nCols, int size, int maxRows, hid_t plistDCreate) {
	hsize_t cdims[2];
	int nRows;

	int nDims = 2;
	if (nCols == 1) {
		nDims = 1;
	}

	//	if ()
}

ssize_t HDF5IOStrategyC::WriteVector(const std::vector<Data_3D>& data, int dims[3], int nProc) {
	char datasetName[80];
	/*
	   hsize_t dimsf[2];
	   hid_t fileSpace;
	   hid_t dsetID;
	   hid_t plistID;
	   hid_t plistDCreate;

	   int nProc = 432;
	   int numPoints = data[0].GetCount();
	   for (int i = 0; i < nProc; i++) {
	//三维数据，二维存储
	dimsf[0] = numPoints;
	dimsf[1] = 3;
	fileSpace = H5Screate_simple(2, dimsf, NULL);

	plistID = H5Pcreate(H5P_LINK_CREATE);
	H5Pset_create_intermediate_group(plistID, 1);

	plistDCreate = H5Pcreate(H5P_DATASET_CREATE);
	DatasetSetChunks(3, sizeof(double), numPoints, plistDCreate);

	sprintf(datasetName, "MESH/processor%d/POINTS", i);

	dsetID = H5Dcreate2(fileid_, datasetName, H5T_NATIVE_DOUBLE, fileSpace, plistID, plistDCreate, H5P_DEFAULT);

	H5Dclose(dsetID);
	H5Pclose(plistID);
	H5Pclose(plistDCreate);
	H5Sclose(fileSpace);
	}
	double pointList[numPoints][3];
	double *x = data[0].pData_;
	double *y = data[1].pData_;
	double *z = data[2].pData_;
	for (int i = 0; i < numPoints; i++) {
	pointList[i][0] = x[i];
	pointList[i][1] = y[i];
	pointList[i][2] = z[i];
	}
	*/


	hsize_t dimsf[3];
	hid_t fileSpace;
	hid_t dsetID;
	hid_t plistID;
	hid_t plistDCreate;

	int numPoints = data[0].GetCount();

	for (int j = 0; j < data.size(); j++) {
		for (int i = 0; i < nProc; i++) {
			dimsf[0] = dims[0];
			dimsf[1] = dims[1];
			dimsf[2] = dims[2];


			fileSpace = H5Screate_simple(3, dimsf, NULL);

			plistID = H5Pcreate(H5P_LINK_CREATE);
			H5Pset_create_intermediate_group(plistID, 1);

			plistDCreate = H5Pcreate(H5P_DATASET_CREATE);
			DatasetSetChunks(3, sizeof(double), numPoints, plistDCreate);

			sprintf(datasetName, "Mesh/%s", data[j].name_.c_str());

			dsetID = H5Dcreate2(fileid_, datasetName, H5T_NATIVE_DOUBLE, fileSpace, plistID, plistDCreate, H5P_DEFAULT);

			H5Dclose(dsetID);
			H5Pclose(plistID);
			H5Pclose(plistDCreate);
			H5Sclose(fileSpace);
		}
	}

/*
	if (rank_ == 0) {
		std::cout << "start write" << std::endl;
	}
	*/

	for (int i = 0; i < data.size(); i++ ) {
		sprintf(datasetName, "Mesh/%s", data[i].name_.c_str());
		/*
		if (rank_ == 0) {
			std::cout << "start open " << datasetName << std::endl;
		}
		*/
		dsetID = H5Dopen2(fileid_, datasetName, H5P_DEFAULT);

		plistID = H5Pcreate(H5P_DATASET_XFER);
		H5Pset_dxpl_mpio(plistID, H5FD_MPIO_COLLECTIVE);
		/*
		if (rank_ == 0) {
			std::cout << "start write " << datasetName << std::endl;
		}
		*/
		H5Dwrite(dsetID, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, plistID, data[i].pData_);

		H5Pclose(plistID);
		H5Dclose(dsetID);
	}
}

//@override Read()
ssize_t HDF5IOStrategyC::Read(Data_3D& data) {

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
	//	char dataFile[80] = "h5Data.h5";

	filename_ = filename;
	hid_t plistID = H5Pcreate(H5P_FILE_ACCESS);
	H5Pset_fapl_mpio(plistID, comm_, MPI_INFO_NULL);

	fileid_ = H5Fcreate(filename.c_str(), H5F_ACC_TRUNC, H5P_DEFAULT, plistID);
	H5Pclose(plistID);
	return 0;
}

//@override Close()
int HDF5IOStrategyC::Close() {

	filename_.clear();
	return 0;
}
//<-------------------------SingleSharedFileOneWrite-------------------------------->



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

