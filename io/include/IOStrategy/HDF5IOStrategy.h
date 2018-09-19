#ifndef HDF5IOSTRATEGY_H
#define HDF5IOSTRATEGY_H

#include "hdf5.h"
#include "Strategy.h"

#include <memory>
#include <vector>

const int g_idim=3;

struct HDF5Data{

	void SetDsetid(hid_t fileid);
	void SetFilespace();
	void SetMemspace();

	hid_t dsetid_;
	hid_t filespace_;
	hid_t memspace_;

	int globaldims_[g_idim];
	int chunkdims_[g_idim];
	int offset_[g_idim];
	int stride_[g_idim];
	int count_[g_idim];
	int block_[g_idim];
};

class HDF5IOStrategy : public Strategy {
public:
	HDF5IOStrategy(const MPI_Comm& comm, int rank);

	//@override Write()
	virtual ssize_t Write(const Data_3D& data);
	virtual ssize_t Write(const Data_3D& data, bool formated);

	//@override Read()
	virtual ssize_t Read(Data_3D& data);
	virtual ssize_t Read(Data_3D& data, bool formated);

	//@override Lseek()
	virtual off_t   Lseek(off_t off);

	//@override Open()
	virtual int     Open(const std::string& filename);

	//@override Close();
	virtual int     Close();

	virtual void SetView(int block_id, const hsize_t dimsf[3], const hsize_t chunk_dims[3], const std::string& name);

	virtual void SetDataspace(const hsize_t dimsf[3], const hsize_t chunk_dims[3]);
	virtual void SetDatasetid(const hsize_t chunk_dims[3], const std::string& dataname);
	virtual void SetChunk(int blockid, const hsize_t chunk_dims[3], const hsize_t chunk_count[3]);

protected:
	MPI_Comm comm_;
	hid_t fileid_;
	std::string filename_;
	int rank_;

	hid_t datasetid_;
	hid_t dataspace_;
	hid_t chunkspace_;
};

//for OneFilePerProcessAllWrite
class HDF5IOStrategyA: public HDF5IOStrategy {
public:
  HDF5IOStrategyA(const MPI_Comm& comm, int rank);

  //@override Write()
  virtual ssize_t Write(const Data_3D& data);
  virtual ssize_t Write(const Data_3D& data, bool formated);

	//@override Read()
	virtual ssize_t Read(Data_3D& data);
	virtual ssize_t Read(Data_3D& data, bool formated);

	//@override Lseek()
	virtual off_t   Lseek(off_t off);

	//@override Open()
	virtual int     Open(const std::string& filename);

	//@override Close();
	virtual int     Close();

private:
//	hid_t datasetid_;
//	hid_t dataspace_;
//	hid_t chunkspace_;
};

//for SingleSharedFileAllWrite
class HDF5IOStrategyB: public HDF5IOStrategy {
  public:
    HDF5IOStrategyB(const MPI_Comm& comm, int rank);

    //@override Write()
    virtual ssize_t Write(const Data_3D& data);
    virtual ssize_t Write(const Data_3D& data, bool formated);

	//@override Read()
	virtual ssize_t Read(Data_3D& data);
	virtual ssize_t Read(Data_3D& data, bool formated);

	//@override Lseek()
	virtual off_t   Lseek(off_t off);

	//@override Open()
	virtual int     Open(const std::string& filename);

	//@override Close();
	virtual int     Close();

//	virtual void SetDimStride(int nx, int ny, int nz);
//	virtual void SetDataspace(const hsize_t dimsf[3], const hsize_t chunk_dims[3]);
//	virtual void SetDatasetid(const hsize_t chunk_dims[3], const std::string& dataname);
//	virtual void SetChunk(int blockid, const hsize_t chunk_dims[3], const hsize_t chunk_count[3]);

  private:
//	int nx_;
//	int ny_;
//	int nz_;
//	hid_t datasetid_;
//	hid_t dataspace_;
//	hid_t chunkspace_;
};

//for SingleSharedFileOneWrites
class HDF5IOStrategyC: public HDF5IOStrategy {
public:
  HDF5IOStrategyC(const MPI_Comm& comm, int rank);

  //@override Write()
  virtual ssize_t Write(const Data_3D& data);
  virtual ssize_t Write(const Data_3D& data, bool formated);
  ssize_t WriteVector(const std::vector<Data_3D>& data);
  void DatasetSetChunks(int nCols, int size, int maxRows, hid_t plistDCreate);

	//@override Read()
	virtual ssize_t Read(Data_3D& data);
	virtual ssize_t Read(Data_3D& data, bool formated);

	//@override Lseek()
	virtual off_t   Lseek(off_t off);

	//@override Open()
	virtual int     Open(const std::string& filename);

	//@override Close();
	virtual int     Close();
};

//for SingleSharedFileSubsetWrites
class HDF5IOStrategyD: public HDF5IOStrategy {
public:
  HDF5IOStrategyD(const MPI_Comm& comm, int rank);

  //@override Write()
  virtual ssize_t Write(const Data_3D& data);
  virtual ssize_t Write(const Data_3D& data, bool formated);

	//@override Read()
	virtual ssize_t Read(Data_3D& data);
	virtual ssize_t Read(Data_3D& data, bool formated);

	//@override Lseek()
	virtual off_t   Lseek(off_t off);

	//@override Open()
	virtual int     Open(const std::string& filename);

	//@override Close();
	virtual int     Close();
};
#endif
