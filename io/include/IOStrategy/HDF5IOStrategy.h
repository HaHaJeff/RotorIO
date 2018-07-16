#ifndef HDF5IOSTRATEGY_H
#define HDF5IOSTRATEGY_H

#include "hdf5.h"
#include "Strategy.h"

#include <memory>

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

protected:
	MPI_Comm comm_;
	hid_t fileid_;
	std::string filename_;
	int rank_;
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
};

//for SingleSharedFileOneWrites
class HDF5IOStrategyC: public HDF5IOStrategy {
public:
  HDF5IOStrategyC(const MPI_Comm& comm, int rank);

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
