#ifndef MPIIOSTRATEGY_H
#define MPIIOSTRATEGY_H

#include <mpi.h>

#include "Strategy.h"

class MPIIOStrategy : public Strategy {
public:
	MPIIOStrategy(const MPI_Comm& comm, int rank);

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

	int 			Split(int groupsize);	
	int 			SetView(int offset, MPI_Datatype etype, MPI_Datatype ftype);
	int 			SetInfo(const std::string& key, const std::string& value);

	const MPI_Comm& GetComm() const {return comm_;}

	const std::string& GetFilename() const {return filename_;}

	const MPI_File& GetFh() const {return fh_;}

protected:
	MPI_Info    info_;
	MPI_Comm 	comm_;
	MPI_File 	fh_;
	std::string filename_;
	int         rank_;
};

class MPIIOStrategyOneFilePerProcessAllWrite : public MPIIOStrategy {
public:

	MPIIOStrategyOneFilePerProcessAllWrite(const MPI_Comm& comm, int rank);

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

class MPIIOStrategySingleSharedFileOneWrites : public MPIIOStrategy {
public:
	MPIIOStrategySingleSharedFileOneWrites(const MPI_Comm& comm, int rank);

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

class MPIIOStrategySingleSharedFileAllWrite: public MPIIOStrategy {
public:
	MPIIOStrategySingleSharedFileAllWrite(const MPI_Comm& comm, int rank);

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

class MPIIOStrategySingleSharedFileSubsetWrite: public MPIIOStrategy {
public:
	MPIIOStrategySingleSharedFileSubsetWrite(const MPI_Comm& comm, int rank);

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
