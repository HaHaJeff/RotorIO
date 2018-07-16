#ifndef HDF5IO_H
#define HDF5IO_h

#include "IOBasic.h"
#include "HDF5IOStrategy.h"
#include "hdf5.h"

class HDF5IO: public IOBasic {
public:
	HDF5IO(const MPI_Comm& comm, int rank);
	virtual HDF5IOStrategy* GetIOStrategy(TYPE type) const;

private:
	MPI_Comm comm_;
	int rank_;
};

#endif
