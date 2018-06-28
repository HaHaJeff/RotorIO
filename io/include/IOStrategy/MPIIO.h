#ifndef MPIIO_H
#define MPIIO_H

#include "IOBasic.h"
#include "MPIIOStrategy.h"

//std::shared_ptr<>
#include <memory>

class MPIIO : public IOBasic {
public:
	MPIIO(const MPI_Comm& comm);
	virtual MPIIOStrategy* GetIOStrategy(TYPE type) const;

private:
	MPI_Comm comm_;
};

#endif
