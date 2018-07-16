#include "HDF5IO.h"

HDF5IO::HDF5IO(const MPI_Comm& comm, int rank) : comm_(comm), rank_(rank) {
}

HDF5IOStrategy* HDF5IO::GetIOStrategy(TYPE type) const {
  HDF5IOStrategy* pRet = NULL;
	switch(type) {
		//for OneFilePerProcessAllWrite
		case 0:
			pRet = new HDF5IOStrategyA(comm_, rank_);
			break;
		//for SingleSharedFileOnelWrites
		case 1:
			pRet = new HDF5IOStrategyB(comm_, rank_);
			break;
		//for SingleSharedFileAlllWrite
		case 2:
			pRet = new HDF5IOStrategyC(comm_, rank_);
			break;
		//for SingleSharedFileSubsetWrite
		case 3:
			pRet = new HDF5IOStrategyD(comm_, rank_);
			break;
		default:
			break;
  }
  return pRet;
}
