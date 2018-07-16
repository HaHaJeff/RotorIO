#include "IOStrategy/MPIIO.h"

MPIIO::MPIIO(const MPI_Comm& comm, int rank) : comm_(comm), rank_(rank){

}

MPIIOStrategy* MPIIO::GetIOStrategy(TYPE type) const {
	MPIIOStrategy* pRet = NULL;
	switch(type) {
		//for OneFilePerProcessAllWrite
		case 0:
			pRet = new MPIIOStrategyA(comm_, rank_);
			break;
		//for SingleSharedFileOnelWrites
		case 1:
			pRet = new MPIIOStrategyB(comm_, rank_);
			break;
		//for SingleSharedFileAlllWrite
		case 2:
			pRet = new MPIIOStrategyC(comm_, rank_);
			break;
		//for SingleSharedFileSubsetWrite
		case 3:
			pRet = new MPIIOStrategyD(comm_, rank_);
			break;
		default:
			break;
	}

	return pRet;
}
