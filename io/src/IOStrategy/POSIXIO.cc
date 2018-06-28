#include "POSIXIO.h"
#include "POSIXIOStrategy.h"

POSIXIO::POSIXIO() {

}

POSIXIOStrategy* POSIXIO::GetIOStrategy(TYPE type) const {
	POSIXIOStrategy *pRet = NULL;

	switch(type) {
		//for OneFilePerProcessAllWrite
		case 0:
			pRet = new POSIXIOStrategyOneFilePerProcessAllWrite();
			break;
		//for SingleSharedFileOnelWrites
		case 1:
			pRet = new POSIXIOStrategySingleSharedFileOneWrites();
			break;
		//for SingleSharedFileAlllWrite
		case 2:
			pRet = new POSIXIOStrategySingleSharedFileAllWrite();
			break;
		//for SingleSharedFileSubsetWrite
		case 3:
			pRet = new POSIXIOStrategySingleSharedFileSubsetWrite();
			break;
		default:
			break;
	}

	return pRet;
}
