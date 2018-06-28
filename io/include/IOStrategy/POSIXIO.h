#ifndef POSIXIO_H
#define POSIXIO_H

#include "IOBasic.h"
#include "POSIXIOStrategy.h"

class POSIXIO : public IOBasic {
public:
	POSIXIO();
	virtual POSIXIOStrategy* GetIOStrategy(TYPE type) const;
	
private:
	POSIXIOStrategy* stategy_;
};

#endif
