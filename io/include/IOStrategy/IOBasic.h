#ifndef IOBASIC_H
#define IOBASIC_H

#include "Strategy.h"
#include "Data.h"

class IOBasic{
public:
	//OneFilePerProcessAllWrite
	//SingleSharedFileOneWrites
	//SingleSharedFileAllWrite
	//SingleSharedFileSubsetWrites
	virtual Strategy* GetIOStrategy(TYPE type) const = 0;
	
//	ssize_t Write(const Data_3D& data);

//	ssize_t Read(Data_3D& data);
};

#endif
