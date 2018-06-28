#ifndef STRATEGY_H
#define STRATEGY_H

#include <string>
#include "Data.h"

enum TYPE {OneFilePerProcessAllWrite, SingleSharedFileOneWrites, SingleSharedFileAllWrite, singleSharedFileSubsetWrite};

//Concrect operation on file
class Strategy {
public:
	virtual ssize_t Write(const Data_3D& data) = 0;
	virtual ssize_t Read(Data_3D& data)  = 0;
	virtual off_t   Lseek(off_t off) = 0;
	virtual int     Open(const std::string& filename)  = 0;
	virtual int     Close() = 0;
};

#endif
