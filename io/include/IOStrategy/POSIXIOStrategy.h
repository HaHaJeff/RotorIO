#ifndef POSIXSTRATEGY_H
#define POSIXSTRATEGY_H

#include "Strategy.h"

//TODO: Add formated write and read
class POSIXIOStrategy : public Strategy {
public:
	POSIXIOStrategy();

	//@override Write()
	virtual ssize_t Write(const Data_3D &data);
	virtual ssize_t Write(const Data_3D &data, bool formated);

	//@override Read()
	virtual ssize_t Read(Data_3D &data);
	virtual ssize_t Read(Data_3D &data, bool formated);

	//@override Lseek()
	virtual off_t   Lseek(off_t off);

	//@override Open()
	virtual int     Open(const std::string& filename);

	//@override Close();
	virtual int     Close();

	ssize_t Write(const int* data, int size);
	ssize_t Read(int* data, int size);
protected:
	int fd_;
	std::string filename_;
};

class POSIXIOStrategyOneFilePerProcessAllWrite : public POSIXIOStrategy {
public:
	
	POSIXIOStrategyOneFilePerProcessAllWrite();
	
	//@override Write()
	virtual ssize_t Write(const Data_3D& data);
	virtual ssize_t Write(const Data_3D &data, bool formated);

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

class POSIXIOStrategySingleSharedFileOneWrites : public POSIXIOStrategy {
public:
	
	POSIXIOStrategySingleSharedFileOneWrites();
	
	//@override Write()
	virtual ssize_t Write(const Data_3D &data);
	virtual ssize_t Write(const Data_3D &data, bool formated);

	//@override Read()
	virtual ssize_t Read(Data_3D &data);
	virtual ssize_t Read(Data_3D &data, bool formated);
	
	//@override Lseek()
	virtual off_t   Lseek(off_t off);

	//@override Open()
	virtual int     Open(const std::string& filename);

	//@override Close();
	virtual int     Close();

};

class POSIXIOStrategySingleSharedFileAllWrite : public POSIXIOStrategy {
public:
	
	POSIXIOStrategySingleSharedFileAllWrite();
	
	//@override Write()
	virtual ssize_t Write(const Data_3D &data);
	virtual ssize_t Write(const Data_3D &data, bool formated);

	//@override Read()
	virtual ssize_t Read(Data_3D &data);
	virtual ssize_t Read(Data_3D &data, bool formated);
	
	//@override Lseek()
	virtual off_t   Lseek(off_t off);

	//@override Open()
	virtual int     Open(const std::string& filename);

	//@override Close();
	virtual int     Close();

};

class POSIXIOStrategySingleSharedFileSubsetWrite : public POSIXIOStrategy {
public:
	
	POSIXIOStrategySingleSharedFileSubsetWrite();
	
	//@override Write()
	virtual ssize_t Write(const Data_3D &data);
	virtual ssize_t Write(const Data_3D &data, bool formated);

	//@override Read()
	virtual ssize_t Read(Data_3D &data);
	virtual ssize_t Read(Data_3D &data, bool formated);
	
	//@override Lseek()
	virtual off_t   Lseek(off_t off);

	//@override Open()
	virtual int     Open(const std::string& filename);

	//@override Close();
	virtual int     Close();

};
#endif
