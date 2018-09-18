#ifndef STRING_H
#define STRING_H
#include <string>

/*
struct Data_3D {
	Data_3D(double *pData, int nx, int ny, int nz, std::string name, int blockid = 0): nx_(nx), ny_(ny),nz_(nz), name_(name), blockid_(blockid), pData_(pData) {}
	std::string name_;
	int GetCount() const {return nx_*ny_*nz_;}
	double *pData_;
	int nx_;
	int ny_;
	int nz_;
	int blockid_;
};
*/

/*
#define Data_3D(T) T##Data_3D

#define _Define(T) \
class T##Data_3D \
{\
	T##Data_3D(T* pData, int count,int block_id, const std::string& name) {}; \
	std::string name_;\
	T* pData_;\
	int count_;\
	int block_id;\
	int GetCount()const {return count_;};\
\
};

_Define(double)
_Define(int)
*/


struct Data_3D {
	explicit Data_3D(double* pData, int count, int block_id, const std::string &name) : \
		pData_(pData), count_(count), blockid_(block_id), name_(name) {}
	
	double* pData_;
	int count_;
	int blockid_;
	std::string name_;

	int GetCount() const {return count_;}
};


#endif
