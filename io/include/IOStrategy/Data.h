#ifndef STRING_H
#define STRING_H

#include <string>

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

#endif
