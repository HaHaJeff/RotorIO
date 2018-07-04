/*---------------------------------------------------------------------------*\


\*---------------------------------------------------------------------------*/

#ifndef writefile_H
#define writefile_H
#include "RKBase.h"

#include <sys/stat.h>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/*---------------------------------------------------------------------------*\
		Class WriteFile Declaration
\*---------------------------------------------------------------------------*/

#include <stdio.h>
#include <stdlib.h>
#include <sys/times.h>

struct Data_3D;

class Time {
	public:
		explicit Time() {
			clktck_ = (double)sysconf(_SC_CLK_TCK);
			start_ = times(&t_start_);
		}

		void print() {
			printf("real: %7.2f, user: %7.2f, system: %7.2f\n", (end_ - start_)/(double)clktck_, (t_end_.tms_utime - t_start_.tms_utime)/(double)clktck_, (t_end_.tms_stime - t_start_.tms_stime)/(double)clktck_ );
		}

		~Time() {
			end_ = times(&t_end_);
			print();
		}
	private:
		struct tms t_start_;
		struct tms t_end_;
		clock_t start_;
		clock_t end_;
	    long clktck_;
};

/*
struct Data_3D {
	Data_3D(double ***pData, int nx, int ny, int nz, string name): nx_(nx), ny_(ny),nz_(nz), name_(name), pData_(pData) {}
	int nx_;
	int ny_;
	int nz_;
	string name_;
	double ***pData_;
};

struct Data_2D {
	Data_2D(double **pData, int nx, int ny, string name): nx_(nx), ny_(ny), name_(name), pData_(pData) {}
	int nx_;
	int ny_;
	string name_;
	double **pData_;
};
*/

class CWriteFile:
	virtual public CRKBase
{
public:
	MPI_Status status;
	string id_m;
	double*** wma;
	fstream f,f1, f2;
	string fOut;
private:
	typedef double ***Dim;
	int DIM[3];
public:
	CWriteFile(CMultiGrid* pgrid, CDictionary* pdict, CGeoField* pfield,int myid):CRKBase(pgrid, pdict, pfield, myid) ,fOut("../output/")	
	{
		id_m = Int_to_string(myid);
		string fError = fOut + "error/";
		string file = fError + "error-" + id_m + "myid.dat";
		f.open(file.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
		
		if (myid == 0)
		{
			string fF1 = fOut + "convergence-inflow.dat";
			string fF2 = fOut + "convergence-outflow.dat";

			f1.open(fF1.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
			f2.open(fF2.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
		}

		Malloc(wma,(pfield->mesh)->nzm1 + 1, (pfield->mesh)->nym1 + 1, (pfield->mesh)->nxm1 + 1);
	}
	~CWriteFile()
	{
		Free(wma);
	}
public:

	void output();
	void flow(int nitt);
	inline void Write_itertime(int iteration,double rmsm,double time_end,double time_begin)
	{
		f << setprecision(17) << iteration << ' ' << rmsm << ' ' << time_end - time_begin << endl;
	}

	inline string Int_to_string(int num)
	{
		ostringstream os;
		os<<num;
		return os.str();
	}

private:
	void InitPressure(int nn);
	void RelativeMa(int nn);
	void lbout();
	void inlout();
	double CubicInterpolation(double* medx, double*medr, double spax, int n1, int n2,int n3);
	void span(int spa);

private:

	void ReadData3D(const string& fileName, Data_3D&, const MPI_Comm& comm);
	void WriteData3D(const string& fileName, const Data_3D&, const MPI_Comm& comm);
//	void WriteData2D(const string& fileName, const Data_2D&);
	void WriteHeader(const string& fileName, const string& zoneName, bool isCell);

	string CreateDir(string dirName, const string& root) {
		string dir = root + dirName;
		int mode = 0777;

		if ( access(dir.c_str(), F_OK) == -1 ) {

			if (mkdir(dir.c_str(), mode) != 0) {	
				printf("%s", "mkdir error!");
				exit(-1);
			}
			
		}

		dir += "/";
		return dir;
	}



};

class CellToPoint : virtual public CRKBase{
public:
	typedef double ***Dim;
	explicit CellToPoint(int x, int y, int z) : CRKBase(NULL, NULL, NULL, 0), nx (x), ny(y), nz(z){
		Malloc(pvxn, nz + 2, ny + 2, nx + 2);
		Malloc(pvyn, nz + 2, ny + 2, nx + 2);
		Malloc(pvzn, nz + 2, ny + 2, nx + 2);
		Malloc(wman, nz + 2, ny + 2, nx + 2);
		Malloc(pn, nz + 2, ny + 2, nx + 2);
		Malloc(q11n, nz + 2, ny + 2, nx + 2);
		Malloc(q16n, nz + 2, ny + 2, nx + 2);

	}



	void Convert(const Dim& pvx, const Dim& pvy, const Dim& pvz, const Dim& wma, const Dim& p, const Dim& q11, const Dim& q16) {

		OneEightCenter_(pvx, pvy, pvz, wma, p, q11, q16);
		SixFace_(pvz, pvy, pvz, wma, p, q11, q16);
		TwevleLine_(pvx, pvy, pvz, wma, p, q11, q16);
		EightPoint_(pvx, pvy, pvz, wma, p, q11, q16);

	}
private:
	void OneEightCenter_(const Dim& pvx, const Dim& pvy, const Dim& pvz, const Dim& wma, const Dim& p, const Dim& q11, const Dim& q16) {
		for (int i = 2; i < nz + 1; i++)
			for (int j = 2; j < ny + 1; j++)
				for (int k = 2; k < nx + 1; k++) {
					pvxn[i][j][k] = 0.125*(pvx[i - 1][j - 1][k - 1] + pvx[i - 1][j - 1][k] + pvx[i - 1][j][k - 1] + pvx[i][j - 1][k - 1] + pvx[i - 1][j][k] + pvx[i][j - 1][k] + pvx[i][j][k - 1] + pvx[i][j][k]);
					pvyn[i][j][k] = 0.125*(pvy[i - 1][j - 1][k - 1] + pvy[i - 1][j - 1][k] + pvy[i - 1][j][k - 1] + pvy[i][j - 1][k - 1] + pvy[i - 1][j][k] + pvy[i][j - 1][k] + pvy[i][j][k - 1] + pvy[i][j][k]);
					pvzn[i][j][k] = 0.125*(pvz[i - 1][j - 1][k - 1] + pvz[i - 1][j - 1][k] + pvz[i - 1][j][k - 1] + pvz[i][j - 1][k - 1] + pvz[i - 1][j][k] + pvz[i][j - 1][k] + pvz[i][j][k - 1] + pvz[i][j][k]);
					wman[i][j][k] = 0.125*(wma[i - 1][j - 1][k - 1] + wma[i - 1][j - 1][k] + wma[i - 1][j][k - 1] + wma[i][j - 1][k - 1] + wma[i - 1][j][k] + wma[i][j - 1][k] + wma[i][j][k - 1] + wma[i][j][k]);
					pn[i][j][k] = 0.125*(p[i - 1][j - 1][k - 1] + p[i - 1][j - 1][k] + p[i - 1][j][k - 1] + p[i][j - 1][k - 1] + p[i - 1][j][k] + p[i][j - 1][k] + p[i][j][k - 1] + p[i][j][k]);
					q11n[i][j][k] = 0.125*(q11[i - 1][j - 1][k - 1] + q11[i - 1][j - 1][k] + q11[i - 1][j][k - 1] + q11[i][j - 1][k - 1] + q11[i - 1][j][k] + q11[i][j - 1][k] + q11[i][j][k - 1] + q11[i][j][k]);
					q16n[i][j][k] = 0.125*(q16[i - 1][j - 1][k - 1] + q16[i - 1][j - 1][k] + q16[i - 1][j][k - 1] + q16[i][j - 1][k - 1] + q16[i - 1][j][k] + q16[i][j - 1][k] + q16[i][j][k - 1] + q16[i][j][k]);
				}
	}

	void SixFace_(const Dim& pvx, const Dim& pvy, const Dim& pvz, const Dim& wma, const Dim& p, const Dim& q11, const Dim& q16) {
		double temp1, temp2, temp3, temp4, temp5, temp6, temp7;
		for (int i = 1; i < nz; i++)
			for (int j = 1; j < ny; j++) {
				pvxn[i + 1][j + 1][1] = 0.5*(pvx[i][j][1] + pvx[i + 1][j][1] + pvx[i][j + 1][1] + pvx[i + 1][j + 1][1]) - pvxn[i + 1][j + 1][2];
				pvyn[i + 1][j + 1][1] = 0.5*(pvy[i][j][1] + pvy[i + 1][j][1] + pvy[i][j + 1][1] + pvy[i + 1][j + 1][1]) - pvyn[i + 1][j + 1][2];
				pvzn[i + 1][j + 1][1] = 0.5*(pvz[i][j][1] + pvz[i + 1][j][1] + pvz[i][j + 1][1] + pvz[i + 1][j + 1][1]) - pvzn[i + 1][j + 1][2];
				wman[i + 1][j + 1][1] = 0.5*(wma[i][j][1] + wma[i + 1][j][1] + wma[i][j + 1][1] + wma[i + 1][j + 1][1]) - wman[i + 1][j + 1][2];
				pn[i + 1][j + 1][1] = 0.5*(p[i][j][1] + p[i + 1][j][1] + p[i][j + 1][1] + p[i + 1][j + 1][1]) - pn[i + 1][j + 1][2];
				q11n[i + 1][j + 1][1] = 0.5*(q11[i][j][1] + q11[i + 1][j][1] + q11[i][j + 1][1] + q11[i + 1][j + 1][1]) - q11n[i + 1][j + 1][2];
				q16n[i + 1][j + 1][1] = 0.5*(q16[i][j][1] + q16[i + 1][j][1] + q16[i][j + 1][1] + q16[i + 1][j + 1][1]) - q16n[i + 1][j + 1][2];

				pvxn[i + 1][j + 1][nx + 1] = 0.5*(pvx[i][j][nx] + pvx[i + 1][j][nx] + pvx[i][j + 1][nx] + pvx[i + 1][j + 1][nx]) - pvxn[i + 1][j + 1][nx];
				pvyn[i + 1][j + 1][nx + 1] = 0.5*(pvy[i][j][nx] + pvy[i + 1][j][nx] + pvy[i][j + 1][nx] + pvy[i + 1][j + 1][nx]) - pvyn[i + 1][j + 1][nx];
				pvzn[i + 1][j + 1][nx + 1] = 0.5*(pvz[i][j][nx] + pvz[i + 1][j][nx] + pvz[i][j + 1][nx] + pvz[i + 1][j + 1][nx]) - pvzn[i + 1][j + 1][nx];
				wman[i + 1][j + 1][nx + 1] = 0.5*(wma[i][j][nx] + wma[i + 1][j][nx] + wma[i][j + 1][nx] + wma[i + 1][j + 1][nx]) - wman[i + 1][j + 1][nx];
				pn[i + 1][j + 1][nx + 1] = 0.5*(p[i][j][nx] + p[i + 1][j][nx] + p[i][j + 1][nx] + p[i + 1][j + 1][nx]) - pn[i + 1][j + 1][nx];
				q11n[i + 1][j + 1][nx + 1] = 0.5*(q11[i][j][nx] + q11[i + 1][j][nx] + q11[i][j + 1][nx] + q11[i + 1][j + 1][nx]) - q11n[i + 1][j + 1][nx];
				q16n[i + 1][j + 1][nx + 1] = 0.5*(q16[i][j][nx] + q16[i + 1][j][nx] + q16[i][j + 1][nx] + q16[i + 1][j + 1][nx]) - q16n[i + 1][j + 1][nx];

			}
		for (int i = 1; i < nz; i++)
			for (int k = 1; k < nx; k++) {
				pvxn[i + 1][1][k + 1] = 0.5*(pvx[i][1][k] + pvx[i + 1][1][k] + pvx[i][1][k + 1] + pvx[i + 1][1][k + 1]) - pvxn[i + 1][2][k + 1];
				pvyn[i + 1][1][k + 1] = 0.5*(pvy[i][1][k] + pvy[i + 1][1][k] + pvy[i][1][k + 1] + pvy[i + 1][1][k + 1]) - pvyn[i + 1][2][k + 1];
				pvzn[i + 1][1][k + 1] = 0.5*(pvz[i][1][k] + pvz[i + 1][1][k] + pvz[i][1][k + 1] + pvz[i + 1][1][k + 1]) - pvzn[i + 1][2][k + 1];
				wman[i + 1][1][k + 1] = 0.5*(wma[i][1][k] + wma[i + 1][1][k] + wma[i][1][k + 1] + wma[i + 1][1][k + 1]) - wman[i + 1][2][k + 1];
				pn[i + 1][1][k + 1] = 0.5*(p[i][1][k] + p[i + 1][1][k] + p[i][1][k + 1] + p[i + 1][1][k + 1]) - pn[i + 1][2][k + 1];
				q11n[i + 1][1][k + 1] = 0.5*(q11[i][1][k] + q11[i + 1][1][k] + q11[i][1][k + 1] + q11[i + 1][1][k + 1]) - q11n[i + 1][2][k + 1];
				q16n[i + 1][1][k + 1] = 0.5*(q16[i][1][k] + q16[i + 1][1][k] + q16[i][1][k + 1] + q16[i + 1][1][k + 1]) - q16n[i + 1][2][k + 1];

				pvxn[i + 1][ny + 1][k + 1] = 0.5*(pvx[i][ny][k] + pvx[i + 1][ny][k] + pvx[i][ny][k + 1] + pvx[i + 1][ny][k + 1]) - pvxn[i + 1][ny][k + 1];
				pvyn[i + 1][ny + 1][k + 1] = 0.5*(pvy[i][ny][k] + pvy[i + 1][ny][k] + pvy[i][ny][k + 1] + pvy[i + 1][ny][k + 1]) - pvyn[i + 1][ny][k + 1];
				pvzn[i + 1][ny + 1][k + 1] = 0.5*(pvz[i][ny][k] + pvz[i + 1][ny][k] + pvz[i][ny][k + 1] + pvz[i + 1][ny][k + 1]) - pvzn[i + 1][ny][k + 1];
				wman[i + 1][ny + 1][k + 1] = 0.5*(wma[i][ny][k] + wma[i + 1][ny][k] + wma[i][ny][k + 1] + wma[i + 1][ny][k + 1]) - wman[i + 1][ny][k + 1];
				pn[i + 1][ny + 1][k + 1] = 0.5*(p[i][ny][k] + p[i + 1][ny][k] + p[i][ny][k + 1] + p[i + 1][ny][k + 1]) - pn[i + 1][ny][k + 1];
				q11n[i + 1][ny + 1][k + 1] = 0.5*(q11[i][ny][k] + q11[i + 1][ny][k] + q11[i][ny][k + 1] + q11[i + 1][ny][k + 1]) - q11n[i + 1][ny][k + 1];
				q16n[i + 1][ny + 1][k + 1] = 0.5*(q16[i][ny][k] + q16[i + 1][ny][k] + q16[i][ny][k + 1] + q16[i + 1][ny][k + 1]) - q16n[i + 1][ny][k + 1];
			}
		for (int j = 1; j < ny; j++)
			for (int k = 1; k < nx; k++) {
				pvxn[1][j + 1][k + 1] = 0.5*(pvx[1][j][k] + pvx[1][j + 1][k] + pvx[1][j][k + 1] + pvx[1][j + 1][k + 1]) - pvxn[2][j + 1][k + 1];
				pvyn[1][j + 1][k + 1] = 0.5*(pvy[1][j][k] + pvy[1][j + 1][k] + pvy[1][j][k + 1] + pvy[1][j + 1][k + 1]) - pvyn[2][j + 1][k + 1];
				pvzn[1][j + 1][k + 1] = 0.5*(pvz[1][j][k] + pvz[1][j + 1][k] + pvz[1][j][k + 1] + pvz[1][j + 1][k + 1]) - pvzn[2][j + 1][k + 1];
				wman[1][j + 1][k + 1] = 0.5*(wma[1][j][k] + wma[1][j + 1][k] + wma[1][j][k + 1] + wma[1][j + 1][k + 1]) - wman[2][j + 1][k + 1];
				pn[1][j + 1][k + 1] = 0.5*(p[1][j][k] + p[1][j + 1][k] + p[1][j][k + 1] + p[1][j + 1][k + 1]) - pn[2][j + 1][k + 1];
				q11n[1][j + 1][k + 1] = 0.5*(q11[1][j][k] + q11[1][j + 1][k] + q11[1][j][k + 1] + q11[1][j + 1][k + 1]) - q11n[2][j + 1][k + 1];
				q16n[1][j + 1][k + 1] = 0.5*(q16[1][j][k] + q16[1][j + 1][k] + q16[1][j][k + 1] + q16[1][j + 1][k + 1]) - q16n[2][j + 1][k + 1];

				pvxn[nz + 1][j + 1][k + 1] = 0.5*(pvx[nz][j][k] + pvx[nz][j + 1][k] + pvx[nz][j][k + 1] + pvx[nz][j + 1][k + 1]) - pvxn[nz][j + 1][k + 1];
				pvyn[nz + 1][j + 1][k + 1] = 0.5*(pvy[nz][j][k] + pvy[nz][j + 1][k] + pvy[nz][j][k + 1] + pvy[nz][j + 1][k + 1]) - pvyn[nz][j + 1][k + 1];
				pvzn[nz + 1][j + 1][k + 1] = 0.5*(pvz[nz][j][k] + pvz[nz][j + 1][k] + pvz[nz][j][k + 1] + pvz[nz][j + 1][k + 1]) - pvzn[nz][j + 1][k + 1];
				wman[nz + 1][j + 1][k + 1] = 0.5*(wma[nz][j][k] + wma[nz][j + 1][k] + wma[nz][j][k + 1] + wma[nz][j + 1][k + 1]) - wman[nz][j + 1][k + 1];
				pn[nz + 1][j + 1][k + 1] = 0.5*(p[nz][j][k] + p[nz][j + 1][k] + p[nz][j][k + 1] + p[nz][j + 1][k + 1]) - pn[nz][j + 1][k + 1];
				q11n[nz + 1][j + 1][k + 1] = 0.5*(q11[nz][j][k] + q11[nz][j + 1][k] + q11[nz][j][k + 1] + q11[nz][j + 1][k + 1]) - q11n[nz][j + 1][k + 1];
				q16n[nz + 1][j + 1][k + 1] = 0.5*(q16[nz][j][k] + q16[nz][j + 1][k] + q16[nz][j][k + 1] + q16[nz][j + 1][k + 1]) - q16n[nz][j + 1][k + 1];
			}

	}

	void TwevleLine_(const Dim& pvx, const Dim& pvy, const Dim& pvz, const Dim& wma, const Dim& p, const Dim& q11, const Dim& q16) {
		for (int i = 2; i < nz + 1; i++) {
			pvxn[i][1][1] = 2 * (pvx[i - 1][1][1] + pvx[i][1][1]) - pvxn[i][1][2] - pvxn[i][2][1] - pvxn[i][2][2];
			pvxn[i][1][nx + 1] = 2 * (pvx[i - 1][1][nx] + pvx[i][1][nx]) - pvxn[i][1][nx] - pvxn[i][2][nx + 1] - pvxn[i][2][nx];
			pvxn[i][ny + 1][1] = 2 * (pvx[i - 1][ny][1] + pvx[i][ny][1]) - pvxn[i][ny][1] - pvxn[i][ny][2] - pvxn[i][ny + 1][2];
			pvxn[i][ny + 1][nx + 1] = 2 * (pvx[i - 1][ny][nx] + pvx[i][ny][nx]) - pvxn[i][ny + 1][nx] - pvxn[i][ny][nx + 1] - pvxn[i][ny][nx];

			pvyn[i][1][1] = 2 * (pvy[i - 1][1][1] + pvy[i][1][1]) - pvyn[i][1][2] - pvyn[i][2][1] - pvyn[i][2][2];
			pvyn[i][1][nx + 1] = 2 * (pvy[i - 1][1][nx] + pvy[i][1][nx]) - pvyn[i][1][nx] - pvyn[i][2][nx + 1] - pvyn[i][2][nx];
			pvyn[i][ny + 1][1] = 2 * (pvy[i - 1][ny][1] + pvy[i][ny][1]) - pvyn[i][ny][1] - pvyn[i][ny][2] - pvyn[i][ny + 1][2];
			pvyn[i][ny + 1][nx + 1] = 2 * (pvy[i - 1][ny][nx] + pvy[i][ny][nx]) - pvyn[i][ny + 1][nx] - pvyn[i][ny][nx + 1] - pvyn[i][ny][nx];

			pvzn[i][1][1] = 2 * (pvz[i - 1][1][1] + pvz[i][1][1]) - pvzn[i][1][2] - pvzn[i][2][1] - pvzn[i][2][2];
			pvzn[i][1][nx + 1] = 2 * (pvz[i - 1][1][nx] + pvz[i][1][nx]) - pvzn[i][1][nx] - pvzn[i][2][nx + 1] - pvzn[i][2][nx];
			pvzn[i][ny + 1][1] = 2 * (pvz[i - 1][ny][1] + pvz[i][ny][1]) - pvzn[i][ny][1] - pvzn[i][ny][2] - pvzn[i][ny + 1][2];
			pvzn[i][ny + 1][nx + 1] = 2 * (pvz[i - 1][ny][nx] + pvz[i][ny][nx]) - pvzn[i][ny + 1][nx] - pvzn[i][ny][nx + 1] - pvzn[i][ny][nx];

			wman[i][1][1] = 2 * (wma[i - 1][1][1] + wma[i][1][1]) - wman[i][1][2] - wman[i][2][1] - wman[i][2][2];
			wman[i][1][nx + 1] = 2 * (wma[i - 1][1][nx] + wma[i][1][nx]) - wman[i][1][nx] - wman[i][2][nx + 1] - wman[i][2][nx];
			wman[i][ny + 1][1] = 2 * (wma[i - 1][ny][1] + wma[i][ny][1]) - wman[i][ny][1] - wman[i][ny][2] - wman[i][ny + 1][2];
			wman[i][ny + 1][nx + 1] = 2 * (wma[i - 1][ny][nx] + wma[i][ny][nx]) - wman[i][ny + 1][nx] - wman[i][ny][nx + 1] - wman[i][ny][nx];

			pn[i][1][1] = 2 * (p[i - 1][1][1] + p[i][1][1]) - pn[i][1][2] - pn[i][2][1] - pn[i][2][2];
			pn[i][1][nx + 1] = 2 * (p[i - 1][1][nx] + p[i][1][nx]) - pn[i][1][nx] - pn[i][2][nx + 1] - pn[i][2][nx];
			pn[i][ny + 1][1] = 2 * (p[i - 1][ny][1] + p[i][ny][1]) - pn[i][ny][1] - pn[i][ny][2] - pn[i][ny + 1][2];
			pn[i][ny + 1][nx + 1] = 2 * (p[i - 1][ny][nx] + p[i][ny][nx]) - pn[i][ny + 1][nx] - pn[i][ny][nx + 1] - pn[i][ny][nx];

			q11n[i][1][1] = 2 * (q11[i - 1][1][1] + q11[i][1][1]) - q11n[i][1][2] - q11n[i][2][1] - q11n[i][2][2];
			q11n[i][1][nx + 1] = 2 * (q11[i - 1][1][nx] + q11[i][1][nx]) - q11n[i][1][nx] - q11n[i][2][nx + 1] - q11n[i][2][nx];
			q11n[i][ny + 1][1] = 2 * (q11[i - 1][ny][1] + q11[i][ny][1]) - q11n[i][ny][1] - q11n[i][ny][2] - q11n[i][ny + 1][2];
			q11n[i][ny + 1][nx + 1] = 2 * (q11[i - 1][ny][nx] + q11[i][ny][nx]) - q11n[i][ny + 1][nx] - q11n[i][ny][nx + 1] - q11n[i][ny][nx];

			q16n[i][1][1] = 2 * (q16[i - 1][1][1] + q16[i][1][1]) - q16n[i][1][2] - q16n[i][2][1] - q16n[i][2][2];
			q16n[i][1][nx + 1] = 2 * (q16[i - 1][1][nx] + q16[i][1][nx]) - q16n[i][1][nx] - q16n[i][2][nx + 1] - q16n[i][2][nx];
			q16n[i][ny + 1][1] = 2 * (q16[i - 1][ny][1] + q16[i][ny][1]) - q16n[i][ny][1] - q16n[i][ny][2] - q16n[i][ny + 1][2];
			q16n[i][ny + 1][nx + 1] = 2 * (q16[i - 1][ny][nx] + q16[i][ny][nx]) - q16n[i][ny + 1][nx] - q16n[i][ny][nx + 1] - q16n[i][ny][nx];
		}
		for (int j = 2; j < ny + 1; j++) {
			pvxn[1][j][1] = 2 * (pvx[1][j - 1][1] + pvx[1][j][1]) - pvxn[1][j][2] - pvxn[2][j][1] - pvxn[2][j][2];
			pvxn[1][j][nx + 1] = 2 * (pvx[1][j - 1][nx] + pvx[1][j][nx]) - pvxn[1][j][nx] - pvxn[2][j][nx + 1] - pvxn[2][j][nx];
			pvxn[nz + 1][j][1] = 2 * (pvx[nz][j - 1][1] + pvx[nz][j][1]) - pvxn[nz][j][1] - pvxn[nz][j][2] - pvxn[nz + 1][j][2];
			pvxn[nz + 1][j][nx + 1] = 2 * (pvx[nz][j - 1][nx] + pvx[nz][j][nx]) - pvxn[nz][j][nx + 1] - pvxn[nz + 1][j][nx] - pvxn[nz][j][nx];

			pvyn[1][j][1] = 2 * (pvy[1][j - 1][1] + pvy[1][j][1]) - pvyn[1][j][2] - pvyn[2][j][1] - pvyn[2][j][2];
			pvyn[1][j][nx + 1] = 2 * (pvy[1][j - 1][nx] + pvy[1][j][nx]) - pvyn[1][j][nx] - pvyn[2][j][nx + 1] - pvyn[2][j][nx];
			pvyn[nz + 1][j][1] = 2 * (pvy[nz][j - 1][1] + pvy[nz][j][1]) - pvyn[nz][j][1] - pvyn[nz][j][2] - pvyn[nz + 1][j][2];
			pvyn[nz + 1][j][nx + 1] = 2 * (pvy[nz][j - 1][nx] + pvy[nz][j][nx]) - pvyn[nz][j][nx + 1] - pvyn[nz + 1][j][nx] - pvyn[nz][j][nx];

			pvzn[1][j][1] = 2 * (pvz[1][j - 1][1] + pvz[1][j][1]) - pvzn[1][j][2] - pvzn[2][j][1] - pvzn[2][j][2];
			pvzn[1][j][nx + 1] = 2 * (pvz[1][j - 1][nx] + pvz[1][j][nx]) - pvzn[1][j][nx] - pvzn[2][j][nx + 1] - pvzn[2][j][nx];
			pvzn[nz + 1][j][1] = 2 * (pvz[nz][j - 1][1] + pvz[nz][j][1]) - pvzn[nz][j][1] - pvzn[nz][j][2] - pvzn[nz + 1][j][2];
			pvzn[nz + 1][j][nx + 1] = 2 * (pvz[nz][j - 1][nx] + pvz[nz][j][nx]) - pvzn[nz][j][nx + 1] - pvzn[nz + 1][j][nx] - pvzn[nz][j][nx];

			wman[1][j][1] = 2 * (wma[1][j - 1][1] + wma[1][j][1]) - wman[1][j][2] - wman[2][j][1] - wman[2][j][2];
			wman[1][j][nx + 1] = 2 * (wma[1][j - 1][nx] + wma[1][j][nx]) - wman[1][j][nx] - wman[2][j][nx + 1] - wman[2][j][nx];
			wman[nz + 1][j][1] = 2 * (wma[nz][j - 1][1] + wma[nz][j][1]) - wman[nz][j][1] - wman[nz][j][2] - wman[nz + 1][j][2];
			wman[nz + 1][j][nx + 1] = 2 * (wma[nz][j - 1][nx] + wma[nz][j][nx]) - wman[nz][j][nx + 1] - wman[nz + 1][j][nx] - wman[nz][j][nx];

			pn[1][j][1] = 2 * (p[1][j - 1][1] + p[1][j][1]) - pn[1][j][2] - pn[2][j][1] - pn[2][j][2];
			pn[1][j][nx + 1] = 2 * (p[1][j - 1][nx] + p[1][j][nx]) - pn[1][j][nx] - pn[2][j][nx + 1] - pn[2][j][nx];
			pn[nz + 1][j][1] = 2 * (p[nz][j - 1][1] + p[nz][j][1]) - pn[nz][j][1] - pn[nz][j][2] - pn[nz + 1][j][2];
			pn[nz + 1][j][nx + 1] = 2 * (p[nz][j - 1][nx] + p[nz][j][nx]) - pn[nz][j][nx + 1] - pn[nz + 1][j][nx] - pn[nz][j][nx];

			q11n[1][j][1] = 2 * (q11[1][j - 1][1] + q11[1][j][1]) - q11n[1][j][2] - q11n[2][j][1] - q11n[2][j][2];
			q11n[1][j][nx + 1] = 2 * (q11[1][j - 1][nx] + q11[1][j][nx]) - q11n[1][j][nx] - q11n[2][j][nx + 1] - q11n[2][j][nx];
			q11n[nz + 1][j][1] = 2 * (q11[nz][j - 1][1] + q11[nz][j][1]) - q11n[nz][j][1] - q11n[nz][j][2] - q11n[nz + 1][j][2];
			q11n[nz + 1][j][nx + 1] = 2 * (q11[nz][j - 1][nx] + q11[nz][j][nx]) - q11n[nz][j][nx + 1] - q11n[nz + 1][j][nx] - q11n[nz][j][nx];

			q16n[1][j][1] = 2 * (q16[1][j - 1][1] + q16[1][j][1]) - q16n[1][j][2] - q16n[2][j][1] - q16n[2][j][2];
			q16n[1][j][nx + 1] = 2 * (q16[1][j - 1][nx] + q16[1][j][nx]) - q16n[1][j][nx] - q16n[2][j][nx + 1] - q16n[2][j][nx];
			q16n[nz + 1][j][1] = 2 * (q16[nz][j - 1][1] + q16[nz][j][1]) - q16n[nz][j][1] - q16n[nz][j][2] - q16n[nz + 1][j][2];
			q16n[nz + 1][j][nx + 1] = 2 * (q16[nz][j - 1][nx] + q16[nz][j][nx]) - q16n[nz][j][nx + 1] - q16n[nz + 1][j][nx] - q16n[nz][j][nx];
		}
		for (int k = 2; k < nx + 1; k++) {
			pvxn[1][1][k] = 2 * (pvx[1][1][k - 1] + pvx[1][1][k]) - pvxn[1][2][k] - pvxn[2][1][k] - pvxn[2][2][k];
			pvxn[1][ny + 1][k] = 2 * (pvx[1][ny][k - 1] + pvx[1][ny][k]) - pvxn[1][ny][k] - pvxn[2][ny + 1][k] - pvxn[2][ny][k];
			pvxn[nz + 1][1][k] = 2 * (pvx[nz][1][k - 1] + pvx[nz][1][k]) - pvxn[nz][1][k] - pvxn[nz][2][k] - pvxn[nz + 1][2][k];
			pvxn[nz + 1][ny + 1][k] = 2 * (pvx[nz][ny][k - 1] + pvx[nz][ny][k]) - pvxn[nz][ny + 1][k] - pvxn[nz + 1][ny][k] - pvxn[nz][ny][k];

			pvyn[1][1][k] = 2 * (pvy[1][1][k - 1] + pvy[1][1][k]) - pvyn[1][2][k] - pvyn[2][1][k] - pvyn[2][2][k];
			pvyn[1][ny + 1][k] = 2 * (pvy[1][ny][k - 1] + pvy[1][ny][k]) - pvyn[1][ny][k] - pvyn[2][ny + 1][k] - pvyn[2][ny][k];
			pvyn[nz + 1][1][k] = 2 * (pvy[nz][1][k - 1] + pvy[nz][1][k]) - pvyn[nz][1][k] - pvyn[nz][2][k] - pvyn[nz + 1][2][k];
			pvyn[nz + 1][ny + 1][k] = 2 * (pvy[nz][ny][k - 1] + pvy[nz][ny][k]) - pvyn[nz][ny + 1][k] - pvyn[nz + 1][ny][k] - pvyn[nz][ny][k];

			pvzn[1][1][k] = 2 * (pvz[1][1][k - 1] + pvz[1][1][k]) - pvzn[1][2][k] - pvzn[2][1][k] - pvzn[2][2][k];
			pvzn[1][ny + 1][k] = 2 * (pvz[1][ny][k - 1] + pvz[1][ny][k]) - pvzn[1][ny][k] - pvzn[2][ny + 1][k] - pvzn[2][ny][k];
			pvzn[nz + 1][1][k] = 2 * (pvz[nz][1][k - 1] + pvz[nz][1][k]) - pvzn[nz][1][k] - pvzn[nz][2][k] - pvzn[nz + 1][2][k];
			pvzn[nz + 1][ny + 1][k] = 2 * (pvz[nz][ny][k - 1] + pvz[nz][ny][k]) - pvzn[nz][ny + 1][k] - pvzn[nz + 1][ny][k] - pvzn[nz][ny][k];

			wman[1][1][k] = 2 * (wma[1][1][k - 1] + wma[1][1][k]) - wman[1][2][k] - wman[2][1][k] - wman[2][2][k];
			wman[1][ny + 1][k] = 2 * (wma[1][ny][k - 1] + wma[1][ny][k]) - wman[1][ny][k] - wman[2][ny + 1][k] - wman[2][ny][k];
			wman[nz + 1][1][k] = 2 * (wma[nz][1][k - 1] + wma[nz][1][k]) - wman[nz][1][k] - wman[nz][2][k] - wman[nz + 1][2][k];
			wman[nz + 1][ny + 1][k] = 2 * (wma[nz][ny][k - 1] + wma[nz][ny][k]) - wman[nz][ny + 1][k] - wman[nz + 1][ny][k] - wman[nz][ny][k];

			pn[1][1][k] = 2 * (p[1][1][k - 1] + p[1][1][k]) - pn[1][2][k] - pn[2][1][k] - pn[2][2][k];
			pn[1][ny + 1][k] = 2 * (p[1][ny][k - 1] + p[1][ny][k]) - pn[1][ny][k] - pn[2][ny + 1][k] - pn[2][ny][k];
			pn[nz + 1][1][k] = 2 * (p[nz][1][k - 1] + p[nz][1][k]) - pn[nz][1][k] - pn[nz][2][k] - pn[nz + 1][2][k];
			pn[nz + 1][ny + 1][k] = 2 * (p[nz][ny][k - 1] + p[nz][ny][k]) - pn[nz][ny + 1][k] - pn[nz + 1][ny][k] - pn[nz][ny][k];

			q11n[1][1][k] = 2 * (q11[1][1][k - 1] + q11[1][1][k]) - q11n[1][2][k] - q11n[2][1][k] - q11n[2][2][k];
			q11n[1][ny + 1][k] = 2 * (q11[1][ny][k - 1] + q11[1][ny][k]) - q11n[1][ny][k] - q11n[2][ny + 1][k] - q11n[2][ny][k];
			q11n[nz + 1][1][k] = 2 * (q11[nz][1][k - 1] + q11[nz][1][k]) - q11n[nz][1][k] - q11n[nz][2][k] - q11n[nz + 1][2][k];
			q11n[nz + 1][ny + 1][k] = 2 * (q11[nz][ny][k - 1] + q11[nz][ny][k]) - q11n[nz][ny + 1][k] - q11n[nz + 1][ny][k] - q11n[nz][ny][k];

			q16n[1][1][k] = 2 * (q16[1][1][k - 1] + q16[1][1][k]) - q16n[1][2][k] - q16n[2][1][k] - q16n[2][2][k];
			q16n[1][ny + 1][k] = 2 * (q16[1][ny][k - 1] + q16[1][ny][k]) - q16n[1][ny][k] - q16n[2][ny + 1][k] - q16n[2][ny][k];
			q16n[nz + 1][1][k] = 2 * (q16[nz][1][k - 1] + q16[nz][1][k]) - q16n[nz][1][k] - q16n[nz][2][k] - q16n[nz + 1][2][k];
			q16n[nz + 1][ny + 1][k] = 2 * (q16[nz][ny][k - 1] + q16[nz][ny][k]) - q16n[nz][ny + 1][k] - q16n[nz + 1][ny][k] - q16n[nz][ny][k];
		}
	}

	void EightPoint_(const Dim& pvx, const Dim& pvy, const Dim& pvz, const Dim& wma, const Dim& p, const Dim& q11, const Dim& q16) {
		pvxn[1][1][1] = 8 * pvx[1][1][1] - pvxn[1][1][2] - pvxn[1][2][1] - pvxn[2][1][1] - pvxn[1][2][2] - pvxn[2][1][2] - pvxn[2][2][1] - pvxn[2][2][2];
		pvxn[1][1][nx + 1] = 8 * pvx[1][1][nx] - pvxn[1][1][nx] - pvxn[1][2][nx + 1] - pvxn[2][1][nx + 1] - pvxn[1][2][nx] - pvxn[2][1][nx] - pvxn[2][2][nx + 1] - pvxn[2][2][nx];
		pvxn[1][ny + 1][1] = 8 * pvx[1][ny][1] - pvxn[1][ny + 1][2] - pvxn[1][ny][1] - pvxn[2][ny + 1][1] - pvxn[1][ny][2] - pvxn[2][ny + 1][2] - pvxn[2][ny][1] - pvxn[2][ny][2];
		pvxn[nz + 1][1][1] = 8 * pvx[nz][1][1] - pvxn[nz + 1][1][2] - pvxn[nz + 1][2][1] - pvxn[nz][1][1] - pvxn[nz + 1][2][2] - pvxn[nz][1][2] - pvxn[nz][2][1] - pvxn[nz][2][2];
		pvxn[1][ny + 1][nx + 1] = 8 * pvx[1][ny][nx] - pvxn[1][ny + 1][nx] - pvxn[1][ny][nx + 1] - pvxn[2][ny + 1][nx + 1] - pvxn[1][ny][nx] - pvxn[2][ny + 1][nx] - pvxn[2][ny][nx + 1] - pvxn[2][ny][nx];
		pvxn[nz + 1][1][nx + 1] = 8 * pvx[nz][1][nx] - pvxn[nz + 1][1][nx] - pvxn[nz + 1][2][nx + 1] - pvxn[nz][1][nx + 1] - pvxn[nz + 1][2][nx] - pvxn[nz][1][nx] - pvxn[nz][2][nx + 1] - pvxn[nz][2][nx];
		pvxn[nz + 1][ny + 1][1] = 8 * pvx[nz][ny][1] - pvxn[nz + 1][nz + 1][2] - pvxn[nz + 1][ny][1] - pvxn[nz][ny + 1][1] - pvxn[nz + 1][ny][2] - pvxn[nz][ny + 1][2] - pvxn[nz][ny][1] - pvxn[nz][ny][2];
		pvxn[nz + 1][ny + 1][nx + 1] = 8 * pvx[nz][ny][nx] - pvxn[nz + 1][ny + 1][nx] - pvxn[nz + 1][ny][nx + 1] - pvxn[nz][ny + 1][nx + 1] - pvxn[nz + 1][ny][nx] - pvxn[nz][ny + 1][nx] - pvxn[nz][ny][nx + 1] - pvxn[nz][ny][nx];

		pvyn[1][1][1] = 8 * pvy[1][1][1] - pvyn[1][1][2] - pvyn[1][2][1] - pvyn[2][1][1] - pvyn[1][2][2] - pvyn[2][1][2] - pvyn[2][2][1] - pvyn[2][2][2];
		pvyn[1][1][nx + 1] = 8 * pvy[1][1][nx] - pvyn[1][1][nx] - pvyn[1][2][nx + 1] - pvyn[2][1][nx + 1] - pvyn[1][2][nx] - pvyn[2][1][nx] - pvyn[2][2][nx + 1] - pvyn[2][2][nx];
		pvyn[1][ny + 1][1] = 8 * pvy[1][ny][1] - pvyn[1][ny + 1][2] - pvyn[1][ny][1] - pvyn[2][ny + 1][1] - pvyn[1][ny][2] - pvyn[2][ny + 1][2] - pvyn[2][ny][1] - pvyn[2][ny][2];
		pvyn[nz + 1][1][1] = 8 * pvy[nz][1][1] - pvyn[nz + 1][1][2] - pvyn[nz + 1][2][1] - pvyn[nz][1][1] - pvyn[nz + 1][2][2] - pvyn[nz][1][2] - pvyn[nz][2][1] - pvyn[nz][2][2];
		pvyn[1][ny + 1][nx + 1] = 8 * pvy[1][ny][nx] - pvyn[1][ny + 1][nx] - pvyn[1][ny][nx + 1] - pvyn[2][ny + 1][nx + 1] - pvyn[1][ny][nx] - pvyn[2][ny + 1][nx] - pvyn[2][ny][nx + 1] - pvyn[2][ny][nx];
		pvyn[nz + 1][1][nx + 1] = 8 * pvy[nz][1][nx] - pvyn[nz + 1][1][nx] - pvyn[nz + 1][2][nx + 1] - pvyn[nz][1][nx + 1] - pvyn[nz + 1][2][nx] - pvyn[nz][1][nx] - pvyn[nz][2][nx + 1] - pvyn[nz][2][nx];
		pvyn[nz + 1][ny + 1][1] = 8 * pvy[nz][ny][1] - pvyn[nz + 1][nz + 1][2] - pvyn[nz + 1][ny][1] - pvyn[nz][ny + 1][1] - pvyn[nz + 1][ny][2] - pvyn[nz][ny + 1][2] - pvyn[nz][ny][1] - pvyn[nz][ny][2];
		pvyn[nz + 1][ny + 1][nx + 1] = 8 * pvy[nz][ny][nx] - pvyn[nz + 1][ny + 1][nx] - pvyn[nz + 1][ny][nx + 1] - pvyn[nz][ny + 1][nx + 1] - pvyn[nz + 1][ny][nx] - pvyn[nz][ny + 1][nx] - pvyn[nz][ny][nx + 1] - pvyn[nz][ny][nx];

		pvzn[1][1][1] = 8 * pvz[1][1][1] - pvzn[1][1][2] - pvzn[1][2][1] - pvzn[2][1][1] - pvzn[1][2][2] - pvzn[2][1][2] - pvzn[2][2][1] - pvzn[2][2][2];
		pvzn[1][1][nx + 1] = 8 * pvz[1][1][nx] - pvzn[1][1][nx] - pvzn[1][2][nx + 1] - pvzn[2][1][nx + 1] - pvzn[1][2][nx] - pvzn[2][1][nx] - pvzn[2][2][nx + 1] - pvzn[2][2][nx];
		pvzn[1][ny + 1][1] = 8 * pvz[1][ny][1] - pvzn[1][ny + 1][2] - pvzn[1][ny][1] - pvzn[2][ny + 1][1] - pvzn[1][ny][2] - pvzn[2][ny + 1][2] - pvzn[2][ny][1] - pvzn[2][ny][2];
		pvzn[nz + 1][1][1] = 8 * pvz[nz][1][1] - pvzn[nz + 1][1][2] - pvzn[nz + 1][2][1] - pvzn[nz][1][1] - pvzn[nz + 1][2][2] - pvzn[nz][1][2] - pvzn[nz][2][1] - pvzn[nz][2][2];
		pvzn[1][ny + 1][nx + 1] = 8 * pvz[1][ny][nx] - pvzn[1][ny + 1][nx] - pvzn[1][ny][nx + 1] - pvzn[2][ny + 1][nx + 1] - pvzn[1][ny][nx] - pvzn[2][ny + 1][nx] - pvzn[2][ny][nx + 1] - pvzn[2][ny][nx];
		pvzn[nz + 1][1][nx + 1] = 8 * pvz[nz][1][nx] - pvzn[nz + 1][1][nx] - pvzn[nz + 1][2][nx + 1] - pvzn[nz][1][nx + 1] - pvzn[nz + 1][2][nx] - pvzn[nz][1][nx] - pvzn[nz][2][nx + 1] - pvzn[nz][2][nx];
		pvzn[nz + 1][ny + 1][1] = 8 * pvz[nz][ny][1] - pvzn[nz + 1][nz + 1][2] - pvzn[nz + 1][ny][1] - pvzn[nz][ny + 1][1] - pvzn[nz + 1][ny][2] - pvzn[nz][ny + 1][2] - pvzn[nz][ny][1] - pvzn[nz][ny][2];
		pvzn[nz + 1][ny + 1][nx + 1] = 8 * pvz[nz][ny][nx] - pvzn[nz + 1][ny + 1][nx] - pvzn[nz + 1][ny][nx + 1] - pvzn[nz][ny + 1][nx + 1] - pvzn[nz + 1][ny][nx] - pvzn[nz][ny + 1][nx] - pvzn[nz][ny][nx + 1] - pvzn[nz][ny][nx];

		wman[1][1][1] = 8 * wma[1][1][1] - wman[1][1][2] - wman[1][2][1] - wman[2][1][1] - wman[1][2][2] - wman[2][1][2] - wman[2][2][1] - wman[2][2][2];
		wman[1][1][nx + 1] = 8 * wma[1][1][nx] - wman[1][1][nx] - wman[1][2][nx + 1] - wman[2][1][nx + 1] - wman[1][2][nx] - wman[2][1][nx] - wman[2][2][nx + 1] - wman[2][2][nx];
		wman[1][ny + 1][1] = 8 * wma[1][ny][1] - wman[1][ny + 1][2] - wman[1][ny][1] - wman[2][ny + 1][1] - wman[1][ny][2] - wman[2][ny + 1][2] - wman[2][ny][1] - wman[2][ny][2];
		wman[nz + 1][1][1] = 8 * wma[nz][1][1] - wman[nz + 1][1][2] - wman[nz + 1][2][1] - wman[nz][1][1] - wman[nz + 1][2][2] - wman[nz][1][2] - wman[nz][2][1] - wman[nz][2][2];
		wman[1][ny + 1][nx + 1] = 8 * wma[1][ny][nx] - wman[1][ny + 1][nx] - wman[1][ny][nx + 1] - wman[2][ny + 1][nx + 1] - wman[1][ny][nx] - wman[2][ny + 1][nx] - wman[2][ny][nx + 1] - wman[2][ny][nx];
		wman[nz + 1][1][nx + 1] = 8 * wma[nz][1][nx] - wman[nz + 1][1][nx] - wman[nz + 1][2][nx + 1] - wman[nz][1][nx + 1] - wman[nz + 1][2][nx] - wman[nz][1][nx] - wman[nz][2][nx + 1] - wman[nz][2][nx];
		wman[nz + 1][ny + 1][1] = 8 * wma[nz][ny][1] - wman[nz + 1][nz + 1][2] - wman[nz + 1][ny][1] - wman[nz][ny + 1][1] - wman[nz + 1][ny][2] - wman[nz][ny + 1][2] - wman[nz][ny][1] - wman[nz][ny][2];
		wman[nz + 1][ny + 1][nx + 1] = 8 * wma[nz][ny][nx] - wman[nz + 1][ny + 1][nx] - wman[nz + 1][ny][nx + 1] - wman[nz][ny + 1][nx + 1] - wman[nz + 1][ny][nx] - wman[nz][ny + 1][nx] - wman[nz][ny][nx + 1] - wman[nz][ny][nx];

		pn[1][1][1] = 8 * p[1][1][1] - pn[1][1][2] - pn[1][2][1] - pn[2][1][1] - pn[1][2][2] - pn[2][1][2] - pn[2][2][1] - pn[2][2][2];
		pn[1][1][nx + 1] = 8 * p[1][1][nx] - pn[1][1][nx] - pn[1][2][nx + 1] - pn[2][1][nx + 1] - pn[1][2][nx] - pn[2][1][nx] - pn[2][2][nx + 1] - pn[2][2][nx];
		pn[1][ny + 1][1] = 8 * p[1][ny][1] - pn[1][ny + 1][2] - pn[1][ny][1] - pn[2][ny + 1][1] - pn[1][ny][2] - pn[2][ny + 1][2] - pn[2][ny][1] - pn[2][ny][2];
		pn[nz + 1][1][1] = 8 * p[nz][1][1] - pn[nz + 1][1][2] - pn[nz + 1][2][1] - pn[nz][1][1] - pn[nz + 1][2][2] - pn[nz][1][2] - pn[nz][2][1] - pn[nz][2][2];
		pn[1][ny + 1][nx + 1] = 8 * p[1][ny][nx] - pn[1][ny + 1][nx] - pn[1][ny][nx + 1] - pn[2][ny + 1][nx + 1] - pn[1][ny][nx] - pn[2][ny + 1][nx] - pn[2][ny][nx + 1] - pn[2][ny][nx];
		pn[nz + 1][1][nx + 1] = 8 * p[nz][1][nx] - pn[nz + 1][1][nx] - pn[nz + 1][2][nx + 1] - pn[nz][1][nx + 1] - pn[nz + 1][2][nx] - pn[nz][1][nx] - pn[nz][2][nx + 1] - pn[nz][2][nx];
		pn[nz + 1][ny + 1][1] = 8 * p[nz][ny][1] - pn[nz + 1][nz + 1][2] - pn[nz + 1][ny][1] - pn[nz][ny + 1][1] - pn[nz + 1][ny][2] - pn[nz][ny + 1][2] - pn[nz][ny][1] - pn[nz][ny][2];
		pn[nz + 1][ny + 1][nx + 1] = 8 * p[nz][ny][nx] - pn[nz + 1][ny + 1][nx] - pn[nz + 1][ny][nx + 1] - pn[nz][ny + 1][nx + 1] - pn[nz + 1][ny][nx] - pn[nz][ny + 1][nx] - pn[nz][ny][nx + 1] - pn[nz][ny][nx];

		q11n[1][1][1] = 8 * q11[1][1][1] - q11n[1][1][2] - q11n[1][2][1] - q11n[2][1][1] - q11n[1][2][2] - q11n[2][1][2] - q11n[2][2][1] - q11n[2][2][2];
		q11n[1][1][nx + 1] = 8 * q11[1][1][nx] - q11n[1][1][nx] - q11n[1][2][nx + 1] - q11n[2][1][nx + 1] - q11n[1][2][nx] - q11n[2][1][nx] - q11n[2][2][nx + 1] - q11n[2][2][nx];
		q11n[1][ny + 1][1] = 8 * q11[1][ny][1] - q11n[1][ny + 1][2] - q11n[1][ny][1] - q11n[2][ny + 1][1] - q11n[1][ny][2] - q11n[2][ny + 1][2] - q11n[2][ny][1] - q11n[2][ny][2];
		q11n[nz + 1][1][1] = 8 * q11[nz][1][1] - q11n[nz + 1][1][2] - q11n[nz + 1][2][1] - q11n[nz][1][1] - q11n[nz + 1][2][2] - q11n[nz][1][2] - q11n[nz][2][1] - q11n[nz][2][2];
		q11n[1][ny + 1][nx + 1] = 8 * q11[1][ny][nx] - q11n[1][ny + 1][nx] - q11n[1][ny][nx + 1] - q11n[2][ny + 1][nx + 1] - q11n[1][ny][nx] - q11n[2][ny + 1][nx] - q11n[2][ny][nx + 1] - q11n[2][ny][nx];
		q11n[nz + 1][1][nx + 1] = 8 * q11[nz][1][nx] - q11n[nz + 1][1][nx] - q11n[nz + 1][2][nx + 1] - q11n[nz][1][nx + 1] - q11n[nz + 1][2][nx] - q11n[nz][1][nx] - q11n[nz][2][nx + 1] - q11n[nz][2][nx];
		q11n[nz + 1][ny + 1][1] = 8 * q11[nz][ny][1] - q11n[nz + 1][nz + 1][2] - q11n[nz + 1][ny][1] - q11n[nz][ny + 1][1] - q11n[nz + 1][ny][2] - q11n[nz][ny + 1][2] - q11n[nz][ny][1] - q11n[nz][ny][2];
		q11n[nz + 1][ny + 1][nx + 1] = 8 * q11[nz][ny][nx] - q11n[nz + 1][ny + 1][nx] - q11n[nz + 1][ny][nx + 1] - q11n[nz][ny + 1][nx + 1] - q11n[nz + 1][ny][nx] - q11n[nz][ny + 1][nx] - q11n[nz][ny][nx + 1] - q11n[nz][ny][nx];

		q16n[1][1][1] = 8 * q16[1][1][1] - q16n[1][1][2] - q16n[1][2][1] - q16n[2][1][1] - q16n[1][2][2] - q16n[2][1][2] - q16n[2][2][1] - q16n[2][2][2];
		q16n[1][1][nx + 1] = 8 * q16[1][1][nx] - q16n[1][1][nx] - q16n[1][2][nx + 1] - q16n[2][1][nx + 1] - q16n[1][2][nx] - q16n[2][1][nx] - q16n[2][2][nx + 1] - q16n[2][2][nx];
		q16n[1][ny + 1][1] = 8 * q16[1][ny][1] - q16n[1][ny + 1][2] - q16n[1][ny][1] - q16n[2][ny + 1][1] - q16n[1][ny][2] - q16n[2][ny + 1][2] - q16n[2][ny][1] - q16n[2][ny][2];
		q16n[nz + 1][1][1] = 8 * q16[nz][1][1] - q16n[nz + 1][1][2] - q16n[nz + 1][2][1] - q16n[nz][1][1] - q16n[nz + 1][2][2] - q16n[nz][1][2] - q16n[nz][2][1] - q16n[nz][2][2];
		q16n[1][ny + 1][nx + 1] = 8 * q16[1][ny][nx] - q16n[1][ny + 1][nx] - q16n[1][ny][nx + 1] - q16n[2][ny + 1][nx + 1] - q16n[1][ny][nx] - q16n[2][ny + 1][nx] - q16n[2][ny][nx + 1] - q16n[2][ny][nx];
		q16n[nz + 1][1][nx + 1] = 8 * q16[nz][1][nx] - q16n[nz + 1][1][nx] - q16n[nz + 1][2][nx + 1] - q16n[nz][1][nx + 1] - q16n[nz + 1][2][nx] - q16n[nz][1][nx] - q16n[nz][2][nx + 1] - q16n[nz][2][nx];
		q16n[nz + 1][ny + 1][1] = 8 * q16[nz][ny][1] - q16n[nz + 1][nz + 1][2] - q16n[nz + 1][ny][1] - q16n[nz][ny + 1][1] - q16n[nz + 1][ny][2] - q16n[nz][ny + 1][2] - q16n[nz][ny][1] - q16n[nz][ny][2];
		q16n[nz + 1][ny + 1][nx + 1] = 8 * q16[nz][ny][nx] - q16n[nz + 1][ny + 1][nx] - q16n[nz + 1][ny][nx + 1] - q16n[nz][ny + 1][nx + 1] - q16n[nz + 1][ny][nx] - q16n[nz][ny + 1][nx] - q16n[nz][ny][nx + 1] - q16n[nz][ny][nx];
	}
public:
	const Dim& GetPvx() const { return pvxn; }
	const Dim& GetPvy() const { return pvyn; }
	const Dim& GetPvz() const { return pvzn; }
	const Dim& GetWma() const { return wman; }
	const Dim& GetP() const { return pn; }
	const Dim& GetQ11() const { return q11n; }
	const Dim& GetQ16() const { return q16n; }
private:
	Dim pvxn;
	Dim pvyn;
	Dim pvzn;
	Dim wman;
	Dim pn;
	Dim q11n;
	Dim q16n;
	int nx;
	int ny;
	int nz;

};


#endif

// ************************************************************************* //
