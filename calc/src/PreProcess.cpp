#include <string>
#include <fstream>
#include <sstream>

#include <sys/stat.h>

#include "Dictionary.h"


string CreateDir(const std::string& root, const string& dirName) {
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

template<class T>
bool Malloc(T* &pArr,int size1)
{
	pArr = (T *)malloc(sizeof(T ) * size1 );
 	if(pArr == NULL)
 	{	
		cout<<" memory allocation error"<<endl;
		return -1;	
 	}
	memset(pArr, 0, size1*sizeof(T));
	return true;
}

template<class T>
bool Malloc(T** &pArr,int size1,int size2)
{
	pArr = (T **)malloc(sizeof(T *) *size1); 
	pArr[0] = (T *)malloc(sizeof(T ) * size1 * size2);

 	if(pArr == NULL)
 	{	
		cout<<" memory allocation error"<<endl;
		return -1;	
 	}

	for(int n=1; n<size1; n++)  
	{  
		pArr[n] = pArr[n-1] + size2; 
	}

	
	memset(pArr[0], 0, size1*size2*sizeof(T));

}


template<class T>
bool Malloc(T*** &pArr,int size1,int size2,int size3)
{
	pArr = (T ***)malloc(sizeof(T **) *size1); 
	pArr[0] = (T **)malloc(sizeof(T *) * size1 * size2);
	pArr[0][0]= (T *)malloc(sizeof(T ) * size1 * size2* size3);

 	if(pArr == NULL)
 	{	
		cout<<" memory allocation error"<<endl;
		return -1;	
 	}

	for(int n=1; n<size1; n++)  
	{  
		pArr[n] = pArr[n-1] + size2;
		pArr[n][0]= pArr[n-1][0]+size2*size3; 
	}

	for(int i=0; i<size1; i++) 
	for(int j=1; j<size2; j++) 
	{  
		pArr[i][j] = pArr[i][j-1] + size3;  
	}	
	memset(pArr[0][0], 0, size1*size2*size3*sizeof(T));

	return true;
}

template<class T>
bool Malloc(T**** &pArr,int size1,int size2,int size3,int size4)
{
	pArr = (T ****)malloc(sizeof(T ***) *size1); 
	pArr[0] = (T ***)malloc(sizeof(T **) * size1 * size2);
	pArr[0][0]= (T **)malloc(sizeof(T *) * size1 * size2* size3);
	pArr[0][0][0]= (T *)malloc(sizeof(T ) * size1 * size2* size3*size4);

	for(int n=1; n<size1; n++)  
	{  
		pArr[n] = pArr[n-1] + size2;
		pArr[n][0]= pArr[n-1][0]+size2*size3; 
		pArr[n][0][0]= pArr[n-1][0][0]+size2*size3*size4;
	}

	for(int n1=0; n1<size1; n1++) 
	for(int n2=1; n2<size2; n2++) 
	{  
		pArr[n1][n2] = pArr[n1][n2-1] + size3;
		pArr[n1][n2][0]= pArr[n1][n2-1][0]+size3*size4;  
	}

	for(int i=0; i<size1; i++) 
	for(int j=0; j<size2; j++) 
	for(int k=1; k<size3; k++)
	{  
		pArr[i][j][k] = pArr[i][j][k-1] + size4;  
	}

    if(pArr == 0)
    {
        return -1;
    }
	
    memset(pArr[0][0][0], 0, size1*size2*size3*size4*sizeof(T));
    return true;

}

template<class T>
bool Free(T* &pArr)
{	if(pArr)
	free(pArr);
	pArr=NULL;
	return true;
}

template<class T>
bool Free(T** &pArr)
{	
	if(pArr[0]) 
	free(pArr[0]);
	if(pArr)
	free(pArr);

	pArr=NULL;
	return true;
}

template<class T>
bool Free(T*** &pArr)
{
	if(pArr[0][0])
	free(pArr[0][0]);
	if(pArr[0])
	free(pArr[0]);
	if(pArr)
	free(pArr); 
	

	pArr=NULL;
	return true;	
}

template<class T>
bool Free(T**** &pArr)
{

	if(pArr[0][0][0])
	free(pArr[0][0][0]);	
	if(pArr[0][0])
	free(pArr[0][0]);
	if(pArr[0])
	free(pArr[0]);
	if(pArr)
	free(pArr); 
		
	pArr=NULL;
	return true;
}

static std::string Int_to_string(int num) {
	ostringstream os;
	os << num;
	return os.str();
}

std::string DirToFileName(const std::string& dir, const std::string& file) {	
	std::string ret;
	
	if(dir[dir.size()-1] != '/') {
		ret = dir + "/" + file;
	} else {
		ret = dir + ret;
	}
	return ret;
}

struct Data_3D {
	Data_3D(int nx, int ny, int nz): nx_(nx), ny_(ny),nz_(nz) {
		Malloc(pData_, nx, ny, nz);
	}
	
	Data_3D(double*** pData, int nx, int ny, int nz, const std::string& name) :  nx_(nx), ny_(ny),nz_(nz), pData_(pData), name_(name) {
		
	}
	
	double** operator[](int i) {
		return pData_[i];
	}

	~Data_3D() {
		if (pData_ != NULL) {
			Free(pData_);
		}
	}
	int nx_;
	int ny_;
	int nz_;
	double ***pData_;
	std::string name_;
};

struct Data_2D {
	Data_2D(int nx, int ny): nx_(nx), ny_(ny) {
		Malloc(pData_, nx, ny);
	}
	Data_2D(double** pData, int nx, int ny, const std::string name) :  nx_(nx), ny_(ny), pData_(pData), name_(name) {
		
	}

	double* operator[](int i) {
		return pData_[i];
	}

	~Data_2D() {
		if (!name_.empty() && pData_ != NULL) {
			Free(pData_);
		}
	}
	int nx_;
	int ny_;
	const std::string name_;
	double **pData_;
};


static void ReadData3D(const string& fileName, Data_3D& data);
static void WriteData3D(const string& fileName, const Data_3D& data);

class CRotorVisualPreProcess {
public:
	CRotorVisualPreProcess(const std::string& dictName, const std::string& root, const std::string& block) : 
		dict(dictName), dirRoot(root), dirBlock(block),
		xx0(dict.nxm, dict.nym, dict.nzm),
		yy0(dict.nxm, dict.nym, dict.nzm),
		zz0(dict.nxm, dict.nym, dict.nzm),
		x(dict.nxm, dict.nym, dict.nzm),
		y(dict.nxm, dict.nym, dict.nzm),
		z(dict.nxm, dict.nym, dict.nzm),
		q11(dict.nxm, dict.nym, dict.nzm),
		q12(dict.nxm, dict.nym, dict.nzm),
		q13(dict.nxm, dict.nym, dict.nzm),
		q14(dict.nxm, dict.nym, dict.nzm),
		q15(dict.nxm, dict.nym, dict.nzm),
		q16(dict.nxm, dict.nym, dict.nzm)	
	{
		Init();
	
	}
	
	~CRotorVisualPreProcess() {
		Destroy();
	}
	
	void ReadField() {

		std::string dirField = DirToFileName(dirRoot, "field");

		ReadData3D(DirToFileName(dirField, "x"), x);
		ReadData3D(DirToFileName(dirField, "y"), y);
		ReadData3D(DirToFileName(dirField, "z"), z);
		ReadData3D(DirToFileName(dirField, "xx0"), xx0);
		ReadData3D(DirToFileName(dirField, "yy0"), yy0);
		ReadData3D(DirToFileName(dirField, "yy0"), zz0);
		ReadData3D(DirToFileName(dirField, "q11"), q11);
		ReadData3D(DirToFileName(dirField, "q12"), q12);
		ReadData3D(DirToFileName(dirField, "q13"), q13);
		ReadData3D(DirToFileName(dirField, "q14"), q14);
		ReadData3D(DirToFileName(dirField, "q15"), q15);
		ReadData3D(DirToFileName(dirField, "q16"), q16);

	}
	
	double CubicInterpolation(double* medx, double*medr, double spax, int n1, int n2,int n3);
	void Span(int spa);
	void Lbout();
	
	void Init() {
		dict.readConfig();
	}
	
	void Destroy() {
	}
	
private:
	
	CDictionary dict;
	std::string dirRoot;
	std::string dirBlock;
	Data_3D x;
	Data_3D y;
	Data_3D z;
	Data_3D xx0;
	Data_3D yy0;
	Data_3D zz0;
	Data_3D q11;
	Data_3D q12;
	Data_3D q13;
	Data_3D q14;
	Data_3D q15;
	Data_3D q16;
};

void ReadData3D(const string& fileName, Data_3D& data) {
	fstream f;

	f.open(fileName.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);

	int count = 0;
	for(int i = 1; i < data.nx_; i++) {
		for(int j = 1; j < data.ny_; j++) {
			for(int k = 1; k < data.nz_; k++) {
				f >> data.pData_[i][j][k];
			}
		}
	}
	f.close();	
}

void WriteData3D(const string& fileName, const Data_3D& data) {

	fstream f;

	f.open(fileName.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);

	int count = 0;
	for(int i = 1; i < data.nx_; i++) {
		for(int j = 1; j < data.ny_; j++) {
			for(int k = 1; k < data.nz_; k++) {
				f << data.pData_[i][j][k] << ' ';
				if (++count == 3) {
					count = 0;
					f << endl;
				}

			}

		}
	}

	f.close();
}

void WriteData2D(const string& fileName, const Data_2D& data) {

	fstream f;
	f.open(fileName.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);

	int count = 0;
	for(int i = 1; i < data.nx_; i++) {
		for(int j = 1; j < data.ny_; j++) {
			f << data.pData_[i][j] << ' ';
			if (++count == 3) {
				count = 0;
				f << endl;
			}

		}
	}

	f.close();
}

double CRotorVisualPreProcess::CubicInterpolation(double* medx, double*medr, double spax, int n1, int n2,int n3)
{
	double* h, *g, *lamt, *miu, *a, *m;

	Malloc(h,n1);
	Malloc(g,n1);
	Malloc(lamt,n1);
	Malloc(miu,n1);
	Malloc(a,n1);
	Malloc(m,n1+1);

	double h1, spar;

	spar = 0;
	m[1] = (medr[2] - medr[1]) / (medx[2] - medx[1]);
	m[n1] = (medr[n1] - medr[n1 - 1]) / (medx[n1] - medx[n1 - 1]);

	for (int i = 1; i<n1; i++)
	{
		h[i] = medx[i + 1] - medx[i];
	}

	for (int i = 2; i<n1; i++)
	{
		a[i] = 2.0;
		lamt[i] = h[i] / (h[i] + h[i - 1]);
		miu[i] = 1 - lamt[i];
		g[i] = 3.0*(lamt[i] * (medr[i] - medr[i - 1]) / h[i - 1] + miu[i] * (medr[i + 1] - medr[i]) / h[i]);
	}

	g[2] = g[2] - lamt[2] * m[1];
	g[n1 - 1] = g[n1 - 1] - miu[n1 - 1] * m[n1];

	double* u,*l,*y;
	Malloc(u,n1 - 1);
	Malloc(l,n1 - 1);
	Malloc(y,n1 - 1);
	u[2] = a[2];
	y[2] = g[2];

	for (int i = 3; i<n1 - 1; i++)
	{
		l[i] = lamt[i] / u[i - 1];
		u[i] = a[i] - l[i] * miu[i - 1];
		y[i] = g[i] - l[i] * y[i - 1];
	}
	m[n1 - 2] = y[n3] / u[n3];
	for (int i = n1 - 2; i>0; i--)
	{
		m[i] = (y[i] - miu[i - 1] * m[i + 1]) / u[i];
	}

	for (int j = 1; j<n2 + 1; j++)
	{
		int i;
		if (spax<medx[n1])
		{
			i = 1;
			while (spax>medx[i + 1])
			{
				i++;
			}
		}
		else
			i = n1 - 1;

		h1 = (spax - medx[i]) / h[i];
		spar = spar + medr[i] * (2.0*h1 + 1.0)*(h1 - 1.0)*(h1 - 1.0) + m[i] * h[i] * h1*(h1 - 1.0)*(h1 - 1.0);
		h1 = (medx[i + 1] - spax) / h[i];
		spar = spar + medr[i + 1] * (2.0*h1 + 1.0)*(h1 - 1.0)*(h1 - 1.0) - m[i + 1] * h[i] * h1*(h1 - 1.0)*(h1 - 1.0);
	}

	Free(u);
	Free(l);
	Free(y);
	Free(h);
	Free(g);
	Free(lamt);
	Free(miu);
	Free(a);
	Free(m);
	return spar;
}

void CRotorVisualPreProcess::Span(int spa)
{
	string nnspan;
	
	double  temp, tem, cvl, qq2, qqw, qqr, a, t1;
	double y1, z1, wx, wy, wz;
	double*  hx, *hy, *hz,* hr, *s,* hu, *hv,* hw,* hp, *hpt,* hh, *hmiu, *hwma;
	double**  hxx, **hyy, **hzz;
	double** huu, **hvv, **hww, **hpp, **hppt, **hht, **hwwma, **hmmiu;
	
	int nx = dict.nxm;
	int ny = dict.nym;
	int nz = dict.nzm;
	
	Malloc(hx,ny+2);
	Malloc(hy,ny+2);
	Malloc(hz,ny+2);
	Malloc(hr,ny+2);
	Malloc(s,ny+2);
	Malloc(hv,ny+2);
	Malloc(hu,ny+1);
	Malloc(hw,ny+1);
	Malloc(hp,ny+1);
	Malloc(hpt,ny+1);
	Malloc(hh,ny+1);
	Malloc(hmiu,ny+1);
	Malloc(hwma,ny+1);

	Malloc(hxx,nz+2,nx+2);
	Malloc(hyy,nz+2,nx+2);
	Malloc(hzz,nz+2,nx+2);

	Malloc(huu,nz + 1,nx + 1);
	Malloc(hvv,nz + 1,nx + 1);
	Malloc(hww,nz + 1,nx + 1);
	Malloc(hpp,nz + 1,nx + 1);
	Malloc(hppt,nz + 1,nx + 1);
	Malloc(hht,nz + 1,nx + 1);
	Malloc(hwwma,nz + 1,nx + 1);
	Malloc(hmmiu,nz + 1,nx + 1);

	temp = double(spa) / 100.0;

	for (int i = 1; i<nz + 2; i++)
		for (int k = 1; k<nx + 2; k++)
		{
			s[1] = 0;
			for (int j = 1; j<ny + 2; j++)
			{
				hx[j] = x[i][j][k];
				hy[j] = y[i][j][k];
				hz[j] = z[i][j][k];
				hr[j] = sqrt(hy[j] * hy[j] + hz[j] * hz[j]);
				if (j>1)
				{
					s[j] = s[j - 1] + sqrt(pow(hx[j] - hx[j - 1], 2) + pow(hr[j] - hr[j - 1], 2));
				}
			}

			t1 = s[ny + 1] * temp;
			
			hxx[i][k] = CubicInterpolation(s, hx, t1, ny + 1, 1, dict.nt + 1);
			hyy[i][k] = CubicInterpolation(s, hy, t1, ny + 1, 1, dict.nt + 1);
			hzz[i][k] = CubicInterpolation(s, hz, t1, ny + 1, 1, dict.nt + 1);
		}

	for (int i = 1; i<nz + 1; i++)
		for (int k = 1; k<nx + 1; k++)
		{
			s[1] = 0;
			for (int j = 1; j<ny + 1; j++)
			{
				hx[j] = xx0[i][j][k];
				hy[j] = yy0[i][j][k];
				hz[j] = zz0[i][j][k];
				
				hr[j] = sqrt(hy[j] * hy[j] + hz[j] * hz[j]);
				
				if (j>1)
				{
					s[j] = s[j - 1] + sqrt(pow(hx[j] - hx[j - 1], 2) + pow(hr[j] - hr[j - 1], 2));
				}
			}
			t1 = s[ny] * temp;

			for (int j = 1; j<ny + 1; j++)
			{
				y1 = yy0[i][j][k];
				z1 = zz0[i][j][k];
				hu[j] = q12[i][j][k] / q11[i][j][k];
				hv[j] = q13[i][j][k] / q11[i][j][k];
				hw[j] = q14[i][j][k] / q11[i][j][k];
				wx = hu[j];
				wy = hv[j] + dict.rpm*z1;
				wz = hw[j] - dict.rpm*y1;
				qqw = wx*wx + wy*wy + wz*wz;
				qqr = dict.rpm*dict.rpm*(z1*z1 + y1*y1);
				qq2 = hu[j] * hu[j] + hv[j] * hv[j] + hw[j] * hw[j];
				hp[j] = 0.40*(q15[i][j][k] - 0.50*q11[i][j][k] * qq2);
				tem = hp[j] / (q11[i][j][k] *dict.rg);		
				cvl = dict.cvl0*pow((tem / dict.t0), 1.5)*(dict.t0 + dict.ts) / (tem + dict.ts);
				a = 1.40*hp[j] / q11[i][j][k];
				hwma[j] = sqrt(qqw / a);
				hh[j] = (q15[i][j][k] + hp[j]) / q11[i][j][k];
				hpt[j] = hp[j] * pow(1.0 - 0.50*qq2 / hh[j], -3.5);
				hmiu[j] = q16[i][j][k] / cvl;
			}

			//nn = 1 , nn < nt + 1, ++nn  (s, hu, t1, ny, 1, nn)
			huu[i][k] = CubicInterpolation(s, hu, t1, ny, 1, 1);
			hvv[i][k] = CubicInterpolation(s, hv, t1, ny, 1, 1);
			hww[i][k] = CubicInterpolation(s, hw, t1, ny, 1, 1);
			hpp[i][k] = CubicInterpolation(s, hp, t1, ny, 1, 1);
			hppt[i][k] = CubicInterpolation(s, hpt, t1, ny, 1, 1);
			hht[i][k] = CubicInterpolation(s, hh, t1, ny, 1, 1);
			hwwma[i][k] = CubicInterpolation(s, hwma, t1, ny, 1,1);
			hmmiu[i][k] = CubicInterpolation(s, hmiu, t1, ny, 1,1);
		
		}

		nnspan = Int_to_string(spa);

        Data_2D hxx_2D(hxx, nz + 2, nz + 2, "x");
        Data_2D hyy_2D(hyy, nz + 2, nx + 2, "y");
        Data_2D hzz_2D(hzz, nz + 2, nz + 2, "z");
        Data_2D huu_2D(huu, nz + 1, nx + 1, "u");
        Data_2D hvv_2D(hvv, nz + 1, nx + 1, "v");
        Data_2D hww_2D(hww, nz + 1, nx + 1, "w");
        Data_2D hwwma_2D(hwwma, nz + 1, nx + 1, "wma");
        Data_2D hpp_2D(hpp, nz + 1, nx + 1, "pressure");
        Data_2D hppt_2D(hppt, nz + 1, nx + 1, "pt");
        Data_2D hht_2D(hht, nz + 1, nx + 1, "zht");
        Data_2D hmmiu_2D(hmmiu, nz + 1, nx + 1, "ut");

		std::string dirVisual = dirRoot;
		std::string dirSpan = CreateDir(dirVisual, "span");
		std::string dirSpanBlock = CreateDir(dirSpan, dirBlock);

        WriteData2D(dirSpanBlock + hxx_2D.name_, hxx_2D);
        WriteData2D(dirSpanBlock + hyy_2D.name_, hyy_2D);
        WriteData2D(dirSpanBlock + hzz_2D.name_, hzz_2D);
        WriteData2D(dirSpanBlock + huu_2D.name_, huu_2D);
        WriteData2D(dirSpanBlock + hvv_2D.name_, hvv_2D);
        WriteData2D(dirSpanBlock + hww_2D.name_, hww_2D);
        WriteData2D(dirSpanBlock + hwwma_2D.name_, hwwma_2D);
        WriteData2D(dirSpanBlock + hpp_2D.name_, hpp_2D);
        WriteData2D(dirSpanBlock + hppt_2D.name_, hppt_2D);
        WriteData2D(dirSpanBlock + hht_2D.name_, hht_2D);
        WriteData2D(dirSpanBlock + hmmiu_2D.name_, hmmiu_2D);


	Free(hx);
	Free(hy);
	Free(hz);
	Free(hr);
	Free(s);
	Free(hu);
	Free(hv);
	Free(hw);
	Free(hp);
	Free(hpt);
	Free(hh);
	Free(hmiu);
	Free(hwma);
	Free(hxx);
	Free(hyy);
	Free(hzz);
	Free(huu);
	Free(hvv);
	Free(hww);
	Free(hpp);
	Free(hppt);
	Free(hht);
	Free(hwwma);
	Free(hmmiu);
}

void CRotorVisualPreProcess::Lbout() {
	double qq2, cvl, rr;
	double vx, vy, vz, y1, z1, dim;
	double qqw, a;
	double wx, wy, wz;

	int nx = dict.nxm;
	int ny = dict.nym;
	int nz = dict.nzm;


	double ***t;
	double ***pvx, ***pvy, ***pvz, ***wma, ***p;
	
	for (int k = 1; k<nz + 1; k++)
		for (int j = 1; j<ny + 1; j++)
			for (int i = 1; i<nx + 1; i++)
			{
				dim = q11[k][j][i];
				pvx[k][j][i] = q12[k][j][i] / dim;
				pvy[k][j][i] = q13[k][j][i] / dim;
				pvz[k][j][i] = q14[k][j][i] / dim;
				vx = pvx[k][j][i];
				vy = pvy[k][j][i];
				vz = pvz[k][j][i];
				y1 = yy0[k][j][i];
				z1 = zz0[k][j][i];
				rr = sqrt(y1*y1 + z1*z1);

				qq2 = vx*vx + vy*vy + vz*vz;
				p[k][j][i] = 0.4*(q15[k][j][i] - 0.5*dim*qq2);
				t[k][j][i] = p[k][j][i] / (dim*dict.rg);
				cvl = dict.cvl0*pow(t[k][j][i] / dict.t0, 1.5)*(dict.t0 + dict.ts) / (t[k][j][i] + dict.ts);
				q16[k][j][i] = max(q16[k][j][i], pow(10.0, -4)*cvl);
				
				y1 = yy0[i][j][k];
				z1 = zz0[i][j][k];
				wx = vx;
				wy = vy + dict.rpm*z1;
				wz = vz - dict.rpm*y1;
				qqw = wx*wx + wy*wy + wz*wz;
				a = 1.40*p[i][j][k] / q11[i][j][k];
				wma[i][j][k] = sqrt(qqw / a);
			}

	Data_3D pvx_3D(pvx, nz + 1, ny + 1, nx + 1, "pvx");
	Data_3D pvy_3D(pvy, nz + 1, ny + 1, nx + 1, "pvy");
	Data_3D pvz_3D(pvz, nz + 1, ny + 1, nx + 1, "pvz");
	Data_3D wma_3D(wma, nz + 1, ny + 1, nx + 1, "wma");
	Data_3D p_3D(p, nz + 1, ny + 1, nx + 1, "p");
	Data_3D q11_3D(q11.pData_, nz + 1, ny + 1, nx + 1, "q11");
	Data_3D q16_3D(q16.pData_, nz + 1, ny + 1, nx + 1, "q16");

	std::string dirVisual = dirRoot;
	std::string dirLbout = CreateDir(dirRoot, "lbout");
	std::string dirLboutBlock = CreateDir(dirLbout, dirBlock);

	WriteData3D(dirLboutBlock + pvx_3D.name_, pvx_3D);
	WriteData3D(dirLboutBlock + pvy_3D.name_, pvy_3D);
	WriteData3D(dirLboutBlock + pvz_3D.name_, pvz_3D);
	WriteData3D(dirLboutBlock + wma_3D.name_, wma_3D);
	WriteData3D(dirLboutBlock + p_3D.name_, p_3D);
	WriteData3D(dirLboutBlock + q11_3D.name_, q11_3D);
	WriteData3D(dirLboutBlock + q16_3D.name_, q16_3D);
	
	Free(pvx);
	Free(pvy);	
	Free(pvz);	
	Free(wma);
	Free(p);	
}
