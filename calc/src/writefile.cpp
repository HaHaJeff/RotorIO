/*---------------------------------------------------------------------------*\


  \*---------------------------------------------------------------------------*/


#include "writefile.h"

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <cassert>
#include <vector>
#include <omp.h>
#include <string.h>
#include <errno.h>

#include "IOStrategy/MPIIOStrategy.h"
#include "IOStrategy/MPIIO.h"

//***************************************************************************


void CWriteFile::ReadData3D(const string& fileName, Data_3D& data, const MPI_Comm& comm) {

	
	int fd = open(fileName.c_str(), O_RDONLY);
	assert(fd != -1);

	double *ptr = data.pData_;
	int left = data.nx_ * data.ny_ * data.nz_ * sizeof(double);
	int done = 0;
	int offset = 0;

	while ( left > 0 && (done = read(fd, ptr + offset, left - done)) != -1 )  {
		offset += done/sizeof(double);
		left -= done;
	}

	if (left > 0) {
		fprintf(stderr, "read '%s', data '%s', bytes '%d' error '%s'\n", fileName.c_str(), data.name_.c_str(), left, strerror(errno));
		return;
	}

	close(fd);
}

void CWriteFile::WriteData3D(const string& fileName, const Data_3D& data, const MPI_Comm& com) {

	int fd = open(fileName.c_str(), O_CREAT | O_WRONLY, 0666);
	assert(fd != -1);

	//FIXME: data[1][1][1] and next need to write
	//double *ptr = &data.pData_[0][0][0];
	//int left = data.nx_ * data.ny_ * data.nz_;
	//double *ptr = &data.pData_[1][1][1];
	double *ptr = data.pData_;
	int left = (data.nx_ - 1) * (data.ny_ - 1) * (data.nz_ - 1)* sizeof(double);
	int done = 0;
	int offset = 0;

	while ( left > 0 && (done = write(fd, ptr + offset, left - done)) != -1 )  {
		offset += done/sizeof(double);
		left -= done;
	}

	if (left > 0) {
		fprintf(stderr, "write '%s', data '%s', bytes '%d' error '%s'\n", fileName.c_str(), data.name_.c_str(), left, strerror(errno));
		return;
	}

	close(fd);
}

/*
void CWriteFile::WriteData2D(const string& fileName, const Data_2D& data) {

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
*/

void CWriteFile::WriteHeader(const string& fileName, const string& zoneName, bool isCell) {

	fstream f;
	f.open(fileName.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);

	f << "VARIABLES= \"x\", \"y\", \"z\", \"u\", \"v\", \"w\", \"wma\", \"pressure\", \"density\", \"ut\" " << endl;
	f << "ZONE T=" << "\"" <<zoneName << "\""  << endl;

	if(isCell) {
		f << "I = " << nx << ",J = " << ny << ",K = " << nz << endl;
		f << "DATAPACKING=BLOCK" << endl;
		//f << "DATAPACKING=BLOCK, VARLOCATION=([4-10]=CELLCENTERED),SOLUTIONTIME=" << setw(8) << setprecision(4) << endl;
	} else {
		f << "I = " << nx + 1 << ",J = " << ny + 1 << ",K = " << nz + 1 << endl;
		f << "DATAPACKING=BLOCK" << endl;
	}

	f.close();
}


void CWriteFile::InitPressure(int n1)
{
	double qq2, cvl, sir, cor;
	double vx, vy, vz, y1, z1, rr, dim;
	for (int k = 1; k<nz + 1; k++)
		for (int j = 1; j<ny + 1; j++)
			for (int i = 1; i<nx + 1; i++)
			{
				dim = pfield->q11[n1][k][j][i];
				pfield->pvx[k][j][i] = pfield->q12[n1][k][j][i] / dim;
				pfield->pvy[k][j][i] = pfield->q13[n1][k][j][i] / dim;
				pfield->pvz[k][j][i] = pfield->q14[n1][k][j][i] / dim;
				vx = pfield->pvx[k][j][i];
				vy = pfield->pvy[k][j][i];
				vz = pfield->pvz[k][j][i];
				y1 = pgrid->yy0[k][j][i];
				z1 = pgrid->zz0[k][j][i];
				rr = sqrt(y1*y1 + z1*z1);
				sir = z1 / rr;
				cor = y1 / rr;
				pfield->vth[k][j][i] = vz*cor - vy*sir;    //\D6\DC\CF\F2\CBٶ\C8
				pfield->vre[k][j][i] = vz*sir + vy*cor;    //\BE\B6\CF\F2\CBٶ\C8
				qq2 = vx*vx + vy*vy + vz*vz;
				pfield->p[k][j][i] = 0.4*(pfield->q15[n1][k][j][i] - 0.5*dim*qq2);
				pfield->t[k][j][i] = pfield->p[k][j][i] / (dim*pdict->rg);
				cvl = pdict->cvl0*pow(pfield->t[k][j][i] / pdict->t0, 1.5)*(pdict->t0 + pdict->ts) / (pfield->t[k][j][i] + pdict->ts);
				pfield->q16[n1][k][j][i] = max(pfield->q16[n1][k][j][i], pow(10.0, -4)*cvl);
			}
}

void CWriteFile::inlout()
{
	string fInlout = fOut + "inlout/";
	int mml=1;
	double  vf1, vf2, qin, qout1, qq2, h0, t2, p2, pout, tout, d1, d2, eff, temp, ftt, fpp;
	double* rr0, *qout, *fp, *ft, *qin0, *qout0, *d10, *d20, *eff0;
	double qinn, qoutt, d11, d22, efff;
	double pp, vx, vy, vz, y1, z1;
	fstream f4, f5, f8, f9, f11;
	string zonename;
	if (myid == 0)
	{
		string fF4 = fInlout + "inflow.dat";
		string fF5 = fInlout + "outflow.dat";
		string fF8 = fInlout + "pbi.dat";
		string fF9 = fInlout + "tbi.dat";
		string fF11 = fInlout + "eff.dat";
		LogInfo::Log("Output flow files");
		f4.open(fF4.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
		f5.open(fF5.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
		f8.open(fF8.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
		f9.open(fF9.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
		f11.open(fF11.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
	}

	for (int nn = 1; nn<pdict->nt + 1; nn++)
	{
		qin = 0;
		for (int j = 1; j<nz + 1; j++)
			for (int k = 1; k<ny + 1; k++)
			{
				vf1 = -(pgrid->s2x[j][k][1] * pfield->q12[nn][j][k][1] + pgrid->s2y[j][k][1] * pfield->q13[nn][j][k][1] + pgrid->s2z[j][k][1] * pfield->q14[nn][j][k][1]);
				qin = qin + vf1;
			}

		if (myid == 0)
		{
			Malloc(qin0,pdict->lbb);
			qin0[0] = qin;
			qinn = 0;
			for (int j = 1; j<pdict->lbb; j++)
				MPI_Recv(&qin0[j], 1, MPI_DOUBLE, j, 3, MPI_COMM_WORLD, &status);

			for (int j = 0; j<pdict->lbb; j++)
				qinn = qinn + qin0[j];

			f4 << setw(6) << setprecision(17) << double(nn - 1) / double(pdict->nt);
			f4 << setw(15) << setprecision(17) << qinn;
			Free(qin0);
		}
		else
			MPI_Send(&qin, 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);

	}

	Malloc(rr0,ny+1);
	Malloc(qout,ny+1);
	Malloc(fp,ny+1);
	Malloc(ft,ny+1);

	for (int j = 1; j<ny + 1; j++)
	{
		y1 = pgrid->yy0[1][j][nx];
		z1 = pgrid->zz0[1][j][nx];
		rr0[j] = sqrt(y1*y1 + z1*z1);
	}

	for (int nn = 1; nn<pdict->nt + 1; nn++)
	{
		for (int j = 1; j<ny + 1; j++)
		{
			qout[j] = 0.0;
			fp[j] = 0.0;
			ft[j] = 0.0;
		}

		qout1 = 0.0;
		fpp = 0.0;
		ftt = 0.0;

		zonename = Int_to_string(nn);

		fstream f44, f55, f66;
		string file = fInlout + "pbispan-timel" + zonename + ".dat";
		f44.open(file.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
		file = fInlout + "tbispan-timel" + zonename + ".dat";
		f55.open(file.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
		file = fInlout + "effspan-timel" + zonename + ".dat";
		f66.open(file.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);

		for (int j = 1; j<ny + 1; j++)
		{
			temp = (rr0[j] - rr0[1]) / (rr0[ny] - rr0[1])*100.0;
			for (int k = 1; k<nz + 1; k++)
			{
				vf2 = -(pgrid->s2x[k][j][nx + 1] * pfield->q12[nn][k][j][nx] + pgrid->s2y[k][j][nx + 1] * pfield->q13[nn][k][j][nx] + pgrid->s2z[k][j][nx + 1] * pfield->q14[nn][k][j][nx]);
				vx = pfield->q12[nn][k][j][nx] / pfield->q11[nn][k][j][nx];
				vy = pfield->q13[nn][k][j][nx] / pfield->q11[nn][k][j][nx];
				vz = pfield->q14[nn][k][j][nx] / pfield->q11[nn][k][j][nx];
				qq2 = vx*vx + vy*vy + vz*vz;
				pp = 0.40*(pfield->q15[nn][k][j][nx] - 0.50*pfield->q11[nn][k][j][nx] * qq2);
				h0 = (pfield->q15[nn][k][j][nx] + pp) / pfield->q11[nn][k][j][nx];
				t2 = h0 / pdict->cp;
				p2 = pp*pow(1.0 - 0.50*qq2 / h0, -3.5);
				ft[j] = ft[j] + t2*vf2;
				fp[j] = fp[j] + p2*vf2;
				qout[j] = qout[j] + vf2;
			}

			pout = fp[j] / qout[j];
			tout = ft[j] / qout[j];
			d1 = pout / pdict->pt;
			d2 = tout / (pdict->ht / pdict->cp);
			eff = (pow(d1, 2.0 / 7.0) - 1.0) / (d2 - 1.0)*100.0;

			if (myid == 0)
			{

				Malloc(d10,pdict->lbb);
				Malloc(d20,pdict->lbb);
				Malloc(eff0,pdict->lbb);

				d10[0] = d1;
				d20[0] = d2;
				eff0[0] = eff;

				d11 = 0.0;
				d22 = 0.0;
				efff = 0.0;

				for (int i = 1; i<pdict->lbb; i++)
				{
					MPI_Recv(&d10[i], 1, MPI_DOUBLE, i, 3, MPI_COMM_WORLD, &status);
					MPI_Recv(&d20[i], 1, MPI_DOUBLE, i, 4, MPI_COMM_WORLD, &status);
					MPI_Recv(&eff0[i], 1, MPI_DOUBLE, i, 5, MPI_COMM_WORLD, &status);
				}

				for (int i = 0; i<pdict->lbb; i++)
				{
					d11 = d11 + d10[i];
					d22 = d22 + d20[i];
					efff = efff + eff0[i];
				}
				f44 << setw(15) << setprecision(6) << d11 / (mml*pdict->lbb) << "         ";
				f44 << setw(15) << setprecision(6) << temp << endl;
				f55 << setw(15) << setprecision(6) << d11 / (mml*pdict->lbb) << "        ";
				f55 << setw(15) << setprecision(6) << temp << endl;
				f66 << setw(15) << setprecision(6) << d11 / (mml*pdict->lbb) << "        ";
				f66 << setw(15) << setprecision(6) << temp << endl;
				Free(d10);
				Free(d20);
				Free(eff0);
			}
			else
			{
				MPI_Send(&d1, 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
				MPI_Send(&d2, 1, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);
				MPI_Send(&eff, 1, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);
			}

			ftt = ftt + ft[j];
			fpp = fpp + fp[j];
			qout1 = qout1 + qout[j];

		}
		f44.close();
		f55.close();
		f66.close();

		pout = fpp / qout1;
		tout = ftt / qout1;
		d1 = pout / pdict->pt;
		d2 = tout / (pdict->ht / pdict->cp);
		eff = (pow(d1, 2.0 / 7.0) - 1.0) / (d2 - 1.0)*100.0;

		if (myid == 0)
		{

			Malloc(qout0,pdict->lbb);
			Malloc(d10,pdict->lbb);
			Malloc(d20,pdict->lbb);
			Malloc(eff0,pdict->lbb);


			qout0[0] = qout1;
			d10[0] = d1;
			d20[0] = d2;
			eff0[0] = eff;

			qoutt = 0.0;
			d11 = 0.0;
			d22 = 0.0;
			efff = 0.0;

			for (int j = 1; j<pdict->lbb; j++)
			{
				MPI_Recv(&qout0[j], 1, MPI_DOUBLE, j, 2, MPI_COMM_WORLD, &status);
				MPI_Recv(&d10[j], 1, MPI_DOUBLE, j, 3, MPI_COMM_WORLD, &status);
				MPI_Recv(&d20[j], 1, MPI_DOUBLE, j, 4, MPI_COMM_WORLD, &status);
				MPI_Recv(&eff0[j], 1, MPI_DOUBLE, j, 5, MPI_COMM_WORLD, &status);
			}

			for (int j = 0; j<pdict->lbb; j++)
			{
				qoutt = qoutt + qout0[j];
				d11 = d11 + d10[j];
				d22 = d22 + d20[j];
				efff = efff + eff0[j];
			}

			f5 << setw(6) << setprecision(3) << double(nn - 1) / double(pdict->nt);
			f5 << setw(15) << setprecision(6) << qoutt;
			f8 << setw(6) << setprecision(3) << double(nn - 1) / double(pdict->nt);
			f8 << setw(15) << setprecision(6) << d11 / pdict->lbb;
			f9 << setw(6) << setprecision(3) << double(nn - 1) / double(pdict->nt);
			f9 << setw(15) << setprecision(6) << d22 / pdict->lbb;
			f11 << setw(6) << setprecision(3) << double(nn - 1) / double(pdict->nt);
			f11 << setw(15) << setprecision(6) << efff / pdict->lbb;
			Free(d10);
			Free(d20);
			Free(eff0);
			Free(qout0);
		}
		else
		{
			MPI_Send(&qout1, 1, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD);
			MPI_Send(&d1, 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
			MPI_Send(&d2, 1, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);
			MPI_Send(&eff, 1, MPI_DOUBLE, 0, 5, MPI_COMM_WORLD);
		}

	}
	Free(rr0);
	Free(qout);
	Free(fp);
	Free(ft);


}
void CWriteFile::RelativeMa(int n1)
{
	double qqw, a;
	double y1, z1, vx, vy, vz, wx, wy, wz;
	for (int i = 1; i<nz + 1; i++)
		for (int j = 1; j<ny + 1; j++)
			for (int k = 1; k<nx + 1; k++)
			{
				y1 = pgrid->yy0[i][j][k];
				z1 = pgrid->zz0[i][j][k];
				vx = pfield->pvx[i][j][k];
				vy = pfield->pvy[i][j][k];
				vz = pfield->pvz[i][j][k];
				wx = vx;
				wy = vy + pdict->rpm*z1;
				wz = vz - pdict->rpm*y1;
				qqw = wx*wx + wy*wy + wz*wz;
				a = 1.40*pfield->p[i][j][k] / pfield->q11[n1][i][j][k];
				wma[i][j][k] = sqrt(qqw / a);
			}
}

void CWriteFile::lbout()
{
	double ***pvxn, ***pvyn, ***pvzn, ***wman, ***pn, ***q11n, ***q16n;

	MPI_Comm comm = MPI_COMM_WORLD;
	
//	MPI_Comm_split(MPI_COMM_WORLD, myid / 1, myid % 1, &comm);

	string zonename;

	if(myid==0)LogInfo::Log("Output out-* files");

	double tStart, tEnd;
	tStart = MPI_Wtime();
	
	if (myid == 0) {
		std::cout << "start" << std::endl;
	}

	for (int nn = 1; nn<pdict->nt + 1; nn++)
	{
		InitPressure(nn);
		RelativeMa(nn);

		{
		//	Time *pTime = NULL;
		//	if (myid == 0)
		//		pTime = new Time();
			string fLbout = fOut + "lbout-cell/";
			string fTimeStepDir = fLbout + Int_to_string(nn) + "/"; 
			string fBlockDir = CreateDir("block-" + Int_to_string(myid), fTimeStepDir);

			if(myid == 0) {
				cout << fBlockDir << endl;
			}

			zonename = Int_to_string(nn);
			zonename = "out_" + zonename;

			// format is nz ny nx
			Data_3D x_3D(&(pfield->mesh->x[1][1][1]), nz + 2, ny + 2, nx + 2, "x");
			Data_3D y_3D(&(pfield->mesh->y[1][1][1]), nz + 2, ny + 2, nx + 2, "y");
			Data_3D z_3D(&(pfield->mesh->z[1][1][1]), nz + 2, ny + 2, nx + 2, "z");
			Data_3D pvx_3D(&(pfield->pvx[1][1][1]), nz + 1, ny + 1, nx + 1, "pvx");
			Data_3D pvy_3D(&(pfield->pvy[1][1][1]), nz + 1, ny + 1, nx + 1, "pvy");
			Data_3D pvz_3D(&(pfield->pvz[1][1][1]), nz + 1, ny + 1, nx + 1, "pvz");
			Data_3D wma_3D(&wma[1][1][1], nz + 1, ny + 1, nx + 1, "wma");
			Data_3D p_3D(&(pfield->p[1][1][1]), nz + 1, ny + 1, nx + 1, "p");
			Data_3D q11_3D((&(pfield->q11)[nn][1][1][1]), nz + 1, ny + 1, nx + 1, "q11");
			Data_3D q16_3D((&(pfield->q16)[nn][1][1][1]), nz + 1, ny + 1, nx + 1, "q16");

			std::vector<Data_3D> vec;
			vec.push_back(x_3D);
			vec.push_back(y_3D);
			vec.push_back(z_3D);
			vec.push_back(pvx_3D);
			vec.push_back(pvy_3D);
			vec.push_back(pvz_3D);
			vec.push_back(wma_3D);
			vec.push_back(p_3D);
			vec.push_back(q11_3D);
			vec.push_back(q16_3D);

			for (int i = 0; i < vec.size(); i++) {
				WriteData3D(fBlockDir + vec[i].name_, vec[i], comm);
			}

/*
			MPIIO mpiio(comm);
			MPIIOStrategy* pStrategy = mpiio.GetIOStrategy(static_cast<TYPE>(0));
//#pragma omp parallel for
			for (int i = 0; i < vec.size(); i++) {
//				WriteData3D(fBlockDir + vec[i].name_, vec[i]);
				pStrategy->Open(fBlockDir + vec[i].name_);
				pStrategy->Write(vec[i]);
				pStrategy->Close();
			}

			delete(pStrategy);
			pStrategy = NULL;
*/

	//		if (myid == 0)
	//			delete(pTime);

		}
	}

	tEnd = MPI_Wtime();

	if (myid == 0) 
		std::cout << "Lbout time: " << tEnd - tStart << endl;

}

double CWriteFile::CubicInterpolation(double* medx, double*medr, double spax, int n1, int n2,int n3)
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
void CWriteFile::span(int spa)
{
	string fSpan = fOut + "span/";
	if(myid==0)LogInfo::Log("Output "+Int_to_string(spa)+"%span-* files");
	string nnspan, zonename;
	double  temp, tem, cvl, qq2, qqw, qqr, a, t1;
	double y1, z1, wx, wy, wz;
	double*  hx, *hy, *hz,* hr, *s,* hu, *hv,* hw,* hp, *hpt,* hh, *hmiu, *hwma;
	double**  hxx, **hyy, **hzz;
	double*** huu, ***hvv, ***hww, ***hpp, ***hppt, ***hht, ***hwwma, ***hmmiu;
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

	Malloc(huu,pdict->nt + 1,nz + 1,nx + 1);
	Malloc(hvv,pdict->nt + 1,nz + 1,nx + 1);
	Malloc(hww,pdict->nt + 1,nz + 1,nx + 1);
	Malloc(hpp,pdict->nt + 1,nz + 1,nx + 1);
	Malloc(hppt,pdict->nt + 1,nz + 1,nx + 1);
	Malloc(hht,pdict->nt + 1,nz + 1,nx + 1);
	Malloc(hwwma,pdict->nt + 1,nz + 1,nx + 1);
	Malloc(hmmiu,pdict->nt + 1,nz + 1,nx + 1);


	temp = double(spa) / 100.0;

	for (int i = 1; i<nz + 2; i++)
		for (int k = 1; k<nx + 2; k++)
		{
			s[1] = 0;
			for (int j = 1; j<ny + 2; j++)
			{
				hx[j] = (pfield->mesh->x)[i][j][k];
				hy[j] = (pfield->mesh->y)[i][j][k];
				hz[j] = (pfield->mesh->z)[i][j][k];
				hr[j] = sqrt(hy[j] * hy[j] + hz[j] * hz[j]);
				if (j>1)
				{
					s[j] = s[j - 1] + sqrt(pow(hx[j] - hx[j - 1], 2) + pow(hr[j] - hr[j - 1], 2));
				}
			}

			t1 = s[ny + 1] * temp;
			hxx[i][k] = CubicInterpolation(s, hx, t1, ny + 1, 1,pdict->nt + 1);
			hyy[i][k] = CubicInterpolation(s, hy, t1, ny + 1, 1,pdict->nt + 1);
			hzz[i][k] = CubicInterpolation(s, hz, t1, ny + 1, 1,pdict->nt + 1);
		}

	for (int i = 1; i<nz + 1; i++)
		for (int k = 1; k<nx + 1; k++)
		{
			s[1] = 0;
			for (int j = 1; j<ny + 1; j++)
			{
				hx[j] = pgrid->xx0[i][j][k];
				hy[j] = pgrid->yy0[i][j][k];
				hz[j] = pgrid->zz0[i][j][k];
				hr[j] = sqrt(hy[j] * hy[j] + hz[j] * hz[j]);
				if (j>1)
				{
					s[j] = s[j - 1] + sqrt(pow(hx[j] - hx[j - 1], 2) + pow(hr[j] - hr[j - 1], 2));
				}
			}
			t1 = s[ny] * temp;
			for (int nn = 1; nn<pdict->nt + 1; nn++)
			{
				for (int j = 1; j<ny + 1; j++)
				{
					y1 = pgrid->yy0[i][j][k];
					z1 = pgrid->zz0[i][j][k];
					hu[j] = pfield->q12[nn][i][j][k] / pfield->q11[nn][i][j][k];
					hv[j] = pfield->q13[nn][i][j][k] / pfield->q11[nn][i][j][k];
					hw[j] = pfield->q14[nn][i][j][k] / pfield->q11[nn][i][j][k];
					wx = hu[j];
					wy = hv[j] + pdict->rpm*z1;
					wz = hw[j] - pdict->rpm*y1;
					qqw = wx*wx + wy*wy + wz*wz;
					qqr = pdict->rpm*pdict->rpm*(z1*z1 + y1*y1);
					qq2 = hu[j] * hu[j] + hv[j] * hv[j] + hw[j] * hw[j];
					hp[j] = 0.40*(pfield->q15[nn][i][j][k] - 0.50*pfield->q11[nn][i][j][k] * qq2);
					tem = hp[j] / (pfield->q11[nn][i][j][k] * pdict->rg);
					cvl = pdict->cvl0*pow((tem / pdict->t0), 1.5)*(pdict->t0 + pdict->ts) / (tem + pdict->ts);
					a = 1.40*hp[j] / pfield->q11[nn][i][j][k];
					hwma[j] = sqrt(qqw / a);
					hh[j] = (pfield->q15[nn][i][j][k] + hp[j]) / pfield->q11[nn][i][j][k];
					hpt[j] = hp[j] * pow(1.0 - 0.50*qq2 / hh[j], -3.5);
					hmiu[j] = pfield->q16[nn][i][j][k] / cvl;
				}

				huu[nn][i][k] = CubicInterpolation(s, hu, t1, ny, 1,nn);
				hvv[nn][i][k] = CubicInterpolation(s, hv, t1, ny, 1,nn);
				hww[nn][i][k] = CubicInterpolation(s, hw, t1, ny, 1,nn);
				hpp[nn][i][k] = CubicInterpolation(s, hp, t1, ny, 1,nn);
				hppt[nn][i][k] = CubicInterpolation(s, hpt, t1, ny, 1,nn);
				hht[nn][i][k] = CubicInterpolation(s, hh, t1, ny, 1,nn);
				hwwma[nn][i][k] = CubicInterpolation(s, hwma, t1, ny, 1,nn);
				hmmiu[nn][i][k] = CubicInterpolation(s, hmiu, t1, ny, 1,nn);
			}
		}

	//	*****************

	nnspan = Int_to_string(spa);

	fstream f21;
	string file = fSpan + nnspan + "%span-" + id_m + "myid.dat";
	f21.open(file.c_str(), std::fstream::in | std::fstream::out | std::fstream::app);
	f21 << "VARIABLES= \"x\", \"y\", \"z\", \"u\", \"v\", \"w\", \"wma\", \"pressure\", \"density\", \"ut\" " << endl;


	for (int nn = 1; nn < pdict->nt + 1; nn++)
	{

		zonename = Int_to_string(nn);
		zonename = "out_" + zonename;
		int count = 0;
		if (nn == 1)
		{
			f21 << "ZONE T=" << "\"" <<zonename << "\""  << endl;
			f21 << "I = " << nx + 1 << ",J = " << ny + 1 << ",K = " << nz + 1 << endl;
			f21 << "DATAPACKING=BLOCK, VARLOCATION=([4-10]=CELLCENTERED),SOLUTIONTIME=" << setw(8) << setprecision(4) << double(nn - 1) / double(pdict->nt) << endl;

			for (int i = 1; i<nz + 2; i++)
				for (int k = 1; k<nx + 2; k++)
				{
					f21 << setprecision(15) << hxx[i][k] << "      ";
					count++;
					if (count == 3)
					{
						f21 << endl;
						count = 0;
					}
					//							    if(myid==0)cout<<i<<' '<<k<<endl;
				}

			count = 0;
			for (int i = 1; i<nz + 2; i++)
				for (int k = 1; k<nx + 2; k++)
				{
					f21 << setprecision(15) << hyy[i][k] << "      ";
					count++;
					if (count == 3)
					{
						f21 << endl;
						count = 0;
					}
				}
			count = 0;
			for (int i = 1; i<nz + 2; i++)
				for (int k = 1; k<nx + 2; k++)
				{
					f21 << setprecision(15) << hzz[i][k] << "      ";
					count++;
					if (count == 3)
					{
						f21 << endl;
						count = 0;
					}
				}

			count = 0;
			for (int i = 1; i<nz + 1; i++)
				for (int k = 1; k<nx + 1; k++)
				{
					f21 << setprecision(15) << huu[nn][i][k] << "      ";
					count++;
					if (count == 3)
					{
						f21 << endl;
						count = 0;
					}
				}
			count = 0;
			for (int i = 1; i<nz + 1; i++)
				for (int k = 1; k<nx + 1; k++)
				{
					f21 << setprecision(15) << hvv[nn][i][k] << "      ";
					count++;
					if (count == 3)
					{
						f21 << endl;
						count = 0;
					}
				}
			count = 0;
			for (int i = 1; i<nz + 1; i++)
				for (int k = 1; k<nx + 1; k++)
				{
					count++;
					if (count == 3)
					{
						f21 << endl;
						count = 0;
					}
					f21 << setprecision(15) << hww[nn][i][k] << "      ";
					count = 0;
				}

			for (int i = 1; i<nz + 1; i++)
				for (int k = 1; k<nx + 1; k++)
				{
					f21 << setprecision(15) << hwwma[nn][i][k] << "      ";
					count++;
					if (count == 3)
					{
						f21 << endl;
						count = 0;
					}
				}
			count = 0;
			for (int i = 1; i<nz + 1; i++)
				for (int k = 1; k<nx + 1; k++)
				{
					f21 << setprecision(15) << hpp[nn][i][k] << "      ";
					count++;
					if (count == 3)
					{
						f21 << endl;
						count = 0;
					}
				}
			count = 0;
			for (int i = 1; i<nz + 1; i++)
				for (int k = 1; k<nx + 1; k++)
				{
					f21 << setprecision(15) << hppt[nn][i][k] << "      ";
					count++;
					if (count == 3)
					{
						f21 << endl;
						count = 0;
					}
				}
			count = 0;
			for (int i = 1; i<nz + 1; i++)
				for (int k = 1; k<nx + 1; k++)
				{
					f21 << setprecision(15) << hht[nn][i][k] << "      ";
					count++;
					if (count == 3)
					{
						f21 << endl;
						count = 0;
					}
				}
			count = 0;
			for (int i = 1; i<nz + 1; i++)
				for (int k = 1; k<nx + 1; k++)
				{
					f21 << setprecision(15) << hmmiu[nn][i][k] << "      ";
					count++;
					if (count == 3)
					{
						f21 << endl;
						count = 0;
					}
				}
		}
		else
		{

			f21 << "ZONE T=" << "\"" <<zonename << "\""  << endl;
			f21 << "I = " << nx + 1 << ",J = " << ny + 1 << ",K = " << nz + 1 << endl;
			f21 << "DATAPACKING=BLOCK, VARLOCATION=([4-10]=CELLCENTERED),SOLUTIONTIME=" << setw(8) << setprecision(4) << double(nn - 1) / double(pdict->nt) << endl;

			//	f21 << "ZONE T=" << zonename << "I=" << nx + 1 << "J=" << ny + 1 << "K=" << nz + 1 << "DATAPACKING=BLOCK, VARSHARELIST=([1-3]=1),VARLOCATION=([4-10]=CELLCENTERED),SOLUTIONTIME="
			//<< setw(8) << setprecision(4) << double(nn - 1) / double(pdict->nt);

			for (int i = 1; i<nz + 1; i++)
				for (int k = 1; k<nx + 1; k++)
				{
					f21 << setprecision(15) << huu[nn][i][k] << "      ";
					count++;
					if (count == 3)
					{
						f21 << endl;
						count = 0;
					}
				}
			count = 0;
			for (int i = 1; i<nz + 1; i++)
				for (int k = 1; k<nx + 1; k++)
				{
					f21 << setprecision(15) << hvv[nn][i][k] << "      ";
					count++;
					if (count == 3)
					{
						f21 << endl;
						count = 0;
					}

				}
			count = 0;
			for (int i = 1; i<nz + 1; i++)
				for (int k = 1; k<nx + 1; k++)
				{
					f21 << setprecision(15) << hww[nn][i][k] << "      ";
					count++;
					if (count == 3)
					{
						f21 << endl;
						count = 0;
					}
				}
			count = 0;
			for (int i = 1; i<nz + 1; i++)
				for (int k = 1; k<nx + 1; k++)
				{
					f21 << setprecision(15) << hwwma[nn][i][k] << "      ";
					count++;
					if (count == 3)
					{
						f21 << endl;
						count = 0;
					}
				}
			count = 0;
			for (int i = 1; i<nz + 1; i++)
				for (int k = 1; k<nx + 1; k++)
				{
					f21 << setprecision(15) << hpp[nn][i][k] << "      ";
					count++;
					if (count == 3)
					{
						f21 << endl;
						count = 0;
					}
				}
			count = 0;
			for (int i = 1; i<nz + 1; i++)
				for (int k = 1; k<nx + 1; k++)
				{
					f21 << setprecision(15) << hppt[nn][i][k] << "      ";
					count++;
					if (count == 3)
					{
						f21 << endl;
						count = 0;
					}
				}
			count = 0;
			for (int i = 1; i<nz + 1; i++)
				for (int k = 1; k<nx + 1; k++)
				{
					f21 << setprecision(15) << hht[nn][i][k] << "      ";
					count++;
					if (count == 3)
					{
						f21 << endl;
						count = 0;
					}
				}
			count = 0;
			for (int i = 1; i<nz + 1; i++)
				for (int k = 1; k<nx + 1; k++)
				{
					f21 << setprecision(15) << hmmiu[nn][i][k] << "      ";
					count++;
					if (count == 3)
					{
						f21 << endl;
						count = 0;
					}
				}
		}
	}
	if(myid==0)LogInfo::Log("Output "+Int_to_string(spa)+"%span-* files end");
	f21.close();

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


void CWriteFile::flow(int nitt)
{
	double  vf1, vf2, qin, qout, qinn, qoutt;
	double* qin0, *qout0;
	qin = 0;
	qout = 0;
	for (int nn = 1; nn<pdict->nt + 1; nn++)
		for (int j = 1; j<nz + 1; j++)
			for (int k = 1; k<ny + 1; k++)			
			{
				vf1 = -(pgrid->s2x[j][k][1] * pfield->q12[nn][j][k][1] + pgrid->s2y[j][k][1] * pfield->q13[nn][j][k][1] + pgrid->s2z[j][k][1] * pfield->q14[nn][j][k][1]);
				qin = qin + vf1;

				vf2 = -(pgrid->s2x[j][k][nx + 1] * pfield->q12[nn][j][k][nx] + pgrid->s2y[j][k][nx + 1] * pfield->q13[nn][j][k][nx] + pgrid->s2z[j][k][nx + 1] * pfield->q14[nn][j][k][nx]);
				qout = qout + vf2;
			}

	qin = qin / double(pdict->nt);
	qout = qout / double(pdict->nt);
	if (myid == 0)
	{

		Malloc(qin0,pdict->lbb);
		Malloc(qout0,pdict->lbb);

		qin0[0] = qin;
		qinn = 0;
		qout0[0] = qout;
		qoutt = 0;
		for (int j = 1; j<pdict->lbb; j++)
		{
			MPI_Recv(&qin0[j], 1, MPI_DOUBLE, j, 3, MPI_COMM_WORLD, &(status));
			MPI_Recv(&qout0[j], 1, MPI_DOUBLE, j, 4, MPI_COMM_WORLD, &(status));
		}

		for (int j = 0; j<pdict->lbb; j++)
		{
			qinn = qinn + qin0[j];
			qoutt = qoutt + qout0[j];
		}
		f1 << setw(5) << nitt;
		f1 << setw(15) << setprecision(6) << qinn;
		f2 << setw(5) << nitt;
		f2 << setw(15) << setprecision(6) << qoutt;
		Free(qin0);
		Free(qout0);

	}
	else
	{
		MPI_Send(&qin, 1, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
		MPI_Send(&qout, 1, MPI_DOUBLE, 0, 4, MPI_COMM_WORLD);
	}

}
void CWriteFile::output()
{
	int  spa;

	//pfield->mesh->x = pfield->mesh->xf;
	//pfield->mesh->y = pfield->mesh->yf;
	//pfield->mesh->z = pfield->mesh->zf;
	memcpy(pfield->mesh->x[0][0] ,pfield->mesh->xf[0][0],((pfield->mesh)->nzm1+1)*( (pfield->mesh)->nym1+1)*((pfield->mesh)->nxm1+1));
	memcpy(pfield->mesh->y[0][0] ,pfield->mesh->yf[0][0],((pfield->mesh)->nzm1+1)*( (pfield->mesh)->nym1+1)*((pfield->mesh)->nxm1+1));
	memcpy(pfield->mesh->z[0][0] ,pfield->mesh->zf[0][0],((pfield->mesh)->nzm1+1)*( (pfield->mesh)->nym1+1)*((pfield->mesh)->nxm1+1));

	inlout();
	lbout();

	spa = 30;
	span(spa);
	spa = 95;
	span(spa);
}








// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
