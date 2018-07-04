/*---------------------------------------------------------------------------*\


\*---------------------------------------------------------------------------*/


#include "Convective.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void CConvective::ResetConvectiveData()
{
	for (int k = 0; k<nz + 2; k++)
	for (int j = 0; j<ny + 2; j++)
	for (int i = 0; i<nx + 2; i++)
	{
		qc1[k][j][i] = 0;
		qc2[k][j][i] = 0;
		qc3[k][j][i] = 0;
		qc4[k][j][i] = 0;
		qc5[k][j][i] = 0;
	}
}
void CConvective::ComputeConvectiveData_i()
{
	double flu1, flu2, flu3, flu4, flu5, vf, rf;
	double vx, vy, vz, dim, pp, en;
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 2; i++)
	{
		vx = 0.50*(pfield->pvx[k][j][i] + pfield->pvx[k][j][i - 1]);
		vy = 0.50*(pfield->pvy[k][j][i] + pfield->pvy[k][j][i - 1]);
		vz = 0.50*(pfield->pvz[k][j][i] + pfield->pvz[k][j][i - 1]);
		vf = vx*pgrid->s2x[k][j][i] + vy*pgrid->s2y[k][j][i] + vz*pgrid->s2z[k][j][i];
		rf = pdict->rpm*pgrid->zz02[k][j][i] * pgrid->s2y[k][j][i] - pdict->rpm*pgrid->yy02[k][j][i] * pgrid->s2z[k][j][i];
		vf = vf + rf;

		dim = 0.50*(pfield->q11[n][k][j][i] + pfield->q11[n][k][j][i - 1]);
		pp = 0.50*(pfield->p[k][j][i] + pfield->p[k][j][i - 1]);
		en = 0.50*(pfield->q15[n][k][j][i] + pfield->q15[n][k][j][i - 1]);

		flu1 = dim*vf;
		flu2 = flu1*vx + pp*pgrid->s2x[k][j][i];
		flu3 = flu1*vy + pp*pgrid->s2y[k][j][i];
		flu4 = flu1*vz + pp*pgrid->s2z[k][j][i];
		flu5 = (en + pp)*vf - pp*rf;

		qc1[k][j][i] = qc1[k][j][i] + flu1;
		qc2[k][j][i] = qc2[k][j][i] + flu2;
		qc3[k][j][i] = qc3[k][j][i] + flu3;
		qc4[k][j][i] = qc4[k][j][i] + flu4;
		qc5[k][j][i] = qc5[k][j][i] + flu5;

		qc1[k][j][i - 1] = qc1[k][j][i - 1] - flu1;
		qc2[k][j][i - 1] = qc2[k][j][i - 1] - flu2;
		qc3[k][j][i - 1] = qc3[k][j][i - 1] - flu3;
		qc4[k][j][i - 1] = qc4[k][j][i - 1] - flu4;
		qc5[k][j][i - 1] = qc5[k][j][i - 1] - flu5;
	}
}
void CConvective::ComputeConvectiveData_j()
{
	double flu1, flu2, flu3, flu4, flu5, vf, rf;
	double vx, vy, vz, dim, pp, en;
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 2; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		vx = 0.50*(pfield->pvx[k][j][i] + pfield->pvx[k][j - 1][i]);
		vy = 0.50*(pfield->pvy[k][j][i] + pfield->pvy[k][j - 1][i]);
		vz = 0.50*(pfield->pvz[k][j][i] + pfield->pvz[k][j - 1][i]);
		vf = vx*pgrid->s3x[k][j][i] + vy*pgrid->s3y[k][j][i] + vz*pgrid->s3z[k][j][i];
		rf = pdict->rpm*pgrid->zz03[k][j][i] * pgrid->s3y[k][j][i] - pdict->rpm*pgrid->yy03[k][j][i] * pgrid->s3z[k][j][i];
		vf = vf + rf;

		dim = 0.50*(pfield->q11[n][k][j][i] + pfield->q11[n][k][j - 1][i]);
		pp = 0.50*(pfield->p[k][j][i] + pfield->p[k][j - 1][i]);
		en = 0.50*(pfield->q15[n][k][j][i] + pfield->q15[n][k][j - 1][i]);

		flu1 = dim*vf;
		flu2 = flu1*vx + pp*pgrid->s3x[k][j][i];
		flu3 = flu1*vy + pp*pgrid->s3y[k][j][i];
		flu4 = flu1*vz + pp*pgrid->s3z[k][j][i];
		flu5 = (en + pp)*vf - pp*rf;

		qc1[k][j][i] = qc1[k][j][i] + flu1;
		qc2[k][j][i] = qc2[k][j][i] + flu2;
		qc3[k][j][i] = qc3[k][j][i] + flu3;
		qc4[k][j][i] = qc4[k][j][i] + flu4;
		qc5[k][j][i] = qc5[k][j][i] + flu5;

		qc1[k][j - 1][i] = qc1[k][j - 1][i] - flu1;
		qc2[k][j - 1][i] = qc2[k][j - 1][i] - flu2;
		qc3[k][j - 1][i] = qc3[k][j - 1][i] - flu3;
		qc4[k][j - 1][i] = qc4[k][j - 1][i] - flu4;
		qc5[k][j - 1][i] = qc5[k][j - 1][i] - flu5;
	}
}
void CConvective::ComputeConvectiveData_k()
{
	double flu1, flu2, flu3, flu4, flu5, vf, rf;
	double vx, vy, vz, dim, pp, en;
	for (int k = 1; k<nz + 2; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		vx = 0.50*(pfield->pvx[k][j][i] + pfield->pvx[k - 1][j][i]);
		vy = 0.50*(pfield->pvy[k][j][i] + pfield->pvy[k - 1][j][i]);
		vz = 0.50*(pfield->pvz[k][j][i] + pfield->pvz[k - 1][j][i]);
		vf = vx*pgrid->s1x[k][j][i] + vy*pgrid->s1y[k][j][i] + vz*pgrid->s1z[k][j][i];
		rf = pdict->rpm*pgrid->zz01[k][j][i] * pgrid->s1y[k][j][i] - pdict->rpm*pgrid->yy01[k][j][i] * pgrid->s1z[k][j][i];
		vf = vf + rf;

		dim = 0.50*(pfield->q11[n][k][j][i] + pfield->q11[n][k - 1][j][i]);
		pp = 0.50*(pfield->p[k][j][i] + pfield->p[k - 1][j][i]);
		en = 0.50*(pfield->q15[n][k][j][i] + pfield->q15[n][k - 1][j][i]);

		flu1 = dim*vf;
		flu2 = flu1*vx + pp*pgrid->s1x[k][j][i];
		flu3 = flu1*vy + pp*pgrid->s1y[k][j][i];
		flu4 = flu1*vz + pp*pgrid->s1z[k][j][i];
		flu5 = (en + pp)*vf - pp*rf;

		qc1[k][j][i] = qc1[k][j][i] + flu1;
		qc2[k][j][i] = qc2[k][j][i] + flu2;
		qc3[k][j][i] = qc3[k][j][i] + flu3;
		qc4[k][j][i] = qc4[k][j][i] + flu4;
		qc5[k][j][i] = qc5[k][j][i] + flu5;

		qc1[k - 1][j][i] = qc1[k - 1][j][i] - flu1;
		qc2[k - 1][j][i] = qc2[k - 1][j][i] - flu2;
		qc3[k - 1][j][i] = qc3[k - 1][j][i] - flu3;
		qc4[k - 1][j][i] = qc4[k - 1][j][i] - flu4;
		qc5[k - 1][j][i] = qc5[k - 1][j][i] - flu5;
	}
}
bool CConvective::ComputeConvectiveData()
{

	ResetConvectiveData();
	ComputeConvectiveData_i();
	ComputeConvectiveData_j();
	ComputeConvectiveData_k();
	return true;
}
void CConvective::ResetConvectiveData_SA()
{
	for (int k = 0; k<nz + 2; k++)
		for (int j = 0; j<ny + 2; j++)
			for (int i = 0; i<nx + 2; i++)
			{
				qc6[k][j][i] = 0;
			}
}
void CConvective::ComputeConvectiveData_SAi()
{
	double flu6, tur, vf, rf;
	double vx, vy, vz, dim;
	for (int k = 1; k<nz + 1; k++)
		for (int j = 1; j<ny + 1; j++)
			for (int i = 1; i<nx + 2; i++)
			{
				vx = 0.50*(pfield->pvx[k][j][i] + pfield->pvx[k][j][i - 1]);
				vy = 0.50*(pfield->pvy[k][j][i] + pfield->pvy[k][j][i - 1]);
				vz = 0.50*(pfield->pvz[k][j][i] + pfield->pvz[k][j][i - 1]);
				vf = vx*pgrid->s2x[k][j][i] + vy*pgrid->s2y[k][j][i] + vz*pgrid->s2z[k][j][i];
				rf = pdict->rpm*pgrid->zz02[k][j][i] * pgrid->s2y[k][j][i] - pdict->rpm*pgrid->yy02[k][j][i] * pgrid->s2z[k][j][i];
				vf = vf + rf;
				tur = 0.50*(pfield->q16[n][k][j][i] + pfield->q16[n][k][j][i - 1]);

				flu6 = tur*vf;
				qc6[k][j][i] = qc6[k][j][i] + flu6;
				qc6[k][j][i - 1] = qc6[k][j][i - 1] - flu6;
			}
}
void CConvective::ComputeConvectiveData_SAj()
{
	double flu6, tur, vf, rf;
	double vx, vy, vz, dim;
	for (int k = 1; k<nz + 1; k++)
		for (int j = 1; j<ny + 2; j++)
			for (int i = 1; i<nx + 1; i++)
			{
				vx = 0.50*(pfield->pvx[k][j][i] + pfield->pvx[k][j - 1][i]);
				vy = 0.50*(pfield->pvy[k][j][i] + pfield->pvy[k][j - 1][i]);
				vz = 0.50*(pfield->pvz[k][j][i] + pfield->pvz[k][j - 1][i]);
				vf = vx*pgrid->s3x[k][j][i] + vy*pgrid->s3y[k][j][i] + vz*pgrid->s3z[k][j][i];
				rf = pdict->rpm*pgrid->zz03[k][j][i] * pgrid->s3y[k][j][i] - pdict->rpm*pgrid->yy03[k][j][i] * pgrid->s3z[k][j][i];
				vf = vf + rf;

				tur = 0.50*(pfield->q16[n][k][j][i] + pfield->q16[n][k][j - 1][i]);
				flu6 = tur*vf;
				qc6[k][j][i] = qc6[k][j][i] + flu6;
				qc6[k][j - 1][i] = qc6[k][j - 1][i] - flu6;
			}
}
void CConvective::ComputeConvectiveData_SAk()
{
	double flu6, tur, vf, rf;
	double vx, vy, vz, dim;
	for (int k = 1; k<nz + 2; k++)
		for (int j = 1; j<ny + 1; j++)
			for (int i = 1; i<nx + 1; i++)
			{
				vx = 0.50*(pfield->pvx[k][j][i] + pfield->pvx[k - 1][j][i]);
				vy = 0.50*(pfield->pvy[k][j][i] + pfield->pvy[k - 1][j][i]);
				vz = 0.50*(pfield->pvz[k][j][i] + pfield->pvz[k - 1][j][i]);
				vf = vx*pgrid->s1x[k][j][i] + vy*pgrid->s1y[k][j][i] + vz*pgrid->s1z[k][j][i];
				rf = pdict->rpm*pgrid->zz01[k][j][i] * pgrid->s1y[k][j][i] - pdict->rpm*pgrid->yy01[k][j][i] * pgrid->s1z[k][j][i];
				vf = vf + rf;

				tur = 0.50*(pfield->q16[n][k][j][i] + pfield->q16[n][k - 1][j][i]);
				flu6 = tur*vf;
				qc6[k][j][i] = qc6[k][j][i] + flu6;
				qc6[k - 1][j][i] = qc6[k - 1][j][i] - flu6;
			}
}
bool CConvective::ComputeConvectiveData_SA()
{

	ResetConvectiveData_SA();
	ComputeConvectiveData_SAi();
	ComputeConvectiveData_SAj();
	ComputeConvectiveData_SAk();
	return true;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
