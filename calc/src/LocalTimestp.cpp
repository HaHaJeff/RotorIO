/*---------------------------------------------------------------------------*\


\*---------------------------------------------------------------------------*/


#include "LocalTimestp.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


void CLocalTimestp::viscosity(double t, double q6, double &cv, double &kc)
{
		
	double cvl, cvt, fv1, tem;
	cvl = (pdict->cvl0)*pow(t / (pdict->t0), 1.5)*((pdict->t0) + (pdict->ts)) / (t + (pdict->ts));
	tem = q6 / cvl;
	fv1 = 1.0 / (1.0 + pow((pdict->cv1) / tem, 3));
	cvt = q6*fv1;
	cv = cvl + cvt;
	kc = (pdict->cp)*(cvl / (pdict->prl) + cvt / (pdict->prt));	
}

void CLocalTimestp::InitTimeData()
{
	double tc, td, aa,cv,kc;
	double vx,vy,vz,wx,wy,wz;
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		// !******x							
		vx = 0.5*(pfield->pvx[k][j][i] + pfield->pvx[k][j][i + 1]);
		vy = 0.5*(pfield->pvy[k][j][i] + pfield->pvy[k][j][i + 1]);
		vz = 0.5*(pfield->pvz[k][j][i] + pfield->pvz[k][j][i + 1]);
		wx = vx;
		wy = vy + (pdict->rpm)*pgrid->zz02[k][j][i + 1];
		wz = vz - (pdict->rpm)*pgrid->yy02[k][j][i + 1];
		tc = abs(wx*pgrid->s2x[k][j][i + 1] + wy*pgrid->s2y[k][j][i + 1] + wz*pgrid->s2z[k][j][i + 1]);

		vx = 0.5*(pfield->pvx[k][j][i] + pfield->pvx[k][j + 1][i]);
		vy = 0.5*(pfield->pvy[k][j][i] + pfield->pvy[k][j + 1][i]);
		vz = 0.5*(pfield->pvz[k][j][i] + pfield->pvz[k][j + 1][i]);
		wx = vx;
		wy = vy + (pdict->rpm)*pgrid->zz03[k][j + 1][i];
		wz = vz - (pdict->rpm)*pgrid->yy03[k][j + 1][i];
		tc = abs(wx*pgrid->s3x[k][j + 1][i] + wy*pgrid->s3y[k][j + 1][i] + wz*pgrid->s3z[k][j + 1][i]) + tc;

		//******z							
		vx = 0.5*(pfield->pvx[k][j][i] + pfield->pvx[k + 1][j][i]);
		vy = 0.5*(pfield->pvy[k][j][i] + pfield->pvy[k + 1][j][i]);
		vz = 0.5*(pfield->pvz[k][j][i] + pfield->pvz[k + 1][j][i]);
		wx = vx;
		wy = vy + (pdict->rpm)*pgrid->zz01[k + 1][j][i];
		wz = vz - (pdict->rpm)*pgrid->yy01[k + 1][j][i];
		tc = abs(wx*pgrid->s1x[k + 1][j][i] + wy*pgrid->s1y[k + 1][j][i] + wz*pgrid->s1z[k + 1][j][i]) + tc;
		aa = sqrt(1.4*pfield->p[k][j][i] / pfield->q11[n][k][j][i]);

		tc = tc + aa*(sqrt(pow(pgrid->s1x[k + 1][j][i], 2) + pow(pgrid->s1y[k + 1][j][i], 2) + pow(pgrid->s1z[k + 1][j][i], 2)) +
					sqrt(pow(pgrid->s2x[k][j][i + 1], 2) + pow(pgrid->s2y[k][j][i + 1], 2) + pow(pgrid->s2z[k][j][i + 1], 2)) +
					sqrt(pow(pgrid->s3x[k][j + 1][i], 2) + pow(pgrid->s3y[k][j + 1][i], 2) + pow(pgrid->s3z[k][j + 1][i], 2)));

		td = pow(pgrid->s1x[k + 1][j][i], 2) + pow(pgrid->s1y[k + 1][j][i], 2) + pow(pgrid->s1z[k + 1][j][i], 2) +
				pow(pgrid->s2x[k][j][i + 1], 2) + pow(pgrid->s2y[k][j][i + 1], 2) + pow(pgrid->s2z[k][j][i + 1], 2) +
				pow(pgrid->s3x[k][j + 1][i], 2) + pow(pgrid->s3y[k][j + 1][i], 2) + pow(pgrid->s3z[k][j + 1][i], 2) +
				2.0*(abs(pgrid->s1x[k + 1][j][i] * pgrid->s2x[k][j][i + 1] + pgrid->s1y[k + 1][j][i] * pgrid->s2y[k][j][i + 1] + 
				pgrid->s1z[k + 1][j][i] * pgrid->s2z[k][j][i + 1]) +abs(pgrid->s1x[k + 1][j][i] * pgrid->s3x[k][j + 1][i] + 
				pgrid->s1y[k + 1][j][i] * pgrid->s3y[k][j + 1][i] + pgrid->s1z[k + 1][j][i] * pgrid->s3z[k][j + 1][i]) +
				abs(pgrid->s3x[k][j + 1][i] * pgrid->s2x[k][j][i + 1] + pgrid->s3y[k][j + 1][i] * pgrid->s2y[k][j][i + 1] + 
				pgrid->s3z[k][j + 1][i] * pgrid->s2z[k][j][i + 1]));

		viscosity(pfield->t[k][j][i], pfield->q16[n][k][j][i], cv,kc);

		td = 8.0*cv*td / pfield->q11[n][k][j][i] / pgrid->vv[k][j][i];
		time[k][j][i] = tc + td;
		sri[k][j][i] = tc + td;
		srj[k][j][i] = tc + td;
		srk[k][j][i] = tc + td;
	}
}

void CLocalTimestp::ComputeTimeData()
{

	double*** sri1,***srj1,***srk1;
	Malloc(sri1,nz+2,ny+2,nx+2);
	Malloc(srj1,nz+2,ny+2,nx+2);
	Malloc(srk1,nz+2,ny+2,nx+2);
	
	
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		sri1[k][j][i] = sri[k][j][i] * (1 + sqrt(srj[k][j][i] / sri[k][j][i]) + sqrt(srk[k][j][i] / sri[k][j][i]));
		srj1[k][j][i] = srj[k][j][i] * (1 + sqrt(sri[k][j][i] / srj[k][j][i]) + sqrt(srk[k][j][i] / srj[k][j][i]));
		srk1[k][j][i] = srk[k][j][i] * (1 + sqrt(sri[k][j][i] / srk[k][j][i]) + sqrt(srj[k][j][i] / srk[k][j][i]));
	}
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	{
		sri1[k][j][0] = sri1[k][j][1];
		srj1[k][j][0] = srj1[k][j][1];
		srk1[k][j][0] = srk1[k][j][1];

		sri1[k][j][nx + 1] = sri1[k][j][nx];
		srj1[k][j][nx + 1] = srj1[k][j][nx];
		srk1[k][j][nx + 1] = srk1[k][j][nx];
	}

	for (int k = 1; k<nz + 1; k++)
	for (int i = 1; i<nx + 1; i++)
	{
		sri1[k][0][i] = sri1[k][1][i];
		srj1[k][0][i] = srj1[k][1][i];
		srk1[k][0][i] = srk1[k][1][i];

		sri1[k][ny + 1][i] = sri1[k][ny][i];
		srj1[k][ny + 1][i] = srj1[k][ny][i];
		srk1[k][ny + 1][i] = srk1[k][ny][i];
	}

	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		sri1[0][j][i] = sri1[1][j][i];
		srj1[0][j][i] = srj1[1][j][i];
		srk1[0][j][i] = srk1[1][j][i];

		sri1[nz + 1][j][i] = sri1[nz][j][i];
		srj1[nz + 1][j][i] = srj1[nz][j][i];
		srk1[nz + 1][j][i] = srk1[nz][j][i];
	}

	for (int k = 0; k<nz + 2; k++)
	for (int j = 0; j<ny + 2; j++)
	for (int i = 0; i<nx + 2; i++)
	{
		sri[k][j][i] = sri1[k][j][i];
		srj[k][j][i] = srj1[k][j][i];
		srk[k][j][i] = srk1[k][j][i];
	}


	Free(sri1);
	Free(srj1);
	Free(srk1);	
}

bool CLocalTimestp::SetLocalTimestep()
{
	InitTimeData();
	ComputeTimeData();
	return true;
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



