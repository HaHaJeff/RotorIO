/*---------------------------------------------------------------------------*\


\*---------------------------------------------------------------------------*/


#include "RSmoothing.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


void CRSmoothing::tdma(double*& a, double*& b, double*& c, double*& d, double*& x, int n)
{

	double* u,*l,*y;
	Malloc(u,n + 1);
	Malloc(l,n + 1);
	Malloc(y,n + 1);


	u[1] = b[1];
	y[1] = d[1];
	for (int i = 2; i<n + 1; i++)
	{
		l[i] = a[i] / u[i - 1];
		u[i] = b[i] - l[i] * c[i - 1];
		y[i] = d[i] - l[i] * y[i - 1];
	}
	x[n] = y[n] / u[n];
	for (int i = n - 1; i>0; i--)
	{
		x[i] = (y[i] - c[i] * x[i + 1]) / u[i];
	}
	Free(u);
	Free(l);
	Free(y);
}
void CRSmoothing::InitResidual()
{
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		rr1[k][j][i] = 0;
		rr2[k][j][i] = 0;
		rr3[k][j][i] = 0;
		rr4[k][j][i] = 0;
		rr5[k][j][i] = 0;
	}
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		rr1[k][j][i] = qp1[k][j][i] - qc1[k][j][i] + av1[k][j][i] - ts1[k][j][i] * pgrid->vv[k][j][i];
		rr2[k][j][i] = qp2[k][j][i] - qc2[k][j][i] + av2[k][j][i] - ts2[k][j][i] * pgrid->vv[k][j][i] + qv2[k][j][i];
		rr3[k][j][i] = qp3[k][j][i] - qc3[k][j][i] + av3[k][j][i] - ts3[k][j][i] * pgrid->vv[k][j][i] + qv3[k][j][i];
		rr4[k][j][i] = qp4[k][j][i] - qc4[k][j][i] + av4[k][j][i] - ts4[k][j][i] * pgrid->vv[k][j][i] + qv4[k][j][i];
		rr5[k][j][i] = qp5[k][j][i] - qc5[k][j][i] + av5[k][j][i] - ts5[k][j][i] * pgrid->vv[k][j][i] + qv5[k][j][i];
	}
}
void CRSmoothing::SmoothResidual(double ta)
{

	double* ax,*ay,*az,*bx,*by,*bz;
	double*  c1, *c2, *c3, *c4, *c5;
	Malloc(ax,nx + 1);
	Malloc(bx,nx + 1);
	Malloc(ay,ny + 1);
	Malloc(by,ny + 1);
	Malloc(az,nz + 1);
	Malloc(bz,nz + 1);

	for (int i = 1; i<nx + 1; i++)
	{
		ax[i] = -ta;
		bx[i] = 1.0 + 2.0 * ta;
	}
	for (int i = 1; i<ny + 1; i++)
	{
		ay[i] = -ta;
		by[i] = 1.0 + 2.0 * ta;
	}
	for (int i = 1; i<nz + 1; i++)
	{
		az[i] = -ta;
		bz[i] = 1.0 + 2.0 * ta;
	}
	//

	Malloc(c1,nx + 1);
	Malloc(c2,nx + 1);
	Malloc(c3,nx + 1);
	Malloc(c4,nx + 1);
	Malloc(c5,nx + 1);
	for (int k = 1; k<nz + 1; k++)
		for (int j = 1; j<ny + 1; j++)
		{
			for (int i = 1; i<nx + 1; i++)
			{
				c1[i] = py1[k][j][i];
				c2[i] = py2[k][j][i];
				c3[i] = py3[k][j][i];
				c4[i] = py4[k][j][i];
				c5[i] = py5[k][j][i];
			}
			tdma(ax, bx, ax, c1, c1, nx);
			tdma(ax, bx, ax, c2, c2, nx);
			tdma(ax, bx, ax, c3, c3, nx);
			tdma(ax, bx, ax, c4, c4, nx);
			tdma(ax, bx, ax, c5, c5, nx);
			for (int i = 1; i<nx + 1; i++)
			{
				py1[k][j][i] = c1[i];
				py2[k][j][i] = c2[i];
				py3[k][j][i] = c3[i];
				py4[k][j][i] = c4[i];
				py5[k][j][i] = c5[i];
			}
		}
	Free(c1);
	Free(c2);
	Free(c3);
	Free(c4);
	Free(c5);

	Malloc(c1,ny + 1);
	Malloc(c2,ny + 1);
	Malloc(c3,ny + 1);
	Malloc(c4,ny + 1);
	Malloc(c5,ny + 1);
	for (int k = 1; k<nz + 1; k++)
		for (int i = 1; i<nx + 1; i++)
		{
			for (int j = 1; j<ny + 1; j++)
			{
				c1[j] = py1[k][j][i];
				c2[j] = py2[k][j][i];
				c3[j] = py3[k][j][i];
				c4[j] = py4[k][j][i];
				c5[j] = py5[k][j][i];
			}
			tdma(ay, by, ay, c1, c1, ny);
			tdma(ay, by, ay, c2, c2, ny);
			tdma(ay, by, ay, c3, c3, ny);
			tdma(ay, by, ay, c4, c4, ny);
			tdma(ay, by, ay, c5, c5, ny);
			for (int j = 1; j<ny + 1; j++)
			{
				py1[k][j][i] = c1[j];
				py2[k][j][i] = c2[j];
				py3[k][j][i] = c3[j];
				py4[k][j][i] = c4[j];
				py5[k][j][i] = c5[j];
			}
		}
	Free(c1);
	Free(c2);
	Free(c3);
	Free(c4);
	Free(c5);

	Malloc(c1,nz + 1);
	Malloc(c2,nz + 1);
	Malloc(c3,nz + 1);
	Malloc(c4,nz + 1);
	Malloc(c5,nz + 1);
	
	for (int j = 1; j<ny + 1; j++)
		for (int i = 1; i<nx + 1; i++)
		{
			for (int k = 1; k<nz + 1; k++)
			{
				c1[k] = py1[k][j][i];
				c2[k] = py2[k][j][i];
				c3[k] = py3[k][j][i];
				c4[k] = py4[k][j][i];
				c5[k] = py5[k][j][i];
			}
			tdma(az, bz, az, c1, c1, nz);
			tdma(az, bz, az, c2, c2, nz);
			tdma(az, bz, az, c3, c3, nz);
			tdma(az, bz, az, c4, c4, nz);
			tdma(az, bz, az, c5, c5, nz);
			for (int k = 1; k<nz + 1; k++)
			{
				py1[k][j][i] = c1[k];
				py2[k][j][i] = c2[k];
				py3[k][j][i] = c3[k];
				py4[k][j][i] = c4[k];
				py5[k][j][i] = c5[k];
			}
		}
	Free(c1);
	Free(c2);
	Free(c3);
	Free(c4);
	Free(c5);
	Free(ax);
	Free(bx);
	Free(ay);
	Free(by);
	Free(az);
	Free(bz);	
}
void CRSmoothing::SmoothResidual_SA(double ta)
{
	
	double* ax,*ay,*az,*bx,*by,*bz;
	double*  c6;
	Malloc(ax,nx + 1);
	Malloc(bx,nx + 1);
	Malloc(ay,ny + 1);
	Malloc(by,ny + 1);
	Malloc(az,nz + 1);
	Malloc(bz,nz + 1);
	Malloc(c6,nx + 1);
	for (int i = 0; i<nx + 1; i++)
	{
		ax[i] = -ta;
		bx[i] = 1.0 + 2.0 * ta;
	}
	for (int i = 0; i<ny + 1; i++)
	{
		ay[i] = -ta;
		by[i] = 1.0 + 2.0 * ta;
	}
	for (int i = 0; i<nz + 1; i++)
	{
		az[i] = -ta;
		bz[i] = 1.0 + 2.0 * ta;
	}


	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	{
		for (int i = 1; i<nx + 1; i++)
		{
			c6[i] = py6[k][j][i];
		}
		tdma(ax, bx, ax, c6, c6, nx);
		for (int i = 1; i<nx + 1; i++)
		{
			py6[k][j][i] = c6[i];
		}
	}
	Free(c6);
	Malloc(c6,ny + 1);

	for (int k = 1; k<nz + 1; k++)
	for (int i = 1; i<nx + 1; i++)
	{
		for (int j = 1; j<ny + 1; j++)
		{
			c6[j] = py6[k][j][i];
		}
		tdma(ay, by, ay, c6, c6, ny);
		for (int j = 1; j<ny + 1; j++)
		{
			py6[k][j][i] = c6[j];
		}
	}
	Free(c6);
	Malloc(c6,nz + 1);

	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		for (int k = 1; k<nz + 1; k++)
		{
			c6[k] = py6[k][j][i];
		}
		tdma(az, bz, az, c6, c6, nz);
		for (int k = 1; k<nz + 1; k++)
		{
			py6[k][j][i] = c6[k];
		}
	}
	Free(c6);
	Free(ax);
	Free(bx);
	Free(ay);
	Free(by);
	Free(az);
	Free(bz);
}
void CRSmoothing::UpdateFieldValue(double ta, double timl,bool smoothing)
{
	double dtime;
	InitResidual();
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		py1[k][j][i] = 0;
		py2[k][j][i] = 0;
		py3[k][j][i] = 0;
		py4[k][j][i] = 0;
		py5[k][j][i] = 0;
	}
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		dtime = timl / time[k][j][i];
		py1[k][j][i] = dtime*rr1[k][j][i];
		py2[k][j][i] = dtime*rr2[k][j][i];
		py3[k][j][i] = dtime*rr3[k][j][i];
		py4[k][j][i] = dtime*rr4[k][j][i];
		py5[k][j][i] = dtime*rr5[k][j][i];
	}
	if (smoothing)
		SmoothResidual(ta);
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		pfield->q11[n][k][j][i] = q01[k][j][i] + py1[k][j][i];
		pfield->q12[n][k][j][i] = q02[k][j][i] + py2[k][j][i];
		pfield->q13[n][k][j][i] = q03[k][j][i] + py3[k][j][i];
		pfield->q14[n][k][j][i] = q04[k][j][i] + py4[k][j][i];
		pfield->q15[n][k][j][i] = q05[k][j][i] + py5[k][j][i];
	}
}
void CRSmoothing::UpdateFieldValue_SA(double ta,double timl,bool smoothing)
{
	double dtime;
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		py6[k][j][i] = 0;
	}

	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		dtime = timl / time[k][j][i];
		py6[k][j][i] = dtime*(-qc6[k][j][i] + av6[k][j][i] - ts6[k][j][i] * pgrid->vv[k][j][i] + qv6[k][j][i]);
	}
	if (smoothing)
		SmoothResidual_SA(ta);
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		pfield->q16[n][k][j][i] = q06[k][j][i] + py6[k][j][i];
	}
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
