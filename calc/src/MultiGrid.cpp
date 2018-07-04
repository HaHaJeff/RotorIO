/*---------------------------------------------------------------------------*\


\*---------------------------------------------------------------------------*/


#include "MultiGrid.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void CMultiGrid::ComputeSurface_CenterPoint(int nng)
{
	int nx = ((mesh->nnx))[nng];
	int ny = ((mesh->nny))[nng];
	int nz = ((mesh->nnz))[nng];
	double x24, y24, z24, x31, y31, z31;	
	for (int k = 1; k<nz + 2; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		x24 = ((mesh->x))[k][j + 1][i] - ((mesh->x))[k][j][i + 1];
		y24 = ((mesh->y))[k][j + 1][i] - ((mesh->y))[k][j][i + 1];
		z24 = ((mesh->z))[k][j + 1][i] - ((mesh->z))[k][j][i + 1];
		x31 = ((mesh->x))[k][j + 1][i + 1] - ((mesh->x))[k][j][i];
		y31 = ((mesh->y))[k][j + 1][i + 1] - ((mesh->y))[k][j][i];
		z31 = ((mesh->z))[k][j + 1][i + 1] - ((mesh->z))[k][j][i];
		s1xn[nng][k][j][i] = 0.5*(y24*z31 - z24*y31);
		s1yn[nng][k][j][i] = 0.5*(z24*x31 - x24*z31);
		s1zn[nng][k][j][i] = 0.5*(x24*y31 - y24*x31);
		xx1[nng][k][j][i] = 0.25*(((mesh->x))[k][j][i] + ((mesh->x))[k][j + 1][i] + ((mesh->x))[k][j][i + 1] + ((mesh->x))[k][j + 1][i + 1]);
		yy1[nng][k][j][i] = 0.25*(((mesh->y))[k][j][i] + ((mesh->y))[k][j + 1][i] + ((mesh->y))[k][j][i + 1] + ((mesh->y))[k][j + 1][i + 1]);
		zz1[nng][k][j][i] = 0.25*(((mesh->z))[k][j][i] + ((mesh->z))[k][j + 1][i] + ((mesh->z))[k][j][i + 1] + ((mesh->z))[k][j + 1][i + 1]);
	}
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 2; i++)
	{
		x24 = ((mesh->x))[k + 1][j + 1][i] - ((mesh->x))[k][j][i];
		y24 = ((mesh->y))[k + 1][j + 1][i] - ((mesh->y))[k][j][i];
		z24 = ((mesh->z))[k + 1][j + 1][i] - ((mesh->z))[k][j][i];
		x31 = ((mesh->x))[k][j + 1][i] - ((mesh->x))[k + 1][j][i];
		y31 = ((mesh->y))[k][j + 1][i] - ((mesh->y))[k + 1][j][i];
		z31 = ((mesh->z))[k][j + 1][i] - ((mesh->z))[k + 1][j][i];
		s2xn[nng][k][j][i] = 0.5*(y24*z31 - z24*y31);
		s2yn[nng][k][j][i] = 0.5*(z24*x31 - x24*z31);
		s2zn[nng][k][j][i] = 0.5*(x24*y31 - y24*x31);
		xx2[nng][k][j][i] = 0.25*(((mesh->x))[k][j][i] + ((mesh->x))[k][j + 1][i] + ((mesh->x))[k+1][j][i] + ((mesh->x))[k + 1][j + 1][i]);
		yy2[nng][k][j][i] = 0.25*(((mesh->y))[k][j][i] + ((mesh->y))[k][j + 1][i] + ((mesh->y))[k+1][j][i] + ((mesh->y))[k + 1][j + 1][i]);
		zz2[nng][k][j][i] = 0.25*(((mesh->z))[k][j][i] + ((mesh->z))[k][j + 1][i] + ((mesh->z))[k+1][j][i] + ((mesh->z))[k + 1][j + 1][i]);
	}
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 2; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		x24 = ((mesh->x))[k + 1][j][i + 1] - ((mesh->x))[k][j][i];
		y24 = ((mesh->y))[k + 1][j][i + 1] - ((mesh->y))[k][j][i];
		z24 = ((mesh->z))[k + 1][j][i + 1] - ((mesh->z))[k][j][i];
		x31 = ((mesh->x))[k + 1][j][i] - ((mesh->x))[k][j][i + 1];
		y31 = ((mesh->y))[k + 1][j][i] - ((mesh->y))[k][j][i + 1];
		z31 = ((mesh->z))[k + 1][j][i] - ((mesh->z))[k][j][i + 1];
		s3xn[nng][k][j][i] = 0.5*(y24*z31 - z24*y31);
		s3yn[nng][k][j][i] = 0.5*(z24*x31 - x24*z31);
		s3zn[nng][k][j][i] = 0.5*(x24*y31 - y24*x31);
		xx3[nng][k][j][i] = 0.25*(((mesh->x))[k][j][i] + ((mesh->x))[k][j][i + 1] + ((mesh->x))[k + 1][j][i] + ((mesh->x))[k + 1][j][i + 1]);
		yy3[nng][k][j][i] = 0.25*(((mesh->y))[k][j][i] + ((mesh->y))[k][j][i + 1] + ((mesh->y))[k + 1][j][i] + ((mesh->y))[k + 1][j][i + 1]);
		zz3[nng][k][j][i] = 0.25*(((mesh->z))[k][j][i] + ((mesh->z))[k][j][i + 1] + ((mesh->z))[k + 1][j][i] + ((mesh->z))[k + 1][j][i + 1]);
	}
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		x24 = ((mesh->x))[k + 1][j + 1][i + 1] - ((mesh->x))[k][j][i];
		y24 = ((mesh->y))[k + 1][j + 1][i + 1] - ((mesh->y))[k][j][i];
		z24 = ((mesh->z))[k + 1][j + 1][i + 1] - ((mesh->z))[k][j][i];
		vvn[nng][k][j][i] = -(x24*(s1xn[nng][k][j][i] + s2xn[nng][k][j][i] + s3xn[nng][k][j][i]) +
					y24*(s1yn[nng][k][j][i] + s2yn[nng][k][j][i] + s3yn[nng][k][j][i]) +
					z24*(s1zn[nng][k][j][i] + s2zn[nng][k][j][i] + s3zn[nng][k][j][i])) / 3.0;
				
		xx[nng][k][j][i] = 0.125*(((mesh->x))[k][j][i] + ((mesh->x))[k][j][i + 1] + ((mesh->x))[k + 1][j][i] + ((mesh->x))[k + 1][j][i + 1] +
					((mesh->x))[k][j + 1][i] + ((mesh->x))[k][j + 1][i + 1] + ((mesh->x))[k + 1][j+1][i] + ((mesh->x))[k + 1][j + 1][i + 1]);
				
		yy[nng][k][j][i] = 0.125*(((mesh->y))[k][j][i] + ((mesh->y))[k][j][i + 1] + ((mesh->y))[k + 1][j][i] + ((mesh->y))[k + 1][j][i + 1] +
					((mesh->y))[k][j + 1][i] + ((mesh->y))[k][j + 1][i + 1] + ((mesh->y))[k + 1][j+1][i] + ((mesh->y))[k + 1][j + 1][i + 1]);
				
		zz[nng][k][j][i] = 0.125*(((mesh->z))[k][j][i] + ((mesh->z))[k][j][i + 1] + ((mesh->z))[k + 1][j][i] + ((mesh->z))[k + 1][j][i + 1] +
					((mesh->z))[k][j + 1][i] + ((mesh->z))[k][j + 1][i + 1] + ((mesh->z))[k + 1][j+1][i] + ((mesh->z))[k + 1][j + 1][i + 1]);
	}
}
void CMultiGrid::V_interpolationPrepare(int nng)
{
	double val1, val2, val3, val4, val5, val6, val7, val8, val;
	int ifine, jfine, kfine,nx,ny,nz;
	for (int ign = nng - 1; ign >= 1; ign--)
	{
		nx = ((mesh->nnx))[ign];
		ny = ((mesh->nny))[ign];
		nz = ((mesh->nnz))[ign];

		for (int k = 1; k<nz + 1; k++)
		for (int j = 1; j<ny + 1; j++)
		for (int i = 1; i<nx + 2; i++)
		{
			kfine = 2 * k - 1;
			jfine = 2 * j - 1;
			ifine = 2 * i - 1;
			
			val1 = sqrt(pow(s2xn[ign + 1][kfine][jfine][ifine], 2) + pow(s2yn[ign + 1][kfine][jfine][ifine], 2) + pow(s2zn[ign + 1][kfine][jfine][ifine], 2));
			val2 = sqrt(pow(s2xn[ign + 1][kfine][jfine + 1][ifine], 2) + pow(s2yn[ign + 1][kfine][jfine + 1][ifine], 2) + pow(s2zn[ign + 1][kfine][jfine + 1][ifine], 2));
			val3 = sqrt(pow(s2xn[ign + 1][kfine + 1][jfine][ifine], 2) + pow(s2yn[ign + 1][kfine + 1][jfine][ifine], 2) + pow(s2zn[ign + 1][kfine + 1][jfine][ifine], 2));
			val4 = sqrt(pow(s2xn[ign + 1][kfine + 1][jfine + 1][ifine], 2) + pow(s2yn[ign + 1][kfine + 1][jfine + 1][ifine], 2) + pow(s2zn[ign + 1][kfine + 1][jfine + 1][ifine], 2));
			val = val1 + val2 + val3 + val4;
			xx2[ign][k][j][i] = (xx2[ign + 1][kfine][jfine][ifine] * val1 + xx2[ign + 1][kfine][jfine + 1][ifine] * val2
						+ xx2[ign + 1][kfine + 1][jfine][ifine] * val3 + xx2[ign + 1][kfine + 1][jfine + 1][ifine] * val4) / val;
						
			yy2[ign][k][j][i] = (yy2[ign + 1][kfine][jfine][ifine] * val1 + yy2[ign + 1][kfine][jfine + 1][ifine] * val2
						+ yy2[ign + 1][kfine + 1][jfine][ifine] * val3 + yy2[ign + 1][kfine + 1][jfine + 1][ifine] * val4) / val;
						
			zz2[ign][k][j][i] = (zz2[ign + 1][kfine][jfine][ifine] * val1 + zz2[ign + 1][kfine][jfine + 1][ifine] * val2
						+ zz2[ign + 1][kfine + 1][jfine][ifine] * val3 + zz2[ign + 1][kfine + 1][jfine + 1][ifine] * val4) / val;
						
						
			s2xn[ign][k][j][i] = s2xn[ign + 1][kfine][jfine][ifine] + s2xn[ign + 1][kfine][jfine + 1][ifine]
						+ s2xn[ign + 1][kfine + 1][jfine][ifine] + s2xn[ign + 1][kfine + 1][jfine + 1][ifine];
						
			s2yn[ign][k][j][i] = s2yn[ign + 1][kfine][jfine][ifine] + s2yn[ign + 1][kfine][jfine + 1][ifine]
						+ s2yn[ign + 1][kfine + 1][jfine][ifine] + s2yn[ign + 1][kfine + 1][jfine + 1][ifine];
						
			s2zn[ign][k][j][i] = s2zn[ign + 1][kfine][jfine][ifine] + s2zn[ign + 1][kfine][jfine + 1][ifine]
						+ s2zn[ign + 1][kfine + 1][jfine][ifine] + s2zn[ign + 1][kfine + 1][jfine + 1][ifine];
						
		}

					
		for (int k = 1; k<nz + 1; k++)
		for (int j = 1; j<ny + 2; j++)
		for (int i = 1; i<nx + 1; i++)
		{
			kfine = 2 * k - 1;
			jfine = 2 * j - 1;
			ifine = 2 * i - 1;
						
			val1 = sqrt(pow(s3xn[ign + 1][kfine][jfine][ifine], 2) + pow(s3yn[ign + 1][kfine][jfine][ifine], 2) + pow(s3zn[ign + 1][kfine][jfine][ifine], 2));
			val2 = sqrt(pow(s3xn[ign + 1][kfine][jfine][ifine + 1], 2) + pow(s3yn[ign + 1][kfine][jfine][ifine + 1], 2) + pow(s3zn[ign + 1][kfine][jfine][ifine + 1], 2));
			val3 = sqrt(pow(s3xn[ign + 1][kfine + 1][jfine][ifine], 2) + pow(s3yn[ign + 1][kfine + 1][jfine][ifine], 2) + pow(s3zn[ign + 1][kfine + 1][jfine][ifine], 2));
			val4 = sqrt(pow(s3xn[ign + 1][kfine + 1][jfine][ifine + 1], 2) + pow(s3yn[ign + 1][kfine + 1][jfine][ifine + 1], 2) + pow(s3zn[ign + 1][kfine + 1][jfine][ifine + 1], 2));
			val = val1 + val2 + val3 + val4;
						
			xx3[ign][k][j][i] = (xx3[ign + 1][kfine][jfine][ifine] * val1 + xx3[ign + 1][kfine][jfine][ifine + 1] * val2
						+ xx3[ign + 1][kfine + 1][jfine][ifine] * val3 + xx3[ign + 1][kfine + 1][jfine][ifine + 1] * val4) / val;
						
			yy3[ign][k][j][i] = (yy3[ign + 1][kfine][jfine][ifine] * val1 + yy3[ign + 1][kfine][jfine][ifine + 1] * val2
						+ yy3[ign + 1][kfine + 1][jfine][ifine] * val3 + yy3[ign + 1][kfine + 1][jfine][ifine + 1] * val4) / val;
						
			zz3[ign][k][j][i] = (zz3[ign + 1][kfine][jfine][ifine] * val1 + zz3[ign + 1][kfine][jfine][ifine + 1] * val2
						+ zz3[ign + 1][kfine + 1][jfine][ifine] * val3 + zz3[ign + 1][kfine + 1][jfine][ifine + 1] * val4) / val;
						
			s3xn[ign][k][j][i] = s3xn[ign + 1][kfine][jfine][ifine] + s3xn[ign + 1][kfine][jfine][ifine + 1]
						+ s3xn[ign + 1][kfine + 1][jfine][ifine] + s3xn[ign + 1][kfine + 1][jfine][ifine + 1];
						
			s3yn[ign][k][j][i] = s3yn[ign + 1][kfine][jfine][ifine] + s3yn[ign + 1][kfine][jfine][ifine + 1]
						+ s3yn[ign + 1][kfine + 1][jfine][ifine] + s3yn[ign + 1][kfine + 1][jfine][ifine + 1];
						
			s3zn[ign][k][j][i] = s3zn[ign + 1][kfine][jfine][ifine] + s3zn[ign + 1][kfine][jfine][ifine + 1]
						+ s3zn[ign + 1][kfine + 1][jfine][ifine] + s3zn[ign + 1][kfine + 1][jfine][ifine + 1];
		}

			
		for (int k = 1; k<nz + 2; k++)
		for (int j = 1; j<ny + 1; j++)
		for (int i = 1; i<nx + 1; i++)
		{
			kfine = 2 * k - 1;
			jfine = 2 * j - 1;
			ifine = 2 * i - 1;
						
			val1 = sqrt(pow(s1xn[ign + 1][kfine][jfine][ifine], 2) + pow(s1yn[ign + 1][kfine][jfine][ifine], 2) + pow(s1zn[ign + 1][kfine][jfine][ifine], 2));
			val2 = sqrt(pow(s1xn[ign + 1][kfine][jfine][ifine + 1], 2) + pow(s1yn[ign + 1][kfine][jfine][ifine + 1], 2) + pow(s1zn[ign + 1][kfine][jfine][ifine + 1], 2));
			val3 = sqrt(pow(s1xn[ign + 1][kfine][jfine + 1][ifine], 2) + pow(s1yn[ign + 1][kfine][jfine + 1][ifine], 2) + pow(s1zn[ign + 1][kfine][jfine + 1][ifine], 2));
			val4 = sqrt(pow(s1xn[ign + 1][kfine][jfine + 1][ifine + 1], 2) + pow(s1yn[ign + 1][kfine][jfine + 1][ifine + 1], 2) + pow(s1zn[ign + 1][kfine][jfine + 1][ifine + 1], 2));
			val = val1 + val2 + val3 + val4;
						
			xx1[ign][k][j][i] = (xx1[ign + 1][kfine][jfine][ifine] * val1 + xx1[ign + 1][kfine][jfine][ifine + 1] * val2
						+ xx1[ign + 1][kfine][jfine + 1][ifine] * val3 + xx1[ign + 1][kfine][jfine + 1][ifine + 1] * val4) / val;
			yy1[ign][k][j][i] = (yy1[ign + 1][kfine][jfine][ifine] * val1 + yy1[ign + 1][kfine][jfine][ifine + 1] * val2
						+ yy1[ign + 1][kfine][jfine + 1][ifine] * val3 + yy1[ign + 1][kfine][jfine + 1][ifine + 1] * val4) / val;
			zz1[ign][k][j][i] = (zz1[ign + 1][kfine][jfine][ifine] * val1 + zz1[ign + 1][kfine][jfine][ifine + 1] * val2
						+ zz1[ign + 1][kfine][jfine + 1][ifine] * val3 + zz1[ign + 1][kfine][jfine + 1][ifine + 1] * val4) / val;
			s1xn[ign][k][j][i] = s1xn[ign + 1][kfine][jfine][ifine] + s1xn[ign + 1][kfine][jfine][ifine + 1]
						+ s1xn[ign + 1][kfine][jfine + 1][ifine] + s1xn[ign + 1][kfine][jfine + 1][ifine + 1];
						
			s1yn[ign][k][j][i] = s1yn[ign + 1][kfine][jfine][ifine] + s1yn[ign + 1][kfine][jfine][ifine + 1]
						+ s1yn[ign + 1][kfine][jfine + 1][ifine] + s1yn[ign + 1][kfine][jfine + 1][ifine + 1];
						
			s1zn[ign][k][j][i] = s1zn[ign + 1][kfine][jfine][ifine] + s1zn[ign + 1][kfine][jfine][ifine + 1]
						+ s1zn[ign + 1][kfine][jfine + 1][ifine] + s1zn[ign + 1][kfine][jfine + 1][ifine + 1];

		}

		for (int k = 1; k<nz + 1; k++)
		for (int j = 1; j<ny + 1; j++)
		for (int i = 1; i<nx + 1; i++)
		{
			kfine = 2 * k - 1;
			jfine = 2 * j - 1;
			ifine = 2 * i - 1;
						
			val1 = vvn[ign + 1][kfine][jfine][ifine];
			val2 = vvn[ign + 1][kfine][jfine][ifine + 1];
			val3 = vvn[ign + 1][kfine][jfine + 1][ifine];
			val4 = vvn[ign + 1][kfine][jfine + 1][ifine + 1];
			val5 = vvn[ign + 1][kfine + 1][jfine][ifine];
			val6 = vvn[ign + 1][kfine + 1][jfine][ifine + 1];
			val7 = vvn[ign + 1][kfine + 1][jfine + 1][ifine];
			val8 = vvn[ign + 1][kfine + 1][jfine + 1][ifine + 1];
			val = val1 + val2 + val3 + val4 + val5 + val6 + val7 + val8;
						
			xx[ign][k][j][i] = (xx[ign + 1][kfine][jfine][ifine] * val1 + xx[ign + 1][kfine][jfine][ifine + 1] * val2
						+ xx[ign + 1][kfine][jfine + 1] [ifine]* val3 + xx[ign + 1][kfine][jfine + 1][ifine + 1] * val4
						+ xx[ign + 1][kfine + 1][jfine][ifine] * val5 + xx[ign + 1][kfine + 1][jfine][ifine + 1] * val6
						+ xx[ign + 1][kfine + 1][jfine + 1][ifine] * val7 + xx[ign + 1][kfine + 1][jfine + 1][ifine + 1] * val8) / val;
						
			yy[ign][k][j][i] = (yy[ign + 1][kfine][jfine][ifine] * val1 + yy[ign + 1][kfine][jfine][ifine + 1] * val2
						+ yy[ign + 1][kfine][jfine + 1][ifine] * val3 + yy[ign + 1][kfine][jfine + 1][ifine + 1] * val4
						+ yy[ign + 1][kfine + 1][jfine][ifine] * val5 + yy[ign + 1][kfine + 1][jfine][ifine + 1] * val6
						+ yy[ign + 1][kfine + 1][jfine + 1][ifine] * val7 + yy[ign + 1][kfine + 1][jfine + 1][ifine + 1] * val8) / val;
						
			zz[ign][k][j][i] = (zz[ign + 1][kfine][jfine][ifine] * val1 + zz[ign + 1][kfine][jfine][ifine + 1] * val2
						+ zz[ign + 1][kfine][jfine + 1][ifine] * val3 + zz[ign + 1][kfine][jfine + 1][ifine + 1] * val4
						+ zz[ign + 1][kfine + 1][jfine][ifine] * val5 + zz[ign + 1][kfine + 1][jfine][ifine + 1] * val6
						+ zz[ign + 1][kfine + 1][jfine + 1][ifine] * val7 + zz[ign + 1][kfine + 1][jfine + 1][ifine + 1] * val8) / val;
						
			vvn[ign][k][j][i] = val;	
		}
	}
}
	
bool CMultiGrid::ComputeGeoData(int nng,int lbb)	
{
	double temp, cor1, sir1;	
	int nx,ny,nz;
	


	if(nng !=1)
	mesh->SetlayerMesh(nng);

	ComputeSurface_CenterPoint(nng);

	V_interpolationPrepare(nng);

	
		
	temp = 2 / double(lbb)*pi;
	cor1 = cos(temp);
	sir1 = sin(temp);
		
	for (int ign = 1; ign<nng + 1; ign++)
	{
		nx = ((mesh->nnx))[ign];
		ny = ((mesh->nny))[ign];
		nz = ((mesh->nnz))[ign];
			
		for (int j = 1; j<ny + 1; j++)
		for (int i = 1; i<nx + 1; i++)
		{
			yy[ign][nz + 1][j][i] = yy[ign][1][j][i] * cor1 - zz[ign][1][j][i] * sir1;
			zz[ign][nz + 1][j][i] = zz[ign][1][j][i] * cor1 + yy[ign][1][j][i] * sir1;
			yy[ign][0][j][i] = yy[ign][nz][j][i] * cor1 + zz[ign][nz][j][i] * sir1;
			zz[ign][0][j][i] = zz[ign][nz][j][i] * cor1 - yy[ign][nz][j][i] * sir1;
		}		
	}
	return nx;	
}
	
bool  CMultiGrid::SetlayerGeoBoundary(int nnng)
{
	
	int nx = ((mesh->nnx))[nnng];
	int ny = ((mesh->nny))[nnng];
	int nz = ((mesh->nnz))[nnng];

	for (int k = 1; k<nz + 2; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		xx01[k][j][i] = xx1[nnng][k][j][i];
		yy01[k][j][i] = yy1[nnng][k][j][i];
		zz01[k][j][i] = zz1[nnng][k][j][i];
		//yy0[k][j][i] = yy[nnng][k][j][i];
		//zz0[k][j][i] = zz[nnng][k][j][i];
		s1x[k][j][i] = s1xn[nnng][k][j][i];
		s1y[k][j][i] = s1yn[nnng][k][j][i];
		s1z[k][j][i] = s1zn[nnng][k][j][i];
	}
	for (int k = 0; k<nz + 2; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		yy0[k][j][i] = yy[nnng][k][j][i];
		zz0[k][j][i] = zz[nnng][k][j][i];
	}

	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 2; i++)
	{
		xx02[k][j][i] = xx2[nnng][k][j][i];
		yy02[k][j][i] = yy2[nnng][k][j][i];
		zz02[k][j][i] = zz2[nnng][k][j][i];
		s2x[k][j][i] = s2xn[nnng][k][j][i];
		s2y[k][j][i] = s2yn[nnng][k][j][i];
		s2z[k][j][i] = s2zn[nnng][k][j][i];
	}
						
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 2; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		xx03[k][j][i] = xx3[nnng][k][j][i];
		yy03[k][j][i] = yy3[nnng][k][j][i];
		zz03[k][j][i] = zz3[nnng][k][j][i];
		s3x[k][j][i] = s3xn[nnng][k][j][i];
		s3y[k][j][i] = s3yn[nnng][k][j][i];
		s3z[k][j][i] = s3zn[nnng][k][j][i];
	}
								
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
	{										
		xx0[k][j][i] = xx[nnng][k][j][i];
		vv[k][j][i] = vvn[nnng][k][j][i];									
	}

																				
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	{
		s2x[k][j][0] = s2x[k][j][1];
		s2y[k][j][0] = s2y[k][j][1];
		s2z[k][j][0] = s2z[k][j][1];
		s2x[k][j][nx + 2] = s2x[k][j][nx + 1];
		s2y[k][j][nx + 2] = s2y[k][j][nx + 1];
		s2z[k][j][nx + 2] = s2z[k][j][nx + 1];
	}
																					
	for (int k = 1; k<nz + 1; k++)
	for (int i = 1; i<nx + 1; i++)
	{
		s3x[k][0][i] = s3x[k][1][i];
		s3y[k][0][i] = s3y[k][1][i];
		s3z[k][0][i] = s3z[k][1][i];
		s3x[k][ny + 2][i] = s3x[k][ny + 1][i];
		s3y[k][ny + 2][i] = s3y[k][ny + 1][i];
		s3z[k][ny + 2][i] = s3z[k][ny + 1][i];
	}
																						
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		s1x[0][j][i] = s1x[1][j][i];
		s1y[0][j][i] = s1y[1][j][i];
		s1z[0][j][i] = s1z[1][j][i];
		s1x[nz + 2][j][i] = s1x[nz + 1][j][i];
		s1y[nz + 2][j][i] = s1y[nz + 1][j][i];
		s1z[nz + 2][j][i] = s1z[nz + 1][j][i];
	}
	return nx;																							
}



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



