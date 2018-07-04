/*---------------------------------------------------------------------------*\


\*---------------------------------------------------------------------------*/


#include "Mesh.h"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool CMesh::getOrgData()
{
		ifstream fin1(file1.c_str());
		fin1 >> nxm1 >> nym1 >> nzm1;


		initvector(nzm1+1,nym1+1,nxm1+1);

		forAll(0,nzm1 + 1,i)
		forAll(0,nym1 + 1,j)
		forAll(0,nxm1 + 1,k)
		{
			if(k==0||j==0||i==0)						
			{						
				(xf)[i][j][k]=(yf)[i][j][k]=(zf)[i][j][k]=0;
			}	
			else						
			{						
				fin1 >> (xf)[i][j][k] >> (yf)[i][j][k] >> (zf)[i][j][k];						
			}
		}
		         


		fin1.close();
	return nxm1&&nym1&&nzm1;
}

void CMesh::GenerateMesh(int myid,int lbb)
{		
	
		if(myid !=0){
		double temp, cor1, sir1, t1, t2;

		temp = 2.0 / double(lbb)*pi*myid;
		cor1 = cos(temp);				
		sir1 = sin(temp);


		forAll(1,nzm1+1,i)
		forAll(1,nym1+1,j)
		forAll(1,nxm1+1,k)
		{
			t1=(yf)[i][j][k];
			t2=(zf)[i][j][k];
			(yf)[i][j][k] = t1*cor1 - t2*sir1;
			(zf)[i][j][k] = t1*sir1 + t2*cor1;
		}
		}

		
}

bool CMesh::initvector(int size1,int size2,int size3)
{
	Malloc(xf, size1, size2, size3);
	Malloc(yf, size1, size2, size3);
	Malloc(zf, size1, size2, size3);
	Malloc(x, size1, size2, size3);
	Malloc(y, size1, size2, size3);
	Malloc(z, size1, size2, size3);
	Malloc(dmini, size1, size2, size3);	
 	return true;
}

bool CMesh::SetMultiGrid(int ibm,int jbm,int itm,int jtm)
{
	double divnum;
	itm = itm - 1;
	jtm = jtm - 1;
	
	forAll(1,ng+1,j)
	{
		divnum=pow(2.0, (ng - j));
		(nnx)[j] = (nxm1-1) / divnum;
		(nny)[j] = (nym1-1) / divnum;
		(nnz)[j] = (nzm1-1) / divnum;
		(nib)[j] = (ibm - 1) / divnum + 1;
		(nit)[j] = itm / divnum;
		(njb)[j] = (jbm - 1) / divnum + 1;
		(njt)[j] = jtm / divnum;
	}
	return divnum;
}

bool CMesh::SetlayerMesh(int LayerRange)
{
	int step, ifine, jfine, kfine;
	int nx = (nnx)[LayerRange];
	int ny = (nny)[LayerRange];
	int nz = (nnz)[LayerRange];
		
	step = pow(2.0, (ng - LayerRange));
	for (int k = 1; k<nz + 2; k++)	
	for (int j = 1; j<ny + 2; j++)
	for (int i = 1; i<nx + 2; i++)
	{
		kfine = step*(k - 1) + 1;
		jfine = step*(j - 1) + 1;
		ifine = step*(i - 1) + 1;
		(x)[k][j][i] = (xf)[kfine][jfine][ifine];
		(y)[k][j][i] = (yf)[kfine][jfine][ifine];
		(z)[k][j][i] = (zf)[kfine][jfine][ifine];
		
	}
	return (x)[nz][ny][nx];
}


double*** CMesh::min_distance()		
{
	int nx = (nnx)[ng];
    	int ny = (nny)[ng];
	int nz = (nnz)[ng];
	int ib = (nib)[ng];
	int it = (nit)[ng];
	int jb = (njb)[ng];
	int jt = (njt)[ng];
	double** xwu,**ywu,**zwu,**xwd,**ywd,**zwd;
	Malloc(xwu,nz+1,nx+1);
	Malloc(ywu,nz+1,nx+1);
	Malloc(zwu,nz+1,nx+1);
	Malloc(xwd,nz+1,nx+1);
	Malloc(ywd,nz+1,nx+1);
	Malloc(zwd,nz+1,nx+1);
	double** xwf,**ywf,**zwf,**xwb,**ywb,**zwb;
	Malloc(xwf,jt + 1,it + 1);
	Malloc(ywf,jt + 1,it + 1);
	Malloc(zwf,jt + 1,it + 1);
	Malloc(xwb,jt + 1,it + 1);
	Malloc(ywb,jt + 1,it + 1);
	Malloc(zwb,jt + 1,it + 1);
	double*** xx00,***yy00,***zz00;
	Malloc(xx00,nz+1,ny+1,nx+1);
	Malloc(yy00,nz+1,ny+1,nx+1);
	Malloc(zz00,nz+1,ny+1,nx+1);


	int ii, jj, kk, iil, iir, kkl, kkr, iib, iit, jjb, jjt;		
	double dy1, dy2, dy, dz1, dz2, dz;
	double d;


	
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		xx00[k][j][i] = 0.125*((xf)[k][j][i] + (xf)[k][j][i + 1] + (xf)[k + 1][j][i] + (xf)[k + 1][j][i + 1]
					+ (xf)[k][j + 1][i] + (xf)[k][j + 1][i + 1] + (xf)[k + 1][j + 1][i] + (xf)[k + 1][j + 1][i + 1]);
					
		yy00[k][j][i] = 0.125*((yf)[k][j][i] + (yf)[k][j][i + 1] + (yf)[k + 1][j][i] + (yf)[k + 1][j][i + 1]
						+ (yf)[k][j + 1][i] + (yf)[k][j + 1][i + 1] + (yf)[k + 1][j + 1][i] + (yf)[k + 1][j + 1][i + 1]);
					
		zz00[k][j][i] = 0.125*((zf)[k][j][i] + (zf)[k][j][i + 1] + (zf)[k + 1][j][i] + (zf)[k + 1][j][i + 1]
						+ (zf)[k][j + 1][i] + (zf)[k][j + 1][i + 1] + (zf)[k + 1][j + 1][i] + (zf)[k + 1][j + 1][i + 1]);
	}

	for (int k = 1; k<nz + 1; k++)
	for (int i = 1; i<nx + 1; i++)
	{
		xwd[k][i] = 0.25*((xf)[k][1][i] + (xf)[k][1][i + 1] + (xf)[k + 1][1][i] + (xf)[k + 1][1][i + 1]);
		ywd[k][i] = 0.25*((yf)[k][1][i] + (yf)[k][1][i + 1] + (yf)[k + 1][1][i] + (yf)[k + 1][1][i + 1]);
		zwd[k][i] = 0.25*((zf)[k][1][i] + (zf)[k][1][i + 1] + (zf)[k + 1][1][i] + (zf)[k + 1][1][i + 1]);

		xwu[k][i] = 0.25*((xf)[k][ny + 1][i] + (xf)[k][ny + 1][i + 1] + (xf)[k + 1][ny + 1][i] + (xf)[k + 1][ny + 1][i + 1]);
		ywu[k][i] = 0.25*((yf)[k][ny + 1][i] + (yf)[k][ny + 1][i + 1] + (yf)[k + 1][ny + 1][i] + (yf)[k + 1][ny + 1][i + 1]);
		zwu[k][i] = 0.25*((zf)[k][ny + 1][i] + (zf)[k][ny + 1][i + 1] + (zf)[k + 1][ny + 1][i] + (zf)[k + 1][ny + 1][i + 1]);
	}
							
	for (int j = jb; j<jt + 1; j++)
	for (int i = ib; i<it + 1; i++)
	{
		xwf[j][i] = 0.25*((xf)[1][j][i] + (xf)[1][j][i + 1] + (xf)[1][j + 1][i] + (xf)[1][j + 1][i + 1]);
		ywf[j][i] = 0.25*((yf)[1][j][i] + (yf)[1][j][i + 1] + (yf)[1][j + 1][i] + (yf)[1][j + 1][i + 1]);
		zwf[j][i] = 0.25*((zf)[1][j][i] + (zf)[1][j][i + 1] + (zf)[1][j + 1][i] + (zf)[1][j + 1][i + 1]);
		xwb[j][i] = 0.25*((xf)[nz + 1][j][i] + (xf)[nz + 1][j][i + 1] + (xf)[nz + 1][j + 1][i] + (xf)[nz + 1][j + 1][i + 1]);
		ywb[j][i] = 0.25*((yf)[nz + 1][j][i] + (yf)[nz + 1][j][i + 1] + (yf)[nz + 1][j + 1][i] + (yf)[nz + 1][j + 1][i + 1]);
		zwb[j][i] = 0.25*((zf)[nz + 1][j][i] + (zf)[nz + 1][j][i + 1] + (zf)[nz + 1][j + 1][i] + (zf)[nz + 1][j + 1][i + 1]);
	}
								
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		dy = pow(10.0, 20);
		iil = max(i - 10, 1);
		iir = min(iil + 20, nx);
		kkl = max(k - 10, 1);
		kkr = min(kkl + 20, nz);
		for (kk = kkl; kk<kkr + 1; kk++)
		for (ii = iil; ii<iir + 1; ii++)
		{
			dy1 = pow((xx00[k][j][i] - xwd[kk][ii]), 2) + pow((yy00[k][j][i] - ywd[kk][ii]), 2) + pow((zz00[k][j][i] - zwd[kk][ii]), 2);
			dy2 = pow((xx00[k][j][i] - xwu[kk][ii]), 2) + pow((yy00[k][j][i] - ywu[kk][ii]), 2) + pow((zz00[k][j][i] - zwu[kk][ii]), 2);
			d = min(dy1, dy2);
			if (d<dy)
				dy = d;
		}

		dz = pow(10.0, 20);
		if (j<jb)
		{
			jjb = jb;
			jjt = jb + 20;
		}
		else if ((j >= jb) && (j <= jt))
		{
			jjb = max(j - 10, jb);
			jjt = min(jjb + 20, jt);
		}
		else
		{
			jjb = jt - 20;
			jjt = jt;
		}

		if (i<ib)
		{
			iib = ib;
			iit = ib + 20;
													
		}
		else if ((i >= ib) && (i <= it))
		{
			iib = max(i - 10, ib);
			iit = min(iib + 20, it);
													
		}
		else
		{
			iib = it - 20;
			iit = it;
		}
	
		for (jj = jjb; jj<jjt + 1; jj++)
		for (ii = iib; ii<iit + 1; ii++)
		{
			dz1 = (pow((xx00[k][j][i] - xwf[jj][ii]), 2) + pow((yy00[k][j][i] - ywf[jj][ii]), 2) + pow((zz00[k][j][i] - zwf[jj][ii]), 2));
			dz2 = (pow((xx00[k][j][i] - xwb[jj][ii]), 2) + pow((yy00[k][j][i] - ywb[jj][ii]), 2) + pow((zz00[k][j][i] - zwb[jj][ii]), 2));
			d = min(dz1, dz2);
			if (d<dz)
				dz = d;
		}
													
		(dmini)[k][j][i] = sqrt(min(dy, dz));

	}

	Free(xwu);
	Free(ywu);
	Free(zwu);
	Free(xwd);
	Free(ywd);
	Free(zwd);
	Free(xwf);
	Free(ywf);
	Free(zwf);
	Free(xwb);
	Free(ywb);
	Free(zwb);
	Free(xx00);
	Free(yy00);
	Free(zz00);
	
	return dmini;
		
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



