/*---------------------------------------------------------------------------*\


\*---------------------------------------------------------------------------*/

#ifndef MultiGrid_H
#define MultiGrid_H
#include "Mesh.h"
#include "Memory.cpp"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //




/*---------------------------------------------------------------------------*\
                           Class MultiGrid Declaration
\*---------------------------------------------------------------------------*/


class CMultiGrid:
public CMemory<double>
{

private:
	CMesh* mesh;
	double**** s1xn, ****s1yn, ****s1zn, ****s2xn, ****s2yn, ****s2zn, ****s3xn, ****s3yn, ****s3zn;
        double**** xx1, ****yy1, ****zz1, ****xx2, ****xx3, ****yy3, ****zz3, ****xx, ****yy, ****zz;
	int ng;

public:


	double ****yy2, ****zz2;
	double**** vvn;
	double*** s1x, ***s1y, ***s1z, ***s2x, ***s2y, ***s2z, ***s3x, ***s3y, ***s3z, ***vv;
	double*** xx01, ***yy01, ***zz01, ***xx02, ***yy02, ***zz02, ***xx03, ***yy03, ***zz03, ***xx0, ***yy0,***zz0;

public:

    // Constructors

        CMultiGrid(CMesh* mesh,int layerNum)
	:mesh(mesh),ng(layerNum)
        {

		Malloc(s1x,mesh->nzm1+2,mesh->nym1,mesh->nxm1);
		Malloc(s1y,mesh->nzm1+2,mesh->nym1,mesh->nxm1);
		Malloc(s1z,mesh->nzm1+2,mesh->nym1,mesh->nxm1);

		Malloc(s3x,mesh->nzm1,mesh->nym1+2,mesh->nxm1);
		Malloc(s3y,mesh->nzm1,mesh->nym1+2,mesh->nxm1);
		Malloc(s3z,mesh->nzm1,mesh->nym1+2,mesh->nxm1);
		
		Malloc(s2x,mesh->nzm1,mesh->nym1,mesh->nxm1+2);
		Malloc(s2y,mesh->nzm1,mesh->nym1,mesh->nxm1+2);
		Malloc(s2z,mesh->nzm1,mesh->nym1,mesh->nxm1+2);
	
		Malloc(xx01,mesh->nzm1+1,mesh->nym1,mesh->nxm1);
		Malloc(yy01,mesh->nzm1+1,mesh->nym1,mesh->nxm1);
		Malloc(zz01,mesh->nzm1+1,mesh->nym1,mesh->nxm1);			

		Malloc(xx03,mesh->nzm1,mesh->nym1+1,mesh->nxm1);
		Malloc(yy03,mesh->nzm1,mesh->nym1+1,mesh->nxm1);
		Malloc(zz03,mesh->nzm1,mesh->nym1+1,mesh->nxm1);
		
		Malloc(xx02,mesh->nzm1,mesh->nym1,mesh->nxm1+1);
		Malloc(yy02,mesh->nzm1,mesh->nym1,mesh->nxm1+1);
		Malloc(zz02,mesh->nzm1,mesh->nym1,mesh->nxm1+1);

		Malloc(xx0,mesh->nzm1+1,mesh->nym1+1,mesh->nxm1+1);
		Malloc(yy0,mesh->nzm1+1,mesh->nym1+1,mesh->nxm1+1);
		Malloc(zz0,mesh->nzm1+1,mesh->nym1+1,mesh->nxm1+1);				

		Malloc(vv,mesh->nzm1,mesh->nym1,mesh->nxm1);

		Malloc(xx1,ng+1,mesh->nzm1+1,mesh->nym1,mesh->nxm1);
		Malloc(yy1,ng+1,mesh->nzm1+1,mesh->nym1,mesh->nxm1);
		Malloc(zz1,ng+1,mesh->nzm1+1,mesh->nym1,mesh->nxm1);
		Malloc(s1xn,ng+1,mesh->nzm1+1,mesh->nym1,mesh->nxm1);
		Malloc(s1yn,ng+1,mesh->nzm1+1,mesh->nym1,mesh->nxm1);
		Malloc(s1zn,ng+1,mesh->nzm1+1,mesh->nym1,mesh->nxm1);
		Malloc(yy,ng+1,mesh->nzm1+1,mesh->nym1,mesh->nxm1);
		Malloc(zz,ng+1,mesh->nzm1+1,mesh->nym1,mesh->nxm1);
	
		Malloc(xx3,ng+1,mesh->nzm1,mesh->nym1+1,mesh->nxm1);	
		Malloc(yy3,ng+1,mesh->nzm1,mesh->nym1+1,mesh->nxm1);
		Malloc(zz3,ng+1,mesh->nzm1,mesh->nym1+1,mesh->nxm1);
		Malloc(s3xn,ng+1,mesh->nzm1,mesh->nym1+1,mesh->nxm1);
		Malloc(s3yn,ng+1,mesh->nzm1,mesh->nym1+1,mesh->nxm1);
		Malloc(s3zn,ng+1,mesh->nzm1,mesh->nym1+1,mesh->nxm1);			
			
		Malloc(xx2,ng+1,mesh->nzm1,mesh->nym1,mesh->nxm1+1);
		Malloc(yy2,ng+1,mesh->nzm1,mesh->nym1,mesh->nxm1+1);
		Malloc(zz2,ng+1,mesh->nzm1,mesh->nym1,mesh->nxm1+1);
		Malloc(s2xn,ng+1,mesh->nzm1,mesh->nym1,mesh->nxm1+1);
		Malloc(s2yn,ng+1,mesh->nzm1,mesh->nym1,mesh->nxm1+1);
		Malloc(s2zn,ng+1,mesh->nzm1,mesh->nym1,mesh->nxm1+1);	
		
		Malloc(xx,ng+1,mesh->nzm1,mesh->nym1,mesh->nxm1);
		Malloc(vvn,ng+1,mesh->nzm1,mesh->nym1,mesh->nxm1);										
	
	}

	~CMultiGrid()
	{
		Free(s1xn);
		Free(s1yn);
		Free(s1zn);
		Free(s2xn);
		Free(s2yn);
		Free(s2zn);
		Free(s3xn);
		Free(s3yn);
		Free(s3zn);
		Free(xx1);
		Free(yy1);
		Free(zz1);
		Free(xx2);
		Free(xx3);
		Free(yy3);
		Free(zz3);
		Free(xx);
		Free(yy);
		Free(zz);
		Free(yy2);
		Free(zz2);
		Free(vvn);
		Free(s1x);
		Free(s1y);
		Free(s1z);
		Free(s2x);
		Free(s2y);
		Free(s2z);
		Free(s3x);
		Free(s3y);
		Free(s3z);
		Free(vv);
		Free(xx01);
		Free(yy01);
		Free(zz01);
		Free(xx02);
		Free(yy02);
		Free(zz02);
		Free(xx03);
		Free(yy03);
		Free(zz03);
		Free(xx0);
		Free(yy0);
		Free(zz0);
	}
	bool ComputeGeoData(int nng,int llb);
	bool SetlayerGeoBoundary(int nnng);

private:
	void ComputeSurface_CenterPoint(int nng);
	void V_interpolationPrepare(int nng);

public:
	friend class CGeoField;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
