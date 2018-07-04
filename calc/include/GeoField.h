/*---------------------------------------------------------------------------*\


\*---------------------------------------------------------------------------*/

#ifndef GeoField_H
#define GeoField_H
#include "Mesh.h"
#include "MultiGrid.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //




/*---------------------------------------------------------------------------*\
                           Class CGeoField Declaration
\*---------------------------------------------------------------------------*/

class CGeoField:
public CMemory<double>
{


public:
	double**** betaxn, ****betayn, ****betazn, ****hatn, ****petn;
	CMesh* mesh;
	int ng,nt;
	double**** q11, ****q12, ****q13, ****q14, ****q15, ****q16;
	double** betax, **betay, **betaz, **pet, **hat;
	double***  pvx, ***pvy, ***pvz, ***vth, ***vre, ***p, ***t;

public:

    // Constructors

        CGeoField(CMesh* mesh,int nt)
	:mesh(mesh),nt(nt),ng(mesh->ng)
        {
		Malloc(q11,nt+1,mesh->nzm1+1,mesh->nym1+1,mesh->nxm1+1);
		Malloc(q12,nt+1,mesh->nzm1+1,mesh->nym1+1,mesh->nxm1+1);
		Malloc(q13,nt+1,mesh->nzm1+1,mesh->nym1+1,mesh->nxm1+1);
		Malloc(q14,nt+1,mesh->nzm1+1,mesh->nym1+1,mesh->nxm1+1);
		Malloc(q15,nt+1,mesh->nzm1+1,mesh->nym1+1,mesh->nxm1+1);
		Malloc(q16,nt+1,mesh->nzm1+1,mesh->nym1+1,mesh->nxm1+1);	

		Malloc(betax,mesh->nzm1+1,mesh->nym1+1);
		Malloc(betay,mesh->nzm1+1,mesh->nym1+1);
		Malloc(betaz,mesh->nzm1+1,mesh->nym1+1);
		Malloc(hat,mesh->nzm1+1,mesh->nym1+1);
		Malloc(pet,mesh->nzm1+1,mesh->nym1+1);
		
		Malloc(betaxn,nt+1,ng+1,mesh->nzm1,mesh->nym1);
		Malloc(betayn,nt+1,ng+1,mesh->nzm1,mesh->nym1);
		Malloc(betazn,nt+1,ng+1,mesh->nzm1,mesh->nym1);
		Malloc(petn,nt+1,ng+1,mesh->nzm1,mesh->nym1);
		Malloc(hatn,nt+1,ng+1,mesh->nzm1,mesh->nym1);
		
		Malloc(pvx,mesh->nzm1+1,mesh->nym1+1,mesh->nxm1+1);
		Malloc(pvy,mesh->nzm1+1,mesh->nym1+1,mesh->nxm1+1);
		Malloc(pvz,mesh->nzm1+1,mesh->nym1+1,mesh->nxm1+1);
		Malloc(vth,mesh->nzm1+1,mesh->nym1+1,mesh->nxm1+1);
		Malloc(vre,mesh->nzm1+1,mesh->nym1+1,mesh->nxm1+1);
		Malloc(p,mesh->nzm1+1,mesh->nym1+1,mesh->nxm1+1);
		Malloc(t,mesh->nzm1+1,mesh->nym1+1,mesh->nxm1+1);		
	}
	~CGeoField()
	{
		Free(betaxn);
		Free(betayn);
		Free(betazn);
		Free(hatn);
		Free(petn);
		Free(q11);
		Free(q12);
		Free(q13);
		Free(q14);
		Free(q15);
		Free(q16);
		Free(betax);
		Free(betay);
		Free(betaz);
		Free(pet);
		Free(hat);
		Free(pvx);
		Free(pvy);
		Free(pvz);
		Free(vth);
		Free(vre);
		Free(p);
		Free(t);
	}
	
	bool SetPhysicalBoundary(int nng0,double vxx,double vrr, double vtt,double vee,double ht,double pt);
	bool FieldInitialization(int nng0,double ma,double cvl0);
	bool GetlayerPhysicalBoundary(int n,int nnng);
public:
	void SetPhysicalBoundary(int nng0, double vxx, double vrr, double vtt, double vee, double ht, double pt, CMultiGrid *p);
private:


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
