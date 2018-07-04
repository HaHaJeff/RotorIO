/*---------------------------------------------------------------------------*\


\*---------------------------------------------------------------------------*/

#ifndef Gradient_H
#define Gradient_H
#include "RKBase.h"
/*---------------------------------------------------------------------------*\
			Class Gradient Declaration
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
class CGradient:
	 public virtual CRKBase
{
private:

public:

	double**** gradfi, ****gradfj, ****gradfk;
	double**** gradc, ****gradcs;

	
public:
	CGradient(CMultiGrid* pgrid,CDictionary* pdict,CGeoField* pfield,int myid):CRKBase(pgrid,pdict,pfield,myid)
	{
		
		Malloc(gradfi,16 ,(pfield->mesh)->nzm1 + 1, (pfield->mesh)->nym1 + 1, (pfield->mesh)->nxm1 + 1);
		Malloc(gradfj,16 ,(pfield->mesh)->nzm1 + 1, (pfield->mesh)->nym1 + 1, (pfield->mesh)->nxm1 + 1);
		Malloc(gradfk,16 ,(pfield->mesh)->nzm1 + 1, (pfield->mesh)->nym1 + 1, (pfield->mesh)->nxm1 + 1);
		
		Malloc(gradc,13, (pfield->mesh)->nzm1 + 1, (pfield->mesh)->nym1 + 1, (pfield->mesh)->nxm1 + 1);
		Malloc(gradcs,10, (pfield->mesh)->nzm1 + 1, (pfield->mesh)->nym1 + 1, (pfield->mesh)->nxm1 + 1);
	}
	~CGradient()
	{
		Free(gradfi);
		Free(gradfj);
		Free(gradfk);
		Free(gradc);
		Free(gradcs);
	}
protected:
	
	void gradsfaceI(double***& si, double***& sj,double***& sk,double***& q,double****& Idqd, int m);
	
	void gradsfaceJ(double***& si, double***& sj,double***& sk,double***& q, double****& Jdqd, int m); 
	
	void gradsfaceK(double***& si, double***& sj,double***& sk,double***& q, double****& Kdqd, int m);
	void dsdt(double****&q, double****& s, int m,double rpm);
	void setBoundary(int nng);
public:
	void ResetGradient();
	void gradsface();
	void gradscentre(int direction, double***& q, double****& dqd, int m);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
