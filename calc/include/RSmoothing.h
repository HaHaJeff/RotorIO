/*---------------------------------------------------------------------------*\


\*---------------------------------------------------------------------------*/

#ifndef RSmoothing_H
#define RSmoothing_H
#include "ArtificialVisc.h"
#include "Diffusion.h"
#include "Convective.h"
#include "TimeDerivation.h"
/*---------------------------------------------------------------------------*\
			Class RSmoothing Declaration
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class CRSmoothing:
	public virtual CDiffusion,public virtual CArtificialVisc,public virtual CConvective,public virtual CTimeDerivation
{
public:

	double*** py1,*** py2, ***py3, ***py4, ***py5, ***py6;
	double*** rr1, ***rr2, ***rr3, ***rr4, ***rr5;
	double*** qp1, ***qp2, ***qp3, ***qp4, ***qp5;	
	double*** q01, ***q02, ***q03, ***q04, ***q05, ***q06;
public:
	CRSmoothing(CMultiGrid* pgrid,CDictionary* pdict,CGeoField* pfield,int myid):
	CDiffusion(pgrid,pdict,pfield, myid), CArtificialVisc(pgrid,pdict,pfield, myid),
	CConvective(pgrid,pdict,pfield, myid), CTimeDerivation(pgrid,pdict,pfield, myid),CGradient(pgrid,pdict,pfield, myid),CLocalTimestp(pgrid,pdict,pfield, myid),CRKBase(pgrid,pdict,pfield, myid)
	{
		
		Malloc(py1, (pfield->mesh)->nzm1 , (pfield->mesh)->nym1 , (pfield->mesh)->nxm1);
		Malloc(py2, (pfield->mesh)->nzm1 , (pfield->mesh)->nym1 , (pfield->mesh)->nxm1);
		Malloc(py3, (pfield->mesh)->nzm1 , (pfield->mesh)->nym1 , (pfield->mesh)->nxm1);
		Malloc(py4, (pfield->mesh)->nzm1 , (pfield->mesh)->nym1 , (pfield->mesh)->nxm1);
		Malloc(py5, (pfield->mesh)->nzm1 , (pfield->mesh)->nym1 , (pfield->mesh)->nxm1);
		Malloc(py6, (pfield->mesh)->nzm1 , (pfield->mesh)->nym1 , (pfield->mesh)->nxm1);
		
		
		Malloc(rr1, (pfield->mesh)->nzm1 , (pfield->mesh)->nym1 , (pfield->mesh)->nxm1 );
		Malloc(rr2, (pfield->mesh)->nzm1 , (pfield->mesh)->nym1 , (pfield->mesh)->nxm1 );
		Malloc(rr3, (pfield->mesh)->nzm1 , (pfield->mesh)->nym1 , (pfield->mesh)->nxm1 );
		Malloc(rr4, (pfield->mesh)->nzm1 , (pfield->mesh)->nym1 , (pfield->mesh)->nxm1 );
		Malloc(rr5, (pfield->mesh)->nzm1 , (pfield->mesh)->nym1 , (pfield->mesh)->nxm1 );
		
		
		Malloc(qp1,(pfield->mesh)->nzm1 , (pfield->mesh)->nym1 , (pfield->mesh)->nxm1 );
		Malloc(qp2,(pfield->mesh)->nzm1 , (pfield->mesh)->nym1 , (pfield->mesh)->nxm1 );
		Malloc(qp3,(pfield->mesh)->nzm1 , (pfield->mesh)->nym1 , (pfield->mesh)->nxm1 );
		Malloc(qp4,(pfield->mesh)->nzm1 , (pfield->mesh)->nym1 , (pfield->mesh)->nxm1 );
		Malloc(qp5,(pfield->mesh)->nzm1 , (pfield->mesh)->nym1 , (pfield->mesh)->nxm1 );

		Malloc(q01, (pfield->mesh)->nzm1 , (pfield->mesh)->nym1 , (pfield->mesh)->nxm1 );
		Malloc(q02, (pfield->mesh)->nzm1 , (pfield->mesh)->nym1 , (pfield->mesh)->nxm1 );
		Malloc(q03, (pfield->mesh)->nzm1 , (pfield->mesh)->nym1 , (pfield->mesh)->nxm1 );
		Malloc(q04, (pfield->mesh)->nzm1 , (pfield->mesh)->nym1 , (pfield->mesh)->nxm1 );
		Malloc(q05, (pfield->mesh)->nzm1 , (pfield->mesh)->nym1 , (pfield->mesh)->nxm1 );
		Malloc(q06, (pfield->mesh)->nzm1 , (pfield->mesh)->nym1 , (pfield->mesh)->nxm1 );
	}
	~CRSmoothing()
	{
		Free(py1);
		Free(py2);
		Free(py3);
		Free(py4);
		Free(py5);
		Free(py6);
		Free(rr1);
		Free(rr2);
		Free(rr3);
		Free(rr4);
		Free(rr5);
		Free(qp1);
		Free(qp2);
		Free(qp3);
		Free(qp4);
		Free(qp5);
		Free(q01);
		Free(q02);
		Free(q03);
		Free(q04);
		Free(q05);
		Free(q06);
	}

	void UpdateFieldValue(double ta,double timl,bool smoothing);
	void UpdateFieldValue_SA(double ta,double timl,bool smoothing);	
private:	
	void InitResidual();
	void SmoothResidual(double ta);
	void SmoothResidual_SA(double ta);
	void tdma(double*& a, double*& b, double*& c, double*& d, double*& x, int n);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
