/*---------------------------------------------------------------------------*\


\*---------------------------------------------------------------------------*/

#ifndef LocalTimestp_H
#define LocalTimestp_H
#include "RKBase.h"

/*---------------------------------------------------------------------------*\
		  	Class Timestep Declaration
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
class CLocalTimestp:
	virtual  public CRKBase
{
private:


public:

	double*** sri, ***srj, ***srk;
	double*** time;

	
public:
	CLocalTimestp(CMultiGrid* pgrid,CDictionary* pdict,CGeoField* pfield,int myid):CRKBase(pgrid,pdict,pfield, myid)
	{
		Malloc(sri,(pfield->mesh)->nzm1 + 1, (pfield->mesh)->nym1+1, (pfield->mesh)->nxm1+1);
		Malloc(srj,(pfield->mesh)->nzm1 + 1, (pfield->mesh)->nym1+1, (pfield->mesh)->nxm1+1);
		Malloc(srk,(pfield->mesh)->nzm1 + 1, (pfield->mesh)->nym1+1, (pfield->mesh)->nxm1+1);
		
		Malloc(time,(pfield->mesh)->nzm1, (pfield->mesh)->nym1, (pfield->mesh)->nxm1);
	}
	~CLocalTimestp()
	{
		Free(sri);
		Free(srj);
		Free(srk);
		Free(time);
	}
	bool SetLocalTimestep();
private:

	void InitTimeData();
	void ComputeTimeData();
	void viscosity(double t, double q6, double &cv, double &kc);
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
