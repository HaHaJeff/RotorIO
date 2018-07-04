/*---------------------------------------------------------------------------*\


\*---------------------------------------------------------------------------*/

#ifndef TimeDerivation_H
#define TimeDerivation_H
#include "RKBase.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //




/*---------------------------------------------------------------------------*\
                           Class TimeDerivation Declaration
\*---------------------------------------------------------------------------*/


class CTimeDerivation:
	virtual public CRKBase
{
private:

public:
	double** dm;
	double*** ts1, ***ts2, ***ts3, ***ts4, ***ts5, ***ts6;

public:

    // Constructors

        CTimeDerivation(CMultiGrid* pgrid,CDictionary* pdict,CGeoField* pfield,int myid):CRKBase(pgrid,pdict,pfield, myid)
	{	
		Malloc(dm,pfield->nt+1,pfield->nt+1);

		Malloc(ts1,pfield->mesh->nzm1,pfield->mesh->nym1,pfield->mesh->nxm1);
		Malloc(ts2,pfield->mesh->nzm1,pfield->mesh->nym1,pfield->mesh->nxm1);
		Malloc(ts3,pfield->mesh->nzm1,pfield->mesh->nym1,pfield->mesh->nxm1);
		Malloc(ts4,pfield->mesh->nzm1,pfield->mesh->nym1,pfield->mesh->nxm1);
		Malloc(ts5,pfield->mesh->nzm1,pfield->mesh->nym1,pfield->mesh->nxm1);
		Malloc(ts6,pfield->mesh->nzm1,pfield->mesh->nym1,pfield->mesh->nxm1);						
	}

	~CTimeDerivation()
	{
		Free(ts1);
		Free(ts2);
		Free(ts3);
		Free(ts4);
		Free(ts5);
		Free(ts6);
		Free(dm);
	}
	bool timedrivationCmpute();
	void SetTimespectrum();
	bool timedrivationCmputeSA();
private:


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
