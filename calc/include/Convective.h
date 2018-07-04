/*---------------------------------------------------------------------------*\


\*---------------------------------------------------------------------------*/

#ifndef Convective_H
#define Convective_H
#include "RKBase.h"
/*---------------------------------------------------------------------------*\
			Class convective Declaration
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

class CConvective:
	virtual public CRKBase
{
private:


public:
	double*** qc1, ***qc2, ***qc3, ***qc4, ***qc5, ***qc6;
public:
	CConvective(CMultiGrid* pgrid,CDictionary* pdict,CGeoField* pfield,int myid):CRKBase(pgrid,pdict,pfield, myid)
	{
		
		Malloc(qc1, (pfield->mesh)->nzm1 + 1, (pfield->mesh)->nym1 + 1, (pfield->mesh)->nxm1 + 1);
		Malloc(qc2, (pfield->mesh)->nzm1 + 1, (pfield->mesh)->nym1 + 1, (pfield->mesh)->nxm1 + 1);
		Malloc(qc3, (pfield->mesh)->nzm1 + 1, (pfield->mesh)->nym1 + 1, (pfield->mesh)->nxm1 + 1);
		Malloc(qc4, (pfield->mesh)->nzm1 + 1, (pfield->mesh)->nym1 + 1, (pfield->mesh)->nxm1 + 1);
		Malloc(qc5, (pfield->mesh)->nzm1 + 1, (pfield->mesh)->nym1 + 1, (pfield->mesh)->nxm1 + 1);
		Malloc(qc6, (pfield->mesh)->nzm1 + 1, (pfield->mesh)->nym1 + 1, (pfield->mesh)->nxm1 + 1);
	}
	~CConvective()
	{
		Free(qc1);
		Free(qc2);
		Free(qc3);
		Free(qc4);
		Free(qc5);
		Free(qc6);
	}

	bool ComputeConvectiveData_SA();
	bool ComputeConvectiveData();
private:
	void ResetConvectiveData();
	void ComputeConvectiveData_i();
	void ComputeConvectiveData_j();
	void ComputeConvectiveData_k();

	void ResetConvectiveData_SA();
	void ComputeConvectiveData_SAi();
	void ComputeConvectiveData_SAj();
	void ComputeConvectiveData_SAk();

public:

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

// ************************************************************************* //
