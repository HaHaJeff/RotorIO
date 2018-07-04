/*---------------------------------------------------------------------------*\


\*---------------------------------------------------------------------------*/

#ifndef ArtificialVisc_H
#define ArtificialVisc_H
#include "LocalTimestp.h"

/*---------------------------------------------------------------------------*\
		  	Class ArtificialVisc Declaration
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
class CArtificialVisc:
	public virtual CLocalTimestp
{
private:

	double a2,a4;

public:


	double*** av1, ***av2, ***av3, ***av4, ***av5, ***av6;


public:
	CArtificialVisc(CMultiGrid* pgrid,CDictionary* pdict,CGeoField* pfield,int myid):
	CLocalTimestp(pgrid,pdict,pfield, myid),CRKBase(pgrid,pdict,pfield, myid),a2(pdict->a2),a4(pdict->a4)
	{
		
		Malloc(av1, (pfield->mesh)->nzm1 + 1, (pfield->mesh)->nym1 + 1, (pfield->mesh)->nxm1 + 1);
		Malloc(av2, (pfield->mesh)->nzm1 + 1, (pfield->mesh)->nym1 + 1, (pfield->mesh)->nxm1 + 1);
		Malloc(av3, (pfield->mesh)->nzm1 + 1, (pfield->mesh)->nym1 + 1, (pfield->mesh)->nxm1 + 1);
		Malloc(av4, (pfield->mesh)->nzm1 + 1, (pfield->mesh)->nym1 + 1, (pfield->mesh)->nxm1 + 1);
		Malloc(av5, (pfield->mesh)->nzm1 + 1, (pfield->mesh)->nym1 + 1, (pfield->mesh)->nxm1 + 1);
		Malloc(av6, (pfield->mesh)->nzm1 + 1, (pfield->mesh)->nym1 + 1, (pfield->mesh)->nxm1 + 1);
		
	}
	~CArtificialVisc()
	{
		Free(av1);
		Free(av2);
		Free(av3);
		Free(av4);
		Free(av5);
		Free(av6);
	}
	
	
public:
	bool step_artificial();
	bool step_artificial_c();
	bool step_artificial_SA();
private:

	void ResetViscousFlux();
	void ResetViscousFlux_SA();

	void ComputeViscousFlux_i();
	void ComputeViscousFlux_j();
	void ComputeViscousFlux_k();
	void ComputeViscousFlux();

	void ComputeViscousFlux_ci();
	void ComputeViscousFlux_cj();
	void ComputeViscousFlux_ck();
	void ComputeViscousFlux_c();

	void ComputeViscousFlux_SAi();
	void ComputeViscousFlux_SAj();
	void ComputeViscousFlux_SAk();
	void ComputeViscousFlux_SA();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
