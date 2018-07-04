/*---------------------------------------------------------------------------*\


\*---------------------------------------------------------------------------*/

#ifndef Diffusion_H
#define Diffusion_H
#include "Gradient.h"
/*---------------------------------------------------------------------------*\
			Class Diffusion Declaration
\*---------------------------------------------------------------------------*/

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
class CDiffusion:
	public virtual CGradient
{
private:

public:

	double*** qv2, ***qv3, ***qv4, ***qv5, ***qv6;


public:
	CDiffusion(CMultiGrid* pgrid,CDictionary* pdict,CGeoField* pfield,int myid):
	CGradient(pgrid,pdict,pfield, myid),CRKBase(pgrid,pdict,pfield, myid)
	{
		Malloc(qv2,(pfield->mesh)->nzm1 + 1, (pfield->mesh)->nym1 + 1, (pfield->mesh)->nxm1 + 1);
		Malloc(qv3,(pfield->mesh)->nzm1 + 1, (pfield->mesh)->nym1 + 1, (pfield->mesh)->nxm1 + 1);
		Malloc(qv4,(pfield->mesh)->nzm1 + 1, (pfield->mesh)->nym1 + 1, (pfield->mesh)->nxm1 + 1);
		Malloc(qv5,(pfield->mesh)->nzm1 + 1, (pfield->mesh)->nym1 + 1, (pfield->mesh)->nxm1 + 1);
		Malloc(qv6,(pfield->mesh)->nzm1 + 1, (pfield->mesh)->nym1 + 1, (pfield->mesh)->nxm1 + 1);

	}
	~CDiffusion()
	{
		Free(qv2);
		Free(qv3);
		Free(qv4);
		Free(qv5);
		Free(qv6);
	}

public:
	void ComputeViscousData();
	void ComputeViscousData_SA();
private:

	void ResetViscousData();
	void viscosity(double temp, double q6, double &cv, double &kc);
	void ComputeViscousData_i();
	void ComputeViscousData_j();
	void ComputeViscousData_k();
	void SAsource();
};






// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
