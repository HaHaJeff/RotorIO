/*---------------------------------------------------------------------------*\


\*---------------------------------------------------------------------------*/


#include "RKBase.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool CRKBase::SetMultigridRange(int nng)
{
	nx = (((pfield->mesh)->nnx))[nng];
	ny = (((pfield->mesh)->nny))[nng];
	nz = (((pfield->mesh)->nnz))[nng];
	ib = ((pfield->mesh->nib))[nng];
	it = ((pfield->mesh->nit))[nng];
	jb = ((pfield->mesh->njb))[nng];
	jt = ((pfield->mesh->njt))[nng];
	//this->n = n;
	return true;
}
int CRKBase::V_Cycle()
{		
	V_stage++;
	if(V_stage==pdict->ng+1)
	return 0;
	else
        return V_stage;
}

int CRKBase::Iteration()
{
	iter++;
	if(iter==((pdict->cg))[V_stage] + 1)
	{
		iter=0;
		return 0;
	}
	else
	{
		iteration++;
		return iter;
	}
}

int CRKBase::Vitual_timeLoop()
{
	n++;
	if(n == pdict->nt+1)
	{
		n=0;
		return 0;
	}
	else
        return n;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



