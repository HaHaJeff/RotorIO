/*---------------------------------------------------------------------------*\


\*---------------------------------------------------------------------------*/

#ifndef BoundaryCondition_H
#define BoundaryCondition_H
#include "RKBase.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                           Class BoundaryCondition Declaration
\*---------------------------------------------------------------------------*/


class CBoundaryCondition:
	public virtual CRKBase
{
private:
	
	double *hr,*hv,*hd,*hp;
	double **turi,**peb;
	int nzm1,nym1;
	int myidl, myidr;
	int numprocs;
public:

	MPI_Status status;

public:

    // Constructors

        CBoundaryCondition(CGeoField* pGeofield,CMultiGrid* pMultigrid,CDictionary* pDictionary,int numprocs,int myid)
	:CRKBase(pMultigrid,pDictionary, pGeofield, myid),nzm1(pGeofield->mesh->nzm1),nym1(pGeofield->mesh->nym1),numprocs(numprocs)
	{	
		myidl=myidr=0;
		Malloc(turi,nzm1,nym1);
		Malloc(peb,nzm1,nym1);		
		Malloc(hr,nym1);
		Malloc(hv,nym1);
		Malloc(hd,nym1);
		Malloc(hp,nym1);		
	}

	~CBoundaryCondition()
	{
		Free(hr);
		Free(hv);
		Free(hd);
		Free(hp);
		Free(turi);
		Free(peb);
	}
	
	bool ExchangeBoundaryCondition();
private:
	void setBC_in();
	void setBC_out();	
	void setBC_root_top();
	void setBC_circumference();
	void update();
	void communication();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
