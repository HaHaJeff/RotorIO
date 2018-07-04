/*---------------------------------------------------------------------------*\


\*---------------------------------------------------------------------------*/

#ifndef RKBase_H
#define RKBase_H
#include "Memory.h"
#include "MultiGrid.h"
#include "GeoField.h"
#include "Dictionary.h"
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //




/*---------------------------------------------------------------------------*\
                           Class RKBase Declaration
\*---------------------------------------------------------------------------*/


class CRKBase:
	public virtual CMemory<double>
{
private:
	
public:


	int nx,ny,nz,n;
	int myid;
	int ib, it, jt, jb;
	int V_stage,iter,iteration;
	CMultiGrid* pgrid;
	CDictionary* pdict;
	CGeoField* pfield;


public:

    // Constructors

        CRKBase(CMultiGrid* pgrid,CDictionary* pdict,CGeoField* pfield,int myid)
	:pgrid(pgrid),pdict(pdict),pfield(pfield),nx(0),ny(0),nz(0),n(0),
	V_stage(0),iter(0),myid(myid),iteration(0)
        {
	};

	~CRKBase()
	{};

	bool SetMultigridRange(int nng);
	int V_Cycle();
	int Iteration();
	int Vitual_timeLoop();

private:


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
