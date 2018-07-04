/*---------------------------------------------------------------------------*\


\*---------------------------------------------------------------------------*/

#ifndef RungeKutta_H
#define RungeKutta_H
#include "RSmoothing.h"
#include "BoundaryCondition.h"
#include "writefile.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

/*---------------------------------------------------------------------------*\
                           Class RKBase Declaration
\*---------------------------------------------------------------------------*/


class CRungeKutta:
	public CBoundaryCondition,public CRSmoothing,public CWriteFile
{
private:
	double ta, timl;	
public:
	double**** q31, ****q32, ****q33, ****q34, ****q35, ****q36;
	double* rms;
public:


public:
	CRungeKutta(CMultiGrid* pgrid, CDictionary* pdict, CGeoField* pfield,int numprocs,int myid):CRKBase(pgrid,pdict,pfield,myid),
		CGradient(pgrid, pdict, pfield,myid), CDiffusion(pgrid, pdict, pfield, myid), CLocalTimestp(pgrid, pdict, pfield, myid),
		CArtificialVisc(pgrid, pdict, pfield, myid), CConvective(pgrid, pdict, pfield, myid), CTimeDerivation(pgrid, pdict, pfield, myid),
		CBoundaryCondition(pfield,pgrid,  pdict, numprocs, myid), CRSmoothing(pgrid, pdict, pfield, myid),CWriteFile(pgrid, pdict, pfield, myid)
	{
	
		Malloc(rms,pdict->nt + 1);

		Malloc(q31,pdict->nt + 1, (pfield->mesh)->nzm1 + 1, (pfield->mesh)->nym1 + 1, (pfield->mesh)->nxm1 + 1);
		Malloc(q32,pdict->nt + 1, (pfield->mesh)->nzm1 + 1, (pfield->mesh)->nym1 + 1, (pfield->mesh)->nxm1 + 1);
		Malloc(q33,pdict->nt + 1, (pfield->mesh)->nzm1 + 1, (pfield->mesh)->nym1 + 1, (pfield->mesh)->nxm1 + 1);
		Malloc(q34,pdict->nt + 1, (pfield->mesh)->nzm1 + 1, (pfield->mesh)->nym1 + 1, (pfield->mesh)->nxm1 + 1);
		Malloc(q35,pdict->nt + 1, (pfield->mesh)->nzm1 + 1, (pfield->mesh)->nym1 + 1, (pfield->mesh)->nxm1 + 1);
		Malloc(q36,pdict->nt + 1, (pfield->mesh)->nzm1 + 1, (pfield->mesh)->nym1 + 1, (pfield->mesh)->nxm1 + 1);

	}
	~CRungeKutta()
	{
		Free(q31);	
		Free(q32);	
		Free(q33);	
		Free(q34);	
		Free(q35);	
		Free(q36);	
		Free(rms);	
	}
public:
	void RungeKutta_I( int nng);
	void RungeKutta_II(  int nng);
	void RungeKutta_SA();
	void FMGCycle( double time_begin,double time_end);
	void residual(double &rmsm);
	void RestrictionOp(int ign);
	void P2hSetting();
	void UpwardCorrection(int ign,int nng);
	void update(int nng);
private:
	void ResetP2h();
	void InitPressure();
	void InitPreField();
	void InitPreField_SA();

	void RK_I();
	void RK_C_I();
	void RK_II();
	void RK_III();

	inline double dp(double a, double b, double c)
	{
		return (3.0*a - 2.0*b - c) / 64.0;
	}

private:


};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
