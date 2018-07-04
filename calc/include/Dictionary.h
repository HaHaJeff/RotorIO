/*---------------------------------------------------------------------------*\


\*---------------------------------------------------------------------------*/

#ifndef Dictionary_H
#define Dictionary_H
#include "Memory.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //




/*---------------------------------------------------------------------------*\
                           Class Dictionary Declaration
\*---------------------------------------------------------------------------*/


class CDictionary:
public CMemory<int>
{
private:
	string file1;


public:
	double cvl0,t0,ts,rg,cp,prl,prt,cfl,a2, a4, ht ,pt, pb1, c2,rmsm0, ma,
		cv1,cv2,kap,cb1,cb2,cw1,cw2,cw3,cr1,cr2,cr3,sigmav,beta1,beta2,vxx,vrr,vtt,vee;
	double rpm;
	int nt,ng, nxm, nym, nzm, ibm,jbm,itm,jtm, lbb;
	


	int* cg;


public:

    // Constructors

        CDictionary(string s1)
	:file1(s1),nt(0),ng(0),beta1(0),beta2(0) ,cfl(0),a2(0), a4(0), ht(0) ,pt(0), pb1(0),
		c2(0),rmsm0(0), nxm(0), nym(0), nzm(0), ibm(0),jbm(0),itm(0),jtm(0), lbb(0), rpm(0), ma(0),vrr(0),vtt(0),vee(0),vxx(1.0)
        {
			
		cvl0 = 1.7161e-5;
		t0 = 273.16;
		ts = 110.4;
		rg = 287.0;
		cp = rg*1.4 / 0.4;
		prl = 0.72;
		prt = 0.9;
		cv1 = 7.1;
		cv2 = 5.0;
		kap = 0.41;
		cb1 = 0.1355;
		cb2 = 0.622;
		sigmav = 2.0 / 3.0;
		cw1 = cb1 / pow(kap, 2) + (1.0 + cb2) / sigmav;
		cw2 = 0.3;
		cw3 = 2.0;
		cr1 = 1.0;
		cr2 = 2.0;
		cr3 = 1.0;


			
	}
	~CDictionary()
	{
		Free(cg);
	}
	bool readConfig();

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
