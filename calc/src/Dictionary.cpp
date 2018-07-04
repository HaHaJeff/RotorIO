/*---------------------------------------------------------------------------*\


\*---------------------------------------------------------------------------*/


#include "Dictionary.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool CDictionary::readConfig()
{
	ifstream fin(file1.c_str());
		fin >> nt >> ng >> beta1 >>beta2 >> cfl >> a2 >> a4 >> ht >> pt >> pb1 >> c2 >> rmsm0;

		
		vtt = vxx*beta1;
		vrr = vxx*beta2;
		vee = sqrt(pow(vxx, 2) + pow(vrr, 2) + pow(vtt, 2));
		
		Malloc(cg,ng+1);
		forAll(1,ng + 1,i)
		{
			fin >> (cg)[i];
		}

		fin >> nxm >> nym >> nzm >> ibm >>itm >> jbm >> jtm >> lbb >> rpm >> ma;


		rpm = rpm * pi / 30.0;
		fin.close(); 
		

	return rpm;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



