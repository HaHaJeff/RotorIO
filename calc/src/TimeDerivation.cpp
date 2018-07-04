/*---------------------------------------------------------------------------*\


\*---------------------------------------------------------------------------*/


#include "TimeDerivation.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void CTimeDerivation::SetTimespectrum()
{	
	double 	period = 2.0*pi / pdict->rpm;
	pfield->nt = pdict->nt;
	if (pfield->nt % 2 == 1)
	{
		for (int i = 1; i<pfield->nt + 1; i++)
		for (int j = 1; j<pfield->nt + 1; j++)
		{
			if (i != j)
			dm[i][j] = (pi / period*pow(-1.0, (i - j)) / sin(pi*double((i - j)) / double(pfield->nt)));
		}
	}
	else
	{
		for (int i = 1; i<pfield->nt + 1; i++)
		for (int j = 1; j<pfield->nt + 1; j++)
		{
			if (i != j)
			dm[i][j] = (pi / period*pow(-1.0, (i - j)) / tan(pi*double((i - j)) / double(pfield->nt)));
		}
					
	}
}

bool CTimeDerivation::timedrivationCmpute()
{


	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
	{
					
		ts1[k][j][i] = 0;
		ts2[k][j][i] = 0;
		ts3[k][j][i] = 0;
		ts4[k][j][i] = 0;
		ts5[k][j][i] = 0;
	}

	for (int nn = 1; nn<pfield->nt + 1; nn++)
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
	{

		ts1[k][j][i] = ts1[k][j][i] + dm[nn][n] * (pfield->q11)[nn][k][j][i];
		ts2[k][j][i] = ts2[k][j][i] + dm[nn][n] * (pfield->q12)[nn][k][j][i];
		ts3[k][j][i] = ts3[k][j][i] + dm[nn][n] * (pfield->q13)[nn][k][j][i];
		ts4[k][j][i] = ts4[k][j][i] + dm[nn][n] * (pfield->q14)[nn][k][j][i];
		ts5[k][j][i] = ts5[k][j][i] + dm[nn][n] * (pfield->q15)[nn][k][j][i];

	}

	return true;					
	
}

bool CTimeDerivation::timedrivationCmputeSA()
{


	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
	{
					
		ts6[k][j][i] = 0;
	}
				
	for (int nn = 1; nn<pfield->nt + 1; nn++)
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		ts6[k][j][i] = ts6[k][j][i] + dm[nn][n] * (pfield->q16)[nn][k][j][i];
	}		
	return 1;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



