/*---------------------------------------------------------------------------*\


\*---------------------------------------------------------------------------*/


#include "ArtificialVisc.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void CArtificialVisc::ResetViscousFlux()
{
	for (int k = 0; k<nz + 2; k++)
	for (int j = 0; j<ny + 2; j++)
	for (int i = 0; i<nx + 2; i++)
	{
		av1[k][j][i] = 0;
		av2[k][j][i] = 0;
		av3[k][j][i] = 0;
		av4[k][j][i] = 0;
		av5[k][j][i] = 0;
		pfield->q15[n][k][j][i] = pfield->q15[n][k][j][i] + pfield->p[k][j][i];
	}
}
void CArtificialVisc::ComputeViscousFlux_i()
{
	double em2, em4;	
	double flu1, flu2, flu3, flu4, flu5, ram;
	double* dp;
	Malloc(dp,nx + 2);
	dp[0] = 0;
	dp[nx + 1] = 0;
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	{
		for (int i = 1; i<nx + 1; i++)
			dp[i] = abs((pfield->p[k][j][i + 1] - 2.0*pfield->p[k][j][i] + pfield->p[k][j][i - 1]) / (pfield->p[k][j][i + 1] + 2.0*pfield->p[k][j][i] + pfield->p[k][j][i - 1]));

		for (int i = 0; i<nx + 1; i++)
		{
			ram = 0.50*(sri[k][j][i] + sri[k][j][i + 1]);
			em2 = a2*ram*max(dp[i], dp[i + 1]);
			em4 = a4*ram;
			em4 = max(em4 - em2, 0.0);

			flu1 = em2*(pfield->q11[n][k][j][i + 1] - pfield->q11[n][k][j][i]);
			flu2 = em2*(pfield->q12[n][k][j][i + 1] - pfield->q12[n][k][j][i]);
			flu3 = em2*(pfield->q13[n][k][j][i + 1] - pfield->q13[n][k][j][i]);
			flu4 = em2*(pfield->q14[n][k][j][i + 1] - pfield->q14[n][k][j][i]);
			flu5 = em2*(pfield->q15[n][k][j][i + 1] - pfield->q15[n][k][j][i]);

			if ((i >= 1) && (i <= nx - 1))
			{
				flu1 = flu1 + em4*(pfield->q11[n][k][j][i - 1] - 3.0*pfield->q11[n][k][j][i] + 3.0*pfield->q11[n][k][j][i + 1] - pfield->q11[n][k][j][i + 2]);
				flu2 = flu2 + em4*(pfield->q12[n][k][j][i - 1] - 3.0*pfield->q12[n][k][j][i] + 3.0*pfield->q12[n][k][j][i + 1] - pfield->q12[n][k][j][i + 2]);
				flu3 = flu3 + em4*(pfield->q13[n][k][j][i - 1] - 3.0*pfield->q13[n][k][j][i] + 3.0*pfield->q13[n][k][j][i + 1] - pfield->q13[n][k][j][i + 2]);
				flu4 = flu4 + em4*(pfield->q14[n][k][j][i - 1] - 3.0*pfield->q14[n][k][j][i] + 3.0*pfield->q14[n][k][j][i + 1] - pfield->q14[n][k][j][i + 2]);
				flu5 = flu5 + em4*(pfield->q15[n][k][j][i - 1] - 3.0*pfield->q15[n][k][j][i] + 3.0*pfield->q15[n][k][j][i + 1] - pfield->q15[n][k][j][i + 2]);
			}

			av1[k][j][i] = av1[k][j][i] + flu1;
			av2[k][j][i] = av2[k][j][i] + flu2;
			av3[k][j][i] = av3[k][j][i] + flu3;
			av4[k][j][i] = av4[k][j][i] + flu4;
			av5[k][j][i] = av5[k][j][i] + flu5;

			av1[k][j][i + 1] = av1[k][j][i + 1] - flu1;
			av2[k][j][i + 1] = av2[k][j][i + 1] - flu2;
			av3[k][j][i + 1] = av3[k][j][i + 1] - flu3;
			av4[k][j][i + 1] = av4[k][j][i + 1] - flu4;
			av5[k][j][i + 1] = av5[k][j][i + 1] - flu5;

		}

	}
	Free(dp);
}
void CArtificialVisc::ComputeViscousFlux_j()
{
	double em2, em4;	
	double flu1, flu2, flu3, flu4, flu5, ram;
	double* dp;
	Malloc(dp,ny + 2);
	dp[0] = 0;
	dp[ny + 1] = 0;

	for (int k = 1; k<nz + 1; k++)
		for (int i = 1; i<nx + 1; i++)
		{
			for (int j = 1; j<ny + 1; j++)
				dp[j] = abs((pfield->p[k][j + 1][i] - 2.0*pfield->p[k][j][i] + pfield->p[k][j - 1][i]) / (pfield->p[k][j + 1][i] + 2.0*pfield->p[k][j][i] + pfield->p[k][j - 1][i]));

			for (int j = 0; j<ny + 1; j++)
			{
				ram = 0.50*(srj[k][j][i] + srj[k][j + 1][i]);
				em2 = a2*ram*max(dp[j], dp[j + 1]);
				em4 = a4*ram;
				em4 = max(em4 - em2, 0.0);
				flu1 = em2*(pfield->q11[n][k][j + 1][i] - pfield->q11[n][k][j][i]);
				flu2 = em2*(pfield->q12[n][k][j + 1][i] - pfield->q12[n][k][j][i]);
				flu3 = em2*(pfield->q13[n][k][j + 1][i] - pfield->q13[n][k][j][i]);
				flu4 = em2*(pfield->q14[n][k][j + 1][i] - pfield->q14[n][k][j][i]);
				flu5 = em2*(pfield->q15[n][k][j + 1][i] - pfield->q15[n][k][j][i]);

				if ((j >= 1) && (j <= ny - 1))
				{
					flu1 = flu1 + em4*(pfield->q11[n][k][j - 1][i] - 3.0*pfield->q11[n][k][j][i] + 3.0*pfield->q11[n][k][j + 1][i] - pfield->q11[n][k][j + 2][i]);
					flu2 = flu2 + em4*(pfield->q12[n][k][j - 1][i] - 3.0*pfield->q12[n][k][j][i] + 3.0*pfield->q12[n][k][j + 1][i] - pfield->q12[n][k][j + 2][i]);
					flu3 = flu3 + em4*(pfield->q13[n][k][j - 1][i] - 3.0*pfield->q13[n][k][j][i] + 3.0*pfield->q13[n][k][j + 1][i] - pfield->q13[n][k][j + 2][i]);
					flu4 = flu4 + em4*(pfield->q14[n][k][j - 1][i] - 3.0*pfield->q14[n][k][j][i] + 3.0*pfield->q14[n][k][j + 1][i] - pfield->q14[n][k][j + 2][i]);
					flu5 = flu5 + em4*(pfield->q15[n][k][j - 1][i] - 3.0*pfield->q15[n][k][j][i] + 3.0*pfield->q15[n][k][j + 1][i] - pfield->q15[n][k][j + 2][i]);
				}

				av1[k][j][i] = av1[k][j][i] + flu1;
				av2[k][j][i] = av2[k][j][i] + flu2;
				av3[k][j][i] = av3[k][j][i] + flu3;
				av4[k][j][i] = av4[k][j][i] + flu4;
				av5[k][j][i] = av5[k][j][i] + flu5;

				av1[k][j + 1][i] = av1[k][j + 1][i] - flu1;
				av2[k][j + 1][i] = av2[k][j + 1][i] - flu2;
				av3[k][j + 1][i] = av3[k][j + 1][i] - flu3;
				av4[k][j + 1][i] = av4[k][j + 1][i] - flu4;
				av5[k][j + 1][i] = av5[k][j + 1][i] - flu5;
			}
		}
	Free(dp);
}
void CArtificialVisc::ComputeViscousFlux_k()
{
	double em2, em4;	
	double flu1, flu2, flu3, flu4, flu5, ram;

	double* dp;
	Malloc(dp,nz + 2);
	dp[0] = 0;
	dp[nz + 1] = 0;

	for (int i = 1; i<nx + 1; i++)
		for (int j = 1; j<ny + 1; j++)
		{
			for (int k = 1; k<nz + 1; k++)
				dp[k] = abs((pfield->p[k + 1][j][i] - 2.0*pfield->p[k][j][i] + pfield->p[k - 1][j][i]) / (pfield->p[k + 1][j][i] + 2.0*pfield->p[k][j][i] + pfield->p[k - 1][j][i]));

			for (int k = 0; k<nz + 1; k++)
			{
				ram = 0.50*(srk[k][j][i] + srk[k + 1][j][i]);
				em2 = a2*ram*max(dp[k], dp[k + 1]);
				em4 = a4*ram;
				em4 = max(em4 - em2, 0.0);

				flu1 = em2*(pfield->q11[n][k + 1][j][i] - pfield->q11[n][k][j][i]);
				flu2 = em2*(pfield->q12[n][k + 1][j][i] - pfield->q12[n][k][j][i]);
				flu3 = em2*(pfield->q13[n][k + 1][j][i] - pfield->q13[n][k][j][i]);
				flu4 = em2*(pfield->q14[n][k + 1][j][i] - pfield->q14[n][k][j][i]);
				flu5 = em2*(pfield->q15[n][k + 1][j][i] - pfield->q15[n][k][j][i]);

				if ((k >= 1) && (k <= nz - 1))
				{
					flu1 = flu1 + em4*(pfield->q11[n][k - 1][j][i] - 3.0*pfield->q11[n][k][j][i] + 3.0*pfield->q11[n][k + 1][j][i] - pfield->q11[n][k + 2][j][i]);
					flu2 = flu2 + em4*(pfield->q12[n][k - 1][j][i] - 3.0*pfield->q12[n][k][j][i] + 3.0*pfield->q12[n][k + 1][j][i] - pfield->q12[n][k + 2][j][i]);
					flu3 = flu3 + em4*(pfield->q13[n][k - 1][j][i] - 3.0*pfield->q13[n][k][j][i] + 3.0*pfield->q13[n][k + 1][j][i] - pfield->q13[n][k + 2][j][i]);
					flu4 = flu4 + em4*(pfield->q14[n][k - 1][j][i] - 3.0*pfield->q14[n][k][j][i] + 3.0*pfield->q14[n][k + 1][j][i] - pfield->q14[n][k + 2][j][i]);
					flu5 = flu5 + em4*(pfield->q15[n][k - 1][j][i] - 3.0*pfield->q15[n][k][j][i] + 3.0*pfield->q15[n][k + 1][j][i] - pfield->q15[n][k + 2][j][i]);
				}

				av1[k][j][i] = av1[k][j][i] + flu1;
				av2[k][j][i] = av2[k][j][i] + flu2;
				av3[k][j][i] = av3[k][j][i] + flu3;
				av4[k][j][i] = av4[k][j][i] + flu4;
				av5[k][j][i] = av5[k][j][i] + flu5;

				av1[k + 1][j][i] = av1[k + 1][j][i] - flu1;
				av2[k + 1][j][i] = av2[k + 1][j][i] - flu2;
				av3[k + 1][j][i] = av3[k + 1][j][i] - flu3;
				av4[k + 1][j][i] = av4[k + 1][j][i] - flu4;
				av5[k + 1][j][i] = av5[k + 1][j][i] - flu5;
			}
		}
	for (int k = 0; k<nz + 2; k++)
	for (int j = 0; j<ny + 2; j++)
	for (int i = 0; i<nx + 2; i++)
		pfield->q15[n][k][j][i] = pfield->q15[n][k][j][i] - pfield->p[k][j][i];
	Free(dp);
}
void CArtificialVisc::ComputeViscousFlux()
{
	ResetViscousFlux();
	ComputeViscousFlux_i();
	ComputeViscousFlux_j();
	ComputeViscousFlux_k();
}
void CArtificialVisc::ComputeViscousFlux_ci()
{
	double em2;
	double flu1, flu2, flu3, flu4, flu5;
	em2 = a2 / 1024.0 * sri[1][1][1];

	for (int k = 1; k<nz + 1; k++)
		for (int j = 1; j<ny + 1; j++)
			for (int i = 0; i<nx + 1; i++)
			{
				flu1 = em2*(pfield->q11[n][k][j][i + 1] - pfield->q11[n][k][j][i]);
				flu2 = em2*(pfield->q12[n][k][j][i + 1] - pfield->q12[n][k][j][i]);
				flu3 = em2*(pfield->q13[n][k][j][i + 1] - pfield->q13[n][k][j][i]);
				flu4 = em2*(pfield->q14[n][k][j][i + 1] - pfield->q14[n][k][j][i]);
				flu5 = em2*(pfield->q15[n][k][j][i + 1] - pfield->q15[n][k][j][i]);

				av1[k][j][i] = av1[k][j][i] + flu1;
				av2[k][j][i] = av2[k][j][i] + flu2;
				av3[k][j][i] = av3[k][j][i] + flu3;
				av4[k][j][i] = av4[k][j][i] + flu4;
				av5[k][j][i] = av5[k][j][i] + flu5;

				av1[k][j][i + 1] = av1[k][j][i + 1] - flu1;
				av2[k][j][i + 1] = av2[k][j][i + 1] - flu2;
				av3[k][j][i + 1] = av3[k][j][i + 1] - flu3;
				av4[k][j][i + 1] = av4[k][j][i + 1] - flu4;
				av5[k][j][i + 1] = av5[k][j][i + 1] - flu5;
			}
}
void CArtificialVisc::ComputeViscousFlux_cj()
{
	double em2;
	double flu1, flu2, flu3, flu4, flu5;
	em2 = a2 / 1024.0 * srj[1][1][1];

	for (int k = 1; k<nz + 1; k++)
		for (int j = 0; j<ny + 1; j++)
			for (int i = 1; i<nx + 1; i++)
			{
				flu1 = em2*(pfield->q11[n][k][j + 1][i] - pfield->q11[n][k][j][i]);
				flu2 = em2*(pfield->q12[n][k][j + 1][i] - pfield->q12[n][k][j][i]);
				flu3 = em2*(pfield->q13[n][k][j + 1][i] - pfield->q13[n][k][j][i]);
				flu4 = em2*(pfield->q14[n][k][j + 1][i] - pfield->q14[n][k][j][i]);
				flu5 = em2*(pfield->q15[n][k][j + 1][i] - pfield->q15[n][k][j][i]);

				av1[k][j][i] = av1[k][j][i] + flu1;
				av2[k][j][i] = av2[k][j][i] + flu2;
				av3[k][j][i] = av3[k][j][i] + flu3;
				av4[k][j][i] = av4[k][j][i] + flu4;
				av5[k][j][i] = av5[k][j][i] + flu5;

				av1[k][j + 1][i] = av1[k][j + 1][i] - flu1;
				av2[k][j + 1][i] = av2[k][j + 1][i] - flu2;
				av3[k][j + 1][i] = av3[k][j + 1][i] - flu3;
				av4[k][j + 1][i] = av4[k][j + 1][i] - flu4;
				av5[k][j + 1][i] = av5[k][j + 1][i] - flu5;

			}
}
void CArtificialVisc::ComputeViscousFlux_ck()
{
	double em2;
	double flu1, flu2, flu3, flu4, flu5;	
	em2 = a2 / 1024.0 * srk[1][1][1];

	for (int k = 0; k<nz + 1; k++)
		for (int j = 1; j<ny + 1; j++)
			for (int i = 1; i<nx + 1; i++)
			{
				flu1 = em2*(pfield->q11[n][k + 1][j][i] - pfield->q11[n][k][j][i]);
				flu2 = em2*(pfield->q12[n][k + 1][j][i] - pfield->q12[n][k][j][i]);
				flu3 = em2*(pfield->q13[n][k + 1][j][i] - pfield->q13[n][k][j][i]);
				flu4 = em2*(pfield->q14[n][k + 1][j][i] - pfield->q14[n][k][j][i]);
				flu5 = em2*(pfield->q15[n][k + 1][j][i] - pfield->q15[n][k][j][i]);

				av1[k][j][i] = av1[k][j][i] + flu1;
				av2[k][j][i] = av2[k][j][i] + flu2;
				av3[k][j][i] = av3[k][j][i] + flu3;
				av4[k][j][i] = av4[k][j][i] + flu4;
				av5[k][j][i] = av5[k][j][i] + flu5;

				av1[k + 1][j][i] = av1[k + 1][j][i] - flu1;
				av2[k + 1][j][i] = av2[k + 1][j][i] - flu2;
				av3[k + 1][j][i] = av3[k + 1][j][i] - flu3;
				av4[k + 1][j][i] = av4[k + 1][j][i] - flu4;
				av5[k + 1][j][i] = av5[k + 1][j][i] - flu5;
			}

	for (int k = 0; k<nz + 2; k++)
		for (int j = 0; j<ny + 2; j++)
			for (int i = 0; i<nx + 2; i++)
			{
				pfield->q15[n][k][j][i] = pfield->q15[n][k][j][i] - pfield->p[k][j][i];
			}
}
void CArtificialVisc::ComputeViscousFlux_c()
{
	ResetViscousFlux();
	ComputeViscousFlux_ci();
	ComputeViscousFlux_cj();
	ComputeViscousFlux_ck();
}
void CArtificialVisc::ResetViscousFlux_SA()
{
	for (int k = 0; k<nz + 2; k++)
	for (int j = 0; j<ny + 2; j++)
	for (int i = 0; i<nx + 2; i++)
	{
		av6[k][j][i] = 0;
	}
}
void CArtificialVisc::ComputeViscousFlux_SAi()
{
	double em2, em4;
	double flu6, ram;

	double* dp;
	Malloc(dp,nx + 2);
	dp[0] = 0;
	dp[nx + 1] = 0;

	for (int k = 1; k<nz + 1; k++)
		for (int j = 1; j<ny + 1; j++)
		{
			for (int i = 1; i<nx + 1; i++)
				dp[i] = abs((pfield->p[k][j][i + 1] - 2.0*pfield->p[k][j][i] + pfield->p[k][j][i - 1]) / (pfield->p[k][j][i + 1] + 2.0*pfield->p[k][j][i] + pfield->p[k][j][i - 1]));

			for (int i = 0; i<nx + 1; i++)
			{
				ram = 0.50*(sri[k][j][i] + sri[k][j][i + 1]);
				em2 = a2*ram*max(dp[i], dp[i + 1]);
				em4 = a4*ram;
				em4 = max(em4 - em2, 0.0);
				flu6 = em2*(pfield->q16[n][k][j][i + 1] - pfield->q16[n][k][j][i]);

				if ((i >= 1) && (i <= nx - 1))
					flu6 = flu6 + em4*(pfield->q16[n][k][j][i - 1] - 3.0*pfield->q16[n][k][j][i] + 3.0*pfield->q16[n][k][j][i + 1] - pfield->q16[n][k][j][i + 2]);

				av6[k][j][i] = av6[k][j][i] + flu6;
				av6[k][j][i + 1] = av6[k][j][i + 1] - flu6;
			}
		}
	Free(dp);
}
void CArtificialVisc::ComputeViscousFlux_SAj()
{
	double em2, em4;
	double flu6, ram;

	double* dp;
	Malloc(dp,ny + 2);	
	dp[0] = 0;
	dp[ny + 1] = 0;

	for (int k = 1; k<nz + 1; k++)
		for (int i = 1; i<nx + 1; i++)
		{
			for (int j = 1; j<ny + 1; j++)
			{
				dp[j] = abs((pfield->p[k][j + 1][i] - 2.0*pfield->p[k][j][i] + pfield->p[k][j - 1][i]) / (pfield->p[k][j + 1][i] + 2.0*pfield->p[k][j][i] + pfield->p[k][j - 1][i]));
			}
			for (int j = 0; j<ny + 1; j++)
			{
				ram = 0.50*(srj[k][j][i] + srj[k][j + 1][i]);
				em2 = a2*ram*max(dp[j], dp[j + 1]);
				em4 = a4*ram;
				em4 = max(em4 - em2, 0.0);
				flu6 = em2*(pfield->q16[n][k][j + 1][i] - pfield->q16[n][k][j][i]);

				if ((j >= 1) && (j <= ny - 1))
					flu6 = flu6 + em4*(pfield->q16[n][k][j - 1][i] - 3.0*pfield->q16[n][k][j][i] + 3.0*pfield->q16[n][k][j + 1][i] - pfield->q16[n][k][j + 2][i]);

				av6[k][j][i] = av6[k][j][i] + flu6;
				av6[k][j + 1][i] = av6[k][j + 1][i] - flu6;
			}
		}
	Free(dp);
}
void CArtificialVisc::ComputeViscousFlux_SAk()
{	
	double em2, em4;
	double flu6, ram;

	double* dp;
	Malloc(dp,nz + 2);	
	dp[0] = 0;
	dp[nz + 1] = 0;

	for (int i = 1; i<nx + 1; i++)
		for (int j = 1; j<ny + 1; j++)
		{
			for (int k = 1; k<nz + 1; k++)
				dp[k] = abs((pfield->p[k + 1][j][i] - 2.0*pfield->p[k][j][i] + pfield->p[k - 1][j][i]) / (pfield->p[k + 1][j][i] + 2.0*pfield->p[k][j][i] + pfield->p[k - 1][j][i]));

			for (int k = 0; k<nz + 1; k++)
			{
				ram = 0.50*(srk[k][j][i] + srk[k + 1][j][i]);
				em2 = a2*ram*max(dp[k], dp[k + 1]);
				em4 = a4*ram;
				em4 = max(em4 - em2, 0.0);
				flu6 = em2*(pfield->q16[n][k + 1][j][i] - pfield->q16[n][k][j][i]);

				if ((k >= 1) && (k <= nz - 1))
				{
					flu6 = flu6 + em4*(pfield->q16[n][k - 1][j][i] - 3.0*pfield->q16[n][k][j][i] + 3.0*pfield->q16[n][k + 1][j][i] - pfield->q16[n][k + 2][j][i]);
				}

				av6[k][j][i] = av6[k][j][i] + flu6;
				av6[k + 1][j][i] = av6[k + 1][j][i] - flu6;
			}
		}
	Free(dp);
}
void CArtificialVisc::ComputeViscousFlux_SA()
{
	ResetViscousFlux_SA();
	ComputeViscousFlux_SAi();
	ComputeViscousFlux_SAj();
	ComputeViscousFlux_SAk();
}

bool CArtificialVisc::step_artificial()
{

	SetLocalTimestep( );
	ComputeViscousFlux();
	return true;
}


bool CArtificialVisc::step_artificial_c()
{

	SetLocalTimestep( );
	ComputeViscousFlux_c();
	return true;
}

bool CArtificialVisc::step_artificial_SA()
{

	SetLocalTimestep( );
	ComputeViscousFlux_SA();
	return true;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



