/*---------------------------------------------------------------------------*\


\*---------------------------------------------------------------------------*/


#include "GeoField.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

bool CGeoField::SetPhysicalBoundary(int nng0,double vxx, double vrr,double vtt,double vee,double ht,double pt)
{
	//int nx = (*(mesh->nnx))[nng0];
	int ny = ((mesh->nny))[nng0];
	int nz = ((mesh->nnz))[nng0];
	double z1,y1,sir,cor,rr;
	//if(nng0==1)
	//mesh->SetlayerMesh(nng0);	
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	{
		y1 = 0.25*((mesh->y)[k][j][1] + (mesh->y)[k][j + 1][1] + (mesh->y)[k + 1][j][1] + (mesh->y)[k + 1][j + 1][1]);
		z1 = 0.25*((mesh->z)[k][j][1] + (mesh->z)[k][j + 1][1] + (mesh->z)[k + 1][j][1] + (mesh->z)[k + 1][j + 1][1]);
		rr = sqrt(y1*y1 + z1*z1);
		sir = z1 / rr;
		cor = y1 / rr;
		for (int i = 1; i<nt + 1; i++)
		{
			(betaxn)[i][nng0][k][j] = vxx / vee;
			(betayn)[i][nng0][k][j] = (vrr*cor - vtt*sir) / vee;
			(betazn)[i][nng0][k][j] = (vrr*sir + vtt*cor) / vee;
		}			
	}		

	for (int k = 1; k<nt + 1; k++)
	for (int j = 1; j<nz + 1; j++)
	for (int i = 1; i<ny + 1; i++)
	{
		(hatn)[k][nng0][j][i] = ht;
		(petn)[k][nng0][j][i] = pt;
	}
	return vee;
}
bool CGeoField::GetlayerPhysicalBoundary(int n,int nnng)
{
	int ny = ((mesh->nny))[nnng];
	int nz = ((mesh->nnz))[nnng];
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	{
		betax[k][j] = betaxn[n][nnng][k][j];
		betay[k][j] = betayn[n][nnng][k][j];
		betaz[k][j] = betazn[n][nnng][k][j];
		hat[k][j] = hatn[n][nnng][k][j];
		pet[k][j] = petn[n][nnng][k][j];
	}
	return 1;
}
bool CGeoField::FieldInitialization(int nng0,double ma,double cvl0)
{
	int nx = ((mesh->nnx))[nng0];
	int ny = ((mesh->nny))[nng0];
	int nz = ((mesh->nnz))[nng0];
	double a,rout,pp,dim,en;
	for (int n = 1; n<nt + 1; n++)
	for (int k = 1; k<nz + 1; k++)							
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		rout = 3.5*(petn)[n][nng0][k][j] / (hatn)[n][nng0][k][j];
		pp = (petn)[n][nng0][k][j] * pow((1.0 + 0.2*ma*ma), -3.5);	
		dim = rout*pow((1.0 + 0.2*ma*ma), -2.5);
		a = sqrt(1.4*pp / dim);
		en = 2.5*pp + 0.5*dim*ma*ma*a*a;
		(q11)[n][k][j][i] = dim;
		(q12)[n][k][j][i] = dim*a*ma*(betaxn)[n][nng0][k][j];
		(q13)[n][k][j][i] = dim*a*ma*(betayn)[n][nng0][k][j];
		(q14)[n][k][j][i] = dim*a*ma*(betazn)[n][nng0][k][j];
		(q15)[n][k][j][i] = en;
		(q16)[n][k][j][i] = 200 * cvl0;
	}
	return dim;
}
void CGeoField::SetPhysicalBoundary(int nng0, double vxx, double vrr, double vtt, double vee, double ht, double pt, CMultiGrid *p)
{
	int ny = ((mesh->nny))[nng0];
	int nz = ((mesh->nnz))[nng0];
	double z1, y1, sir, cor, rr;
	//if(nng0==1)
	//mesh->SetlayerMesh(nng0);	
	for (int k = 1; k<nz + 1; k++)
		for (int j = 1; j<ny + 1; j++)
		{
			y1 = p->yy2[nng0][k][j][1];
			z1 = p->zz2[nng0][k][j][1];
			rr = sqrt(y1*y1 + z1*z1);
			sir = z1 / rr;
			cor = y1 / rr;
			for (int i = 1; i<nt + 1; i++)
			{
				(betaxn)[i][nng0][k][j] = vxx / vee;
				(betayn)[i][nng0][k][j] = (vrr*cor - vtt*sir) / vee;
				(betazn)[i][nng0][k][j] = (vrr*sir + vtt*cor) / vee;
			}
		}

	for (int k = 1; k<nt + 1; k++)
		for (int j = 1; j<nz + 1; j++)
			for (int i = 1; i<ny + 1; i++)
			{
				(hatn)[k][nng0][j][i] = ht;
				(petn)[k][nng0][j][i] = pt;
			}
}
	
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



