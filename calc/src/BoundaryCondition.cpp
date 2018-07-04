/*---------------------------------------------------------------------------*\


\*---------------------------------------------------------------------------*/


#include "BoundaryCondition.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void CBoundaryCondition::setBC_in()
{	
	double cvl,cvu,cvlt,vx,vy,vz,en,pp,t1,dim;
	double tem, tur1, tur2, tur3, fv1,sir,cor;
	double qq2, sxn, syn, szn, ds, a;
	double u, v, w, uabs, unorm, rinv, c02, dis, cb, cosa, hb, cc02;
	for (int i = 0; i < ny + 1; i++)
	{
		hr[i] = 0;
		hv[i] = 0;
		hd[i] = 0;
		hp[i] = 0;
	}
	

	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	{
				
		dim = (pfield->q11)[n][k][j][1];
		vx = (pfield->q12)[n][k][j][1] / dim;
		vy = (pfield->q13)[n][k][j][1] / dim;
		vz = (pfield->q14)[n][k][j][1] / dim;
		qq2 = vx*vx + vy*vy + vz*vz;
		en = (pfield->q15)[n][k][j][1];
		pp = 0.4*(en - 0.5*dim*qq2);
		t1 = pp / (dim*(pdict->rg));
		cvlt = (pdict->cvl0)*pow(t1 / (pdict->t0), 1.5)*((pdict->t0) + (pdict->ts)) / (t1 + (pdict->ts));
		cvl = cvlt / dim;
		cvu = (pdict->c2)*cvl;
		tur1 = pow(10.0, -4);
		tur2 = pow(10.0, -6);
		tur3 = pow(10.0, -6);
		while (abs((tur1 - tur3) / tur3)>pow(10.0, -6))
		{
			tem = tur1 / cvl;
			fv1 = 1.0 / (1.0 + pow((pdict->cv1) / tem, 3));
			tur2 = cvu / fv1;
			tur3 = tur1;
			tur1 = tur2;
		}		
		turi[k][j] = tur2;

		ds = sqrt((pgrid->s2x)[k][j][1] * (pgrid->s2x)[k][j][1] + (pgrid->s2y)[k][j][1] * (pgrid->s2y)[k][j][1] + (pgrid->s2z)[k][j][1] * (pgrid->s2z)[k][j][1]);
		sxn = (pgrid->s2x)[k][j][1] / ds;
		syn = (pgrid->s2y)[k][j][1] / ds;
		szn = (pgrid->s2z)[k][j][1] / ds;
		u = (pfield->pvx)[k][j][1];
		v = (pfield->pvy)[k][j][1];
		w = (pfield->pvz)[k][j][1];
		uabs = sqrt(u*u + v*v + w*w);
		unorm = u*sxn + v*syn + w*szn;
		a = sqrt(1.4*(pfield->p)[k][j][1] / (pfield->q11)[n][k][j][1]);
		
		/*cosa=(uabs<pow(10.0, -20))?1.0:(-unorm / uabs);*/
		if (uabs<pow(10.0, -20))
			cosa = 1.0;
		else
			cosa = -unorm / uabs;
		rinv = unorm - 5 * a;
		c02 = a*a + 0.2*uabs*uabs;
		dis = (0.4*cosa*cosa + 2.0)*c02 / (0.4*rinv*rinv) - 0.2;
		if (dis<0)
			dis = pow(10.0, -20);
		cb = -rinv*(0.4 / (0.4*cosa*cosa + 2.0))*(1.0 + cosa*sqrt(dis));
		cc02 = min(cb*cb / c02, 1.0);
		hb = (pfield->hat)[k][j] * cc02;

		(pfield->t)[k][j][0] = hb / (pdict->cp);
		(pfield->p)[k][j][0] = (pfield->pet)[k][j] * pow(cc02, 3.5);
		(pfield->q11)[n][k][j][0] = (pfield->p)[k][j][0] / ((pfield->t)[k][j][0] * (pdict->rg));
		uabs = 2.0*((pfield->hat)[k][j] - hb);
		(pfield->q15)[n][k][j][0] = 2.5*(pfield->p)[k][j][0] + 0.5*(pfield->q11)[n][k][j][0] * uabs;
		(pfield->pvx)[k][j][0] = sqrt(uabs)*(pfield->betax)[k][j];
		(pfield->pvy)[k][j][0] = sqrt(uabs)*(pfield->betay)[k][j];
		(pfield->pvz)[k][j][0] = sqrt(uabs)*(pfield->betaz)[k][j];
		(pfield->q12)[n][k][j][0] = (pfield->q11)[n][k][j][0] * (pfield->pvx)[k][j][0];
		(pfield->q13)[n][k][j][0] = (pfield->q11)[n][k][j][0] * (pfield->pvy)[k][j][0];
		(pfield->q14)[n][k][j][0] = (pfield->q11)[n][k][j][0] * (pfield->pvz)[k][j][0];
		(pfield->q16)[n][k][j][0] = turi[k][j] * (pfield->q11)[n][k][j][0];
	}
}

void CBoundaryCondition::setBC_out()
{
	double ss, sr, sv, sd, s1,y1,z1, ve, v2, dr,dim,rr ,deltp, rrhoc;
	double qq2, sxn, syn, szn, ds, a;
		
     
       	
	for (int j = 1; j<ny + 1; j++)
	{
		ss = 0;
		sr = 0;
		sv = 0;
		sd = 0;
	for (int k = 1; k<nz + 1; k++)
	{
		s1 = sqrt((pgrid->s2x)[k][j][nx + 1] * (pgrid->s2x)[k][j][nx + 1] + (pgrid->s2y)[k][j][nx + 1] * (pgrid->s2y)[k][j][nx + 1] + (pgrid->s2z)[k][j][nx + 1] * (pgrid->s2z)[k][j][nx + 1]);
		y1 = (pgrid->yy02)[k][j][nx + 1];
		z1 = (pgrid->zz02)[k][j][nx + 1];
		rr = sqrt(y1*y1 + z1*z1);
		dim = (pfield->q11)[n][k][j][nx];
		ve = (pfield->vth)[k][j][nx];
		ss = ss + s1;
		sr = sr + rr*s1;
		sd = sd + dim*s1;
		sv = sv + ve*s1;
	}
			
		hr[j] = sr / ss;
		hv[j] = sv / ss;
		hd[j] = sd / ss;

	}
		
	hp[1] = (pdict->pb1);
		
	for (int j = 2; j<ny + 1; j++)
	{
		dim = 0.5*(hd[j - 1] + hd[j]);
		v2 = 0.5*(hv[j - 1] + hv[j]);
		dr = hr[j] - hr[j - 1];
		rr = 0.5*(hr[j] + hr[j - 1]);
		hp[j] = hp[j - 1] + dim*v2*v2 / rr*dr;
	}
		
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	{
		peb[k][j] = hp[j];
	}
			
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	{
		ds = sqrt((pgrid->s2x)[k][j][nx + 1] * (pgrid->s2x)[k][j][nx + 1] + (pgrid->s2y)[k][j][nx + 1] * (pgrid->s2y)[k][j][nx + 1] + (pgrid->s2z)[k][j][nx + 1] * (pgrid->s2z)[k][j][nx + 1]);
		sxn = -(pgrid->s2x)[k][j][nx + 1] / ds;
		syn = -(pgrid->s2y)[k][j][nx + 1] / ds;
		szn = -(pgrid->s2z)[k][j][nx + 1] / ds;
		a = sqrt(1.4*(pfield->p)[k][j][nx] / (pfield->q11)[n][k][j][nx]);
		rrhoc = 1.0 / ((pfield->q11)[n][k][j][nx] * a);
		deltp = (pfield->p)[k][j][nx] - peb[k][j];
		(pfield->p)[k][j][nx + 1] = peb[k][j];
		(pfield->q11)[n][k][j][nx + 1] = (pfield->q11)[n][k][j][nx] - deltp / (a*a);
		(pfield->pvx)[k][j][nx + 1] = (pfield->pvx)[k][j][nx] + sxn*deltp*rrhoc;
		(pfield->pvy)[k][j][nx + 1] = (pfield->pvy)[k][j][nx] + syn*deltp*rrhoc;
		(pfield->pvz)[k][j][nx + 1] = (pfield->pvz)[k][j][nx] + szn*deltp*rrhoc;
		(pfield->q12)[n][k][j][nx + 1] = (pfield->q11)[n][k][j][nx + 1] * (pfield->pvx)[k][j][nx + 1];
		(pfield->q13)[n][k][j][nx + 1] = (pfield->q11)[n][k][j][nx + 1] * (pfield->pvy)[k][j][nx + 1];
		(pfield->q14)[n][k][j][nx + 1] = (pfield->q11)[n][k][j][nx + 1] * (pfield->pvz)[k][j][nx + 1];
		qq2 = (pfield->pvx)[k][j][nx + 1] * (pfield->pvx)[k][j][nx + 1] + (pfield->pvy)[k][j][nx + 1] * (pfield->pvy)[k][j][nx + 1] + (pfield->pvz)[k][j][nx + 1] * (pfield->pvz)[k][j][nx + 1];
		(pfield->q15)[n][k][j][nx + 1] = 2.5*(pfield->p)[k][j][nx + 1] + 0.5*(pfield->q11)[n][k][j][nx + 1] * qq2;
		(pfield->q16)[n][k][j][nx + 1] = (pfield->q16)[n][k][j][nx];
		(pfield->t)[k][j][nx + 1]=(pfield->p)[k][j][nx + 1] / ((pfield->q11)[n][k][j][nx + 1] * (pdict->rg));
	}
}
void CBoundaryCondition::setBC_root_top()
{	
	double qq2;
	

	for (int k = 1; k<nz + 1; k++)
	for (int i = 1; i<nx + 1; i++)
	{
		//*****叶根位置****
						
		(pfield->p)[k][0][i] = (pfield->p)[k][1][i];
		(pfield->t)[k][0][i] = (pfield->t)[k][1][i];
		(pfield->q11)[n][k][0][i] = (pfield->q11)[n][k][1][i];
		(pfield->pvx)[k][0][i] = -(pfield->pvx)[k][1][i];
		(pfield->pvy)[k][0][i] = -2.0*(pdict->rpm)*(pgrid->zz03)[k][1][i] - (pfield->pvy)[k][1][i];
		(pfield->pvz)[k][0][i] = 2.0*(pdict->rpm)*(pgrid->yy03)[k][1][i] - (pfield->pvz)[k][1][i];
		(pfield->q12)[n][k][0][i] = (pfield->q11)[n][k][0][i] * (pfield->pvx)[k][0][i];
		(pfield->q13)[n][k][0][i] = (pfield->q11)[n][k][0][i] * (pfield->pvy)[k][0][i];
		(pfield->q14)[n][k][0][i] = (pfield->q11)[n][k][0][i] * (pfield->pvz)[k][0][i];
		qq2 = (pfield->pvx)[k][0][i] * (pfield->pvx)[k][0][i] + (pfield->pvy)[k][0][i] * (pfield->pvy)[k][0][i] + (pfield->pvz)[k][0][i] * (pfield->pvz)[k][0][i];
		(pfield->q15)[n][k][0][i] = 2.5*(pfield->p)[k][0][i] + 0.5*(pfield->q11)[n][k][0][i] * qq2;
		(pfield->q16)[n][k][0][i] = -(pfield->q16)[n][k][1][i];
						
		//*****叶顶位置****
						
		(pfield->p)[k][ny + 1][i] = (pfield->p)[k][ny][i];
		(pfield->t)[k][ny + 1][i] = (pfield->t)[k][ny][i];
		(pfield->q11)[n][k][ny + 1][i] = (pfield->q11)[n][k][ny][i];
		(pfield->pvx)[k][ny + 1][i] = -(pfield->pvx)[k][ny][i];
		(pfield->pvy)[k][ny + 1][i] = -(pfield->pvy)[k][ny][i];
		(pfield->pvz)[k][ny + 1][i] = -(pfield->pvz)[k][ny][i];
		(pfield->q12)[n][k][ny + 1][i] = (pfield->q11)[n][k][ny + 1][i] * (pfield->pvx)[k][ny + 1][i];
		(pfield->q13)[n][k][ny + 1][i] = (pfield->q11)[n][k][ny + 1][i] * (pfield->pvy)[k][ny + 1][i];
		(pfield->q14)[n][k][ny + 1][i] = (pfield->q11)[n][k][ny + 1][i] * (pfield->pvz)[k][ny + 1][i];
		qq2 = (pfield->pvx)[k][ny + 1][i] * (pfield->pvx)[k][ny + 1][i] + (pfield->pvy)[k][ny + 1][i] * (pfield->pvy)[k][ny + 1][i] + (pfield->pvz)[k][ny + 1][i] * (pfield->pvz)[k][ny + 1][i];
		(pfield->q15)[n][k][ny + 1][i] = 2.5*(pfield->p)[k][ny + 1][i] + 0.5*(pfield->q11)[n][k][ny + 1][i] * qq2;
		(pfield->q16)[n][k][ny + 1][i] = -(pfield->q16)[n][k][ny][i];
	}
}

void CBoundaryCondition::setBC_circumference()
{

     
	//***************前后固壁边界z的计算***************
	for (int j = jb; j<jt + 1; j++)
	for (int i = ib; i<it + 1; i++)
	{
		(pfield->q11)[n][0][j][i] = (pfield->q11)[n][1][j][i];
		(pfield->pvx)[0][j][i] = -(pfield->pvx)[1][j][i];
		(pfield->pvy)[0][j][i] = -2.0*(pdict->rpm)*(pgrid->zz01)[1][j][i] - (pfield->pvy)[1][j][i];
		(pfield->pvz)[0][j][i] = 2.0*(pdict->rpm)*(pgrid->yy01)[1][j][i] - (pfield->pvz)[1][j][i];
		(pfield->p)[0][j][i] = (pfield->p)[1][j][i];
		(pfield->q16)[n][0][j][i] = -(pfield->q16)[n][1][j][i];
		(pfield->q11)[n][nz + 1][j][i] = (pfield->q11)[n][nz][j][i];
		(pfield->pvx)[nz + 1][j][i] = -(pfield->pvx)[nz][j][i];
		(pfield->pvy)[nz + 1][j][i] = -2.0*(pdict->rpm)*(pgrid->zz01)[nz + 1][j][i] - (pfield->pvy)[nz][j][i];
		(pfield->pvz)[nz + 1][j][i] = 2.0*(pdict->rpm)*(pgrid->yy01)[nz + 1][j][i] - (pfield->pvz)[nz][j][i];
		(pfield->p)[nz + 1][j][i] = (pfield->p)[nz][j][i];
		(pfield->q16)[n][nz + 1][j][i] = -(pfield->q16)[n][nz][j][i];
	}
}

void CBoundaryCondition::communication()
{
	myidl = (myid +numprocs- 1)%numprocs;
	myidr = (myid + 1)%numprocs;

	for (int i = 1; i < ny+1; i++)
	{
		MPI_Sendrecv(&(pfield->q11)[n][1][i][1], ib-1, MPI_DOUBLE, myidl, 20, &(pfield->q11)[n][nz+1][i][1] ,ib-1, MPI_DOUBLE, myidr, 20, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&(pfield->pvx)[1][i][1], ib-1, MPI_DOUBLE, myidl, 21, &(pfield->pvx)[nz+1][i][1] ,ib-1, MPI_DOUBLE, myidr, 21, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&(pfield->vth)[1][i][1], ib-1, MPI_DOUBLE, myidl, 22, &(pfield->vth)[nz+1][i][1] ,ib-1, MPI_DOUBLE, myidr, 22, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&(pfield->vre)[1][i][1], ib-1, MPI_DOUBLE, myidl, 23, &(pfield->vre)[nz+1][i][1] ,ib-1, MPI_DOUBLE, myidr, 23, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&(pfield->p)[1][i][1], ib-1, MPI_DOUBLE, myidl, 24, &(pfield->p)[nz+1][i][1] ,ib-1, MPI_DOUBLE, myidr, 24, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&(pfield->q16)[n][1][i][1], ib-1, MPI_DOUBLE, myidl, 25, &(pfield->q16)[n][nz+1][i][1] ,ib-1, MPI_DOUBLE, myidr, 25, MPI_COMM_WORLD, &status);
	}

	for (int i = 1; i < ny+1; i++)
	{
		MPI_Sendrecv(&(pfield->q11)[n][1][i][it+1], nx-it, MPI_DOUBLE, myidl, 30, &(pfield->q11)[n][nz+1][i][it+1], nx-it, MPI_DOUBLE, myidr, 30, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&(pfield->pvx)[1][i][it+1], nx-it, MPI_DOUBLE, myidl, 31, &(pfield->pvx)[nz+1][i][it+1], nx-it, MPI_DOUBLE, myidr, 31, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&(pfield->vth)[1][i][it+1], nx-it, MPI_DOUBLE, myidl, 32, &(pfield->vth)[nz+1][i][it+1], nx-it, MPI_DOUBLE, myidr, 32, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&(pfield->vre)[1][i][it+1], nx-it, MPI_DOUBLE, myidl, 33, &(pfield->vre)[nz+1][i][it+1], nx-it, MPI_DOUBLE, myidr, 33, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&(pfield->p)[1][i][it+1], nx-it, MPI_DOUBLE, myidl, 34, &(pfield->p)[nz+1][i][it+1], nx-it, MPI_DOUBLE, myidr, 34, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&(pfield->q16)[n][1][i][it+1], nx-it, MPI_DOUBLE, myidl, 35, &(pfield->q16)[n][nz+1][i][it+1], nx-it, MPI_DOUBLE, myidr, 35, MPI_COMM_WORLD, &status);
	}

	for(int i=jt+1 ;i<ny+1;i++)
	{
		MPI_Sendrecv(&(pfield->q11)[n][1][i][ib], it-ib+1, MPI_DOUBLE, myidl, 40, &(pfield->q11)[n][nz+1][i][ib] ,it-ib+1, MPI_DOUBLE, myidr, 40, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&(pfield->pvx)[1][i][ib], it-ib+1, MPI_DOUBLE, myidl, 41, &(pfield->pvx)[nz+1][i][ib] ,it-ib+1, MPI_DOUBLE, myidr, 41, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&(pfield->vth)[1][i][ib], it-ib+1, MPI_DOUBLE, myidl, 42, &(pfield->vth)[nz+1][i][ib] ,it-ib+1, MPI_DOUBLE, myidr, 42, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&(pfield->vre)[1][i][ib], it-ib+1, MPI_DOUBLE, myidl, 43, &(pfield->vre)[nz+1][i][ib] ,it-ib+1, MPI_DOUBLE, myidr, 43, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&(pfield->p)[1][i][ib], it-ib+1, MPI_DOUBLE, myidl, 44, &(pfield->p)[nz+1][i][ib] ,it-ib+1, MPI_DOUBLE, myidr, 44, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&(pfield->q16)[n][1][i][ib], it-ib+1, MPI_DOUBLE, myidl, 45, &(pfield->q16)[n][nz+1][i][ib] ,it-ib+1, MPI_DOUBLE, myidr, 45, MPI_COMM_WORLD, &status);
	}


	for(int i=1;i<ny+1;i++)
	{
		MPI_Sendrecv(&(pfield->q11)[n][nz][i][1], ib-1, MPI_DOUBLE, myidr, 50, &(pfield->q11)[n][0][i][1] ,ib-1, MPI_DOUBLE, myidl, 50, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&(pfield->pvx)[nz][i][1], ib-1, MPI_DOUBLE, myidr, 51, &(pfield->pvx)[0][i][1] ,ib-1, MPI_DOUBLE, myidl, 51, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&(pfield->vth)[nz][i][1], ib-1, MPI_DOUBLE, myidr, 52, &(pfield->vth)[0][i][1] ,ib-1, MPI_DOUBLE, myidl, 52, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&(pfield->vre)[nz][i][1], ib-1, MPI_DOUBLE, myidr, 53, &(pfield->vre)[0][i][1] ,ib-1, MPI_DOUBLE, myidl, 53, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&(pfield->p)[nz][i][1], ib-1, MPI_DOUBLE, myidr, 54, &(pfield->p)[0][i][1] ,ib-1, MPI_DOUBLE, myidl, 54, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&(pfield->q16)[n][nz][i][1], ib-1, MPI_DOUBLE, myidr, 55, &(pfield->q16)[n][0][i][1] ,ib-1, MPI_DOUBLE, myidl, 55, MPI_COMM_WORLD, &status);
	}
					 
	for(int i=1;i<ny+1;i++)
	{
		MPI_Sendrecv(&(pfield->q11)[n][nz][i][it+1], nx-it, MPI_DOUBLE, myidr, 60, &(pfield->q11)[n][0][i][it+1] ,nx-it, MPI_DOUBLE, myidl, 60, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&(pfield->pvx)[nz][i][it+1], nx-it, MPI_DOUBLE, myidr, 61, &(pfield->pvx)[0][i][it+1] ,nx-it, MPI_DOUBLE, myidl, 61, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&(pfield->vth)[nz][i][it+1], nx-it, MPI_DOUBLE, myidr, 62, &(pfield->vth)[0][i][it+1] ,nx-it, MPI_DOUBLE, myidl, 62, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&(pfield->vre)[nz][i][it+1], nx-it, MPI_DOUBLE, myidr, 63, &(pfield->vre)[0][i][it+1] ,nx-it, MPI_DOUBLE, myidl, 63, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&(pfield->p)[nz][i][it+1], nx-it, MPI_DOUBLE, myidr, 64, &(pfield->p)[0][i][it+1] ,nx-it, MPI_DOUBLE, myidl, 64, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&(pfield->q16)[n][nz][i][it+1], nx-it, MPI_DOUBLE, myidr, 65, &(pfield->q16)[n][0][i][it+1] ,nx-it, MPI_DOUBLE, myidl, 65, MPI_COMM_WORLD, &status);

	}

	for(int i=jt+1;i<ny+1;i++)
	{
		MPI_Sendrecv(&(pfield->q11)[n][nz][i][ib], it-ib+1, MPI_DOUBLE, myidr, 70, &(pfield->q11)[n][0][i][ib] ,it-ib+1, MPI_DOUBLE, myidl, 70, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&(pfield->pvx)[nz][i][ib], it-ib+1, MPI_DOUBLE, myidr, 71, &(pfield->pvx)[0][i][ib] ,it-ib+1, MPI_DOUBLE, myidl, 71, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&(pfield->vth)[nz][i][ib], it-ib+1, MPI_DOUBLE, myidr, 72, &(pfield->vth)[0][i][ib] ,it-ib+1, MPI_DOUBLE, myidl, 72, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&(pfield->vre)[nz][i][ib], it-ib+1, MPI_DOUBLE, myidr, 73, &(pfield->vre)[0][i][ib] ,it-ib+1, MPI_DOUBLE, myidl, 73, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&(pfield->p)[nz][i][ib], it-ib+1, MPI_DOUBLE, myidr, 74, &(pfield->p)[0][i][ib] ,it-ib+1, MPI_DOUBLE, myidl, 74, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&(pfield->q16)[n][nz][i][ib], it-ib+1, MPI_DOUBLE, myidr, 75, &(pfield->q16)[n][0][i][ib] ,it-ib+1, MPI_DOUBLE, myidl, 75, MPI_COMM_WORLD, &status);
	}
}

void CBoundaryCondition::update()
{
	double y1,z1,rr,sir,cor,qq2;
	
	for (int j=1;j<ny+1;j++)
        for (int i=1;i<nx+1;i++) 
        {
		if (!((i >= ib) && (i <= it) && (j >= jb) && (j <= jt)))
                {
                    	y1 = (pgrid->yy0)[0][j][i];
                        z1 = (pgrid->zz0)[0][j][i];
                        rr = sqrt(y1*y1 + z1*z1);                                                  
                        sir = z1 / rr;                                                        
                        cor = y1 / rr;                                                                                                               
			(pfield->pvy)[0][j][i] = (pfield->vre)[0][j][i] * cor - (pfield->vth)[0][j][i] * sir;                   
			(pfield->pvz)[0][j][i] = (pfield->vre)[0][j][i] * sir + (pfield->vth)[0][j][i] * cor;
                        y1 = (pgrid->yy0)[nz+1][j][i];
                        z1 = (pgrid->zz0)[nz+1][j][i];
                        rr = sqrt(y1*y1 + z1*z1);
                        sir = z1 / rr;
                        cor = y1 / rr;
                        (pfield->pvy)[nz+1][j][i] = (pfield->vre)[nz+1][j][i] * cor - (pfield->vth)[nz+1][j][i] * sir;
                        (pfield->pvz)[nz+1][j][i] = (pfield->vre)[nz+1][j][i] * sir + (pfield->vth)[nz+1][j][i] * cor;  
                  }            
                        (pfield->q12)[n][0][j][i] = (pfield->q11)[n][0][j][i] * (pfield->pvx)[0][j][i];
                        (pfield->q13)[n][0][j][i] = (pfield->q11)[n][0][j][i] * (pfield->pvy)[0][j][i];
                        (pfield->q14)[n][0][j][i] = (pfield->q11)[n][0][j][i] * (pfield->pvz)[0][j][i];
                        qq2 = (pfield->pvx)[0][j][i] * (pfield->pvx)[0][j][i] + (pfield->pvy)[0][j][i] * (pfield->pvy)[0][j][i] + (pfield->pvz)[0][j][i] * (pfield->pvz)[0][j][i];
                        (pfield->q15)[n][0][j][i] = 2.5*(pfield->p)[0][j][i] + 0.5*(pfield->q11)[n][0][j][i] * qq2;
                        (pfield->t)[0][j][i] = (pfield->p)[0][j][i] / ((pfield->q11)[n][0][j][i] * (pdict->rg));
                        (pfield->q12)[n][nz+1][j][i] = (pfield->q11)[n][nz+1][j][i] * (pfield->pvx)[nz+1][j][i];
                        (pfield->q13)[n][nz+1][j][i] = (pfield->q11)[n][nz+1][j][i] * (pfield->pvy)[nz+1][j][i];
                        (pfield->q14)[n][nz+1][j][i] = (pfield->q11)[n][nz+1][j][i] * (pfield->pvz)[nz+1][j][i];
                        qq2 = (pfield->pvx)[nz+1][j][i] * (pfield->pvx)[nz+1][j][i] + (pfield->pvy)[nz+1][j][i] * (pfield->pvy)[nz+1][j][i] + (pfield->pvz)[nz+1][j][i] * (pfield->pvz)[nz+1][j][i];
                        (pfield->q15)[n][nz+1][j][i] = 2.5*(pfield->p)[nz+1][j][i] + 0.5*(pfield->q11)[n][nz+1][j][i] * qq2;
                        (pfield->t)[nz+1][j][i] = (pfield->p)[nz+1][j][i] / ((pfield->q11)[n][nz+1][j][i] * (pdict->rg));
         }
}

bool CBoundaryCondition::ExchangeBoundaryCondition()
{
	//SetMultigridRange(n1, nng);
	setBC_in();
	setBC_out();
	setBC_root_top( );
	setBC_circumference();
	communication();
	update();
	return true;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



