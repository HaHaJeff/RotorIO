/*---------------------------------------------------------------------------*\


\*---------------------------------------------------------------------------*/


#include "Diffusion.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



void CDiffusion::ResetViscousData()
{
	for (int k = 0; k<nz + 2; k++)
	for (int j = 0; j<ny + 2; j++)
	for (int i = 0; i<nx + 2; i++)
	{
		qv2[k][j][i] = 0;
		qv3[k][j][i] = 0;
		qv4[k][j][i] = 0;
		qv5[k][j][i] = 0;
	}
	//Memset(qv2[0][0], 0, (nz + 2)*(ny + 2)*(nx + 2)*sizeof(double));
}
void CDiffusion::viscosity(double temp, double q6, double &cv, double &kc)
{
	double cvl, cvt, fv1, tem;

	cvl = pdict->cvl0*pow(temp / pdict->t0, 1.5)*(pdict->t0 + pdict->ts) / (temp + pdict->ts);
	tem = q6 / cvl;
	fv1 = 1.0 / (1.0 + pow(pdict->cv1 / tem, 3));
	cvt = q6*fv1;
	cv = cvl + cvt;
	kc = pdict->cp*(cvl / pdict->prl + cvt / pdict->prt);
}
void CDiffusion::ComputeViscousData_i()
{
	double flu2, flu3, flu4, flu5, tauxx, tauyy, tauzz, tauxy, tauxz, tauyz, phix, phiy, phiz;
	double two3, uav, vav, wav, q16av, tav, mav, kav;	
	two3 = 2.0 / 3.0;

	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 2; i++)
	{
		uav = 0.50*(pfield->pvx[k][j][i - 1] + pfield->pvx[k][j][i]);
		vav = 0.50*(pfield->pvy[k][j][i - 1] + pfield->pvy[k][j][i]);
		wav = 0.50*(pfield->pvz[k][j][i - 1] + pfield->pvz[k][j][i]);
		//
		q16av = 0.50*(pfield->q16[n][k][j][i - 1] + pfield->q16[n][k][j][i]);
		tav = 0.50*(pfield->t[k][j][i - 1] + pfield->t[k][j][i]);
		viscosity(tav, q16av, mav, kav);
		tauxx = two3*mav*(2.0*gradfi[1][k][j][i] - gradfi[5][k][j][i] - gradfi[9][k][j][i]);
		tauyy = two3*mav*(2.0*gradfi[5][k][j][i] - gradfi[1][k][j][i] - gradfi[9][k][j][i]);
		tauzz = two3*mav*(2.0*gradfi[9][k][j][i] - gradfi[1][k][j][i] - gradfi[5][k][j][i]);
		tauxy = mav*(gradfi[2][k][j][i] + gradfi[4][k][j][i]);
		tauxz = mav*(gradfi[3][k][j][i] + gradfi[7][k][j][i]);
		tauyz = mav*(gradfi[6][k][j][i] + gradfi[8][k][j][i]);
		phix = uav*tauxx + vav*tauxy + wav*tauxz + kav*gradfi[10][k][j][i];
		phiy = uav*tauxy + vav*tauyy + wav*tauyz + kav*gradfi[11][k][j][i];
		phiz = uav*tauxz + vav*tauyz + wav*tauzz + kav*gradfi[12][k][j][i];
		flu2 = pgrid->s2x[k][j][i] * tauxx + pgrid->s2y[k][j][i] * tauxy + pgrid->s2z[k][j][i] * tauxz;
		flu3 = pgrid->s2x[k][j][i] * tauxy + pgrid->s2y[k][j][i] * tauyy + pgrid->s2z[k][j][i] * tauyz;
		flu4 = pgrid->s2x[k][j][i] * tauxz + pgrid->s2y[k][j][i] * tauyz + pgrid->s2z[k][j][i] * tauzz;
		flu5 = pgrid->s2x[k][j][i] * phix + pgrid->s2y[k][j][i] * phiy + pgrid->s2z[k][j][i] * phiz;
		qv2[k][j][i] = qv2[k][j][i] + flu2;
		qv3[k][j][i] = qv3[k][j][i] + flu3;
		qv4[k][j][i] = qv4[k][j][i] + flu4;
		qv5[k][j][i] = qv5[k][j][i] + flu5;
		qv2[k][j][i - 1] = qv2[k][j][i - 1] - flu2;
		qv3[k][j][i - 1] = qv3[k][j][i - 1] - flu3;
		qv4[k][j][i - 1] = qv4[k][j][i - 1] - flu4;
		qv5[k][j][i - 1] = qv5[k][j][i - 1] - flu5;
	}
}
void CDiffusion::ComputeViscousData_j()
{
	double flu2, flu3, flu4, flu5, tauxx, tauyy, tauzz, tauxy, tauxz, tauyz, phix, phiy, phiz;
	double two3, uav, vav, wav, q16av, tav, mav, kav;	
	two3 = 2.0 / 3.0;
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 2; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		uav = 0.50*(pfield->pvx[k][j - 1][i] + pfield->pvx[k][j][i]);
		vav = 0.50*(pfield->pvy[k][j - 1][i] + pfield->pvy[k][j][i]);
		wav = 0.50*(pfield->pvz[k][j - 1][i] + pfield->pvz[k][j][i]);
		q16av = 0.50*(pfield->q16[n][k][j - 1][i] + pfield->q16[n][k][j][i]);
		tav = 0.50*(pfield->t[k][j - 1][i] + pfield->t[k][j][i]);
		viscosity(tav, q16av, mav, kav);
		tauxx = two3*mav*(2.0*gradfj[1][k][j][i] - gradfj[5][k][j][i] - gradfj[9][k][j][i]);
		tauyy = two3*mav*(2.0*gradfj[5][k][j][i] - gradfj[1][k][j][i] - gradfj[9][k][j][i]);
		tauzz = two3*mav*(2.0*gradfj[9][k][j][i] - gradfj[1][k][j][i] - gradfj[5][k][j][i]);
		tauxy = mav*(gradfj[2][k][j][i] + gradfj[4][k][j][i]);
		tauxz = mav*(gradfj[3][k][j][i] + gradfj[7][k][j][i]);
		tauyz = mav*(gradfj[6][k][j][i] + gradfj[8][k][j][i]);
		phix = uav*tauxx + vav*tauxy + wav*tauxz + kav*gradfj[10][k][j][i];
		phiy = uav*tauxy + vav*tauyy + wav*tauyz + kav*gradfj[11][k][j][i];
		phiz = uav*tauxz + vav*tauyz + wav*tauzz + kav*gradfj[12][k][j][i];
		flu2 = pgrid->s3x[k][j][i] * tauxx + pgrid->s3y[k][j][i] * tauxy + pgrid->s3z[k][j][i] * tauxz;
		flu3 = pgrid->s3x[k][j][i] * tauxy + pgrid->s3y[k][j][i] * tauyy + pgrid->s3z[k][j][i] * tauyz;
		flu4 = pgrid->s3x[k][j][i] * tauxz + pgrid->s3y[k][j][i] * tauyz + pgrid->s3z[k][j][i] * tauzz;
		flu5 = pgrid->s3x[k][j][i] * phix + pgrid->s3y[k][j][i] * phiy + pgrid->s3z[k][j][i] * phiz;
		qv2[k][j][i] = qv2[k][j][i] + flu2;
		qv3[k][j][i] = qv3[k][j][i] + flu3;
		qv4[k][j][i] = qv4[k][j][i] + flu4;
		qv5[k][j][i] = qv5[k][j][i] + flu5;
		qv2[k][j - 1][i] = qv2[k][j - 1][i] - flu2;
		qv3[k][j - 1][i] = qv3[k][j - 1][i] - flu3;
		qv4[k][j - 1][i] = qv4[k][j - 1][i] - flu4;
		qv5[k][j - 1][i] = qv5[k][j - 1][i] - flu5;
	}
}
void CDiffusion::ComputeViscousData_k()
{
	double flu2, flu3, flu4, flu5, tauxx, tauyy, tauzz, tauxy, tauxz, tauyz, phix, phiy, phiz;
	double two3, uav, vav, wav, q16av, tav, mav, kav;	
	two3 = 2.0 / 3.0;
	for (int k = 1; k<nz + 2; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		uav = 0.50*(pfield->pvx[k - 1][j][i] + pfield->pvx[k][j][i]);
		vav = 0.50*(pfield->pvy[k - 1][j][i] + pfield->pvy[k][j][i]);
		wav = 0.50*(pfield->pvz[k - 1][j][i] + pfield->pvz[k][j][i]);
		q16av = 0.50*(pfield->q16[n][k - 1][j][i] + pfield->q16[n][k][j][i]);
		tav = 0.50*(pfield->t[k - 1][j][i] + pfield->t[k][j][i]);
		viscosity(tav, q16av, mav, kav);
		tauxx = two3*mav*(2.0*gradfk[1][k][j][i] - gradfk[5][k][j][i] - gradfk[9][k][j][i]);
		tauyy = two3*mav*(2.0*gradfk[5][k][j][i] - gradfk[1][k][j][i] - gradfk[9][k][j][i]);
		tauzz = two3*mav*(2.0*gradfk[9][k][j][i] - gradfk[1][k][j][i] - gradfk[5][k][j][i]);
		tauxy = mav*(gradfk[2][k][j][i] + gradfk[4][k][j][i]);
		tauxz = mav*(gradfk[3][k][j][i] + gradfk[7][k][j][i]);
		tauyz = mav*(gradfk[6][k][j][i] + gradfk[8][k][j][i]);
		phix = uav*tauxx + vav*tauxy + wav*tauxz + kav*gradfk[10][k][j][i];
		phiy = uav*tauxy + vav*tauyy + wav*tauyz + kav*gradfk[11][k][j][i];
		phiz = uav*tauxz + vav*tauyz + wav*tauzz + kav*gradfk[12][k][j][i];
		flu2 = pgrid->s1x[k][j][i] * tauxx + pgrid->s1y[k][j][i] * tauxy + pgrid->s1z[k][j][i] * tauxz;
		flu3 = pgrid->s1x[k][j][i] * tauxy + pgrid->s1y[k][j][i] * tauyy + pgrid->s1z[k][j][i] * tauyz;
		flu4 = pgrid->s1x[k][j][i] * tauxz + pgrid->s1y[k][j][i] * tauyz + pgrid->s1z[k][j][i] * tauzz;
		flu5 = pgrid->s1x[k][j][i] * phix + pgrid->s1y[k][j][i] * phiy + pgrid->s1z[k][j][i] * phiz;
		qv2[k][j][i] = qv2[k][j][i] + flu2;
		qv3[k][j][i] = qv3[k][j][i] + flu3;
		qv4[k][j][i] = qv4[k][j][i] + flu4;
		qv5[k][j][i] = qv5[k][j][i] + flu5;
		qv2[k - 1][j][i] = qv2[k - 1][j][i] - flu2;
		qv3[k - 1][j][i] = qv3[k - 1][j][i] - flu3;
		qv4[k - 1][j][i] = qv4[k - 1][j][i] - flu4;
		qv5[k - 1][j][i] = qv5[k - 1][j][i] - flu5;
	}
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		// 
		qv3[k][j][i] = qv3[k][j][i] + pdict->rpm*pfield->q14[n][k][j][i] * pgrid->vv[k][j][i];
		qv4[k][j][i] = qv4[k][j][i] - pdict->rpm*pfield->q13[n][k][j][i] * pgrid->vv[k][j][i];
	}
}
void CDiffusion::ComputeViscousData()
{
	gradsface();
	ResetViscousData();
	ComputeViscousData_i();
	ComputeViscousData_j();
	ComputeViscousData_k();
}
void CDiffusion::SAsource()
{
	double*** tur,***pwx,***pwy,***pwz;

	Malloc(tur,nz + 2,ny + 2,nx + 2);
	Malloc(pwx,nz + 2,ny + 2,nx + 2);
	Malloc(pwy,nz + 2,ny + 2,nx + 2);
	Malloc(pwz,nz + 2,ny + 2,nx + 2);


	double cvl, tem, fv1, fv2, fv3, rp1, rp2, rp3, w12, w13, w23, w21, w31, w32, ww2, ww, svot, vm, gv;
	double s11, s22, s33, s12, s13, s23, ss2, ss, dd, fr1r1, fr1r2, fr1, yv, ra, ga, fw, dv;
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	{
		pgrid->yy0[k][j][0] = pgrid->yy0[k][j][1];
		pgrid->zz0[k][j][0] = pgrid->zz0[k][j][1];
		pgrid->yy0[k][j][nx + 1] = pgrid->yy0[k][j][nx];
		pgrid->zz0[k][j][nx + 1] = pgrid->zz0[k][j][nx];
	}
	for (int k = 1; k<nz + 1; k++)
	for (int i = 1; i<nx + 1; i++)
	{
		pgrid->yy0[k][0][i] = pgrid->yy0[k][1][i];
		pgrid->zz0[k][0][i] = pgrid->zz0[k][1][i];
		pgrid->yy0[k][ny + 1][i] = pgrid->yy0[k][ny][i];
		pgrid->zz0[k][ny + 1][i] = pgrid->zz0[k][ny][i];
	}
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		pgrid->yy0[0][j][i] = pgrid->yy0[1][j][i];
		pgrid->zz0[0][j][i] = pgrid->zz0[1][j][i];
		pgrid->yy0[nz + 1][j][i] = pgrid->yy0[nz][j][i];
		pgrid->zz0[nz + 1][j][i] = pgrid->zz0[nz][j][i];
	}
	for (int k = 0; k<nz + 2; k++)
	for (int j = 0; j<ny + 2; j++)
	for (int i = 0; i<nx + 2; i++)
	{
		tur[k][j][i] = pfield->q16[n][k][j][i] / pfield->q11[n][k][j][i];
		pwx[k][j][i] = pfield->pvx[k][j][i];
		pwy[k][j][i] = pfield->pvy[k][j][i] + pdict->rpm*pgrid->zz0[k][j][i];
		pwz[k][j][i] = pfield->pvz[k][j][i] - pdict->rpm*pgrid->yy0[k][j][i];
	}


	//******************
	for (int l = 1; l<13; l++)
	for (int k = 0; k<nz + 2; k++)
	for (int j = 0; j<ny + 2; j++)
	for (int i = 0; i<nx + 2; i++)
	{
		gradc[l][k][j][i] = 0;
	}
	gradscentre(1, pwx, gradc, 1);
	gradscentre(2, pwx, gradc, 2);
	gradscentre(3, pwx, gradc, 3);
	gradscentre(1, pwy, gradc, 4);
	gradscentre(2, pwy, gradc, 5);
	gradscentre(3, pwy, gradc, 6);
	gradscentre(1, pwz, gradc, 7);
	gradscentre(2, pwz, gradc, 8);
	gradscentre(3, pwz, gradc, 9);
	gradscentre(1, tur, gradc, 10);
	gradscentre(2, tur, gradc, 11);
	gradscentre(3, tur, gradc, 12);

	for (int l = 1; l<10; l++)
	for (int k = 0; k<nz + 2; k++)
	for (int j = 0; j<ny + 2; j++)
	for (int i = 0; i<nx + 2; i++)
	{
		gradcs[l][k][j][i] = 0;
	}

	dsdt(gradc, gradcs, 1,pdict->rpm);
	dsdt(gradc, gradcs, 2, pdict->rpm);
	dsdt(gradc, gradcs, 3, pdict->rpm);
	dsdt(gradc, gradcs, 4, pdict->rpm);
	dsdt(gradc, gradcs, 5, pdict->rpm);
	dsdt(gradc, gradcs, 6, pdict->rpm);
	dsdt(gradc, gradcs, 7, pdict->rpm);
	dsdt(gradc, gradcs, 8, pdict->rpm);
	dsdt(gradc, gradcs, 9, pdict->rpm);


	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		dv = pdict->cb2 / pdict->sigmav*(pow(gradc[10][k][j][i], 2) + pow(gradc[11][k][j][i], 2) + pow(gradc[12][k][j][i], 2))*pfield->q11[n][k][j][i];
		cvl = pdict->cvl0*pow(pfield->t[k][j][i] / pdict->t0, 1.5)*(pdict->t0 + pdict->ts) / (pfield->t[k][j][i] + pdict->ts);
		tem = pfield->q16[n][k][j][i] / cvl;
		tem = max(tem, pow(10.0, -4));
		fv1 = 1.0 / (1.0 + pow(pdict->cv1 / tem, 3));
		fv2 = pow(1.0 + tem / pdict->cv2, -3);
		fv3 = (1.0 + tem*fv1)*(1.0 - fv2) / tem;
		w12 = 0.50*(gradc[2][k][j][i] - gradc[4][k][j][i]);
		w13 = 0.50*(gradc[3][k][j][i] - gradc[7][k][j][i]);
		w23 = 0.50*(gradc[6][k][j][i] - gradc[8][k][j][i]);
		ww2 = 4.0*(w12*w12 + w13*w13 + w23*w23);
		ww = sqrt(ww2);
		vm = tur[k][j][i] / (pdict->kap*pdict->kap*(((pfield->mesh)->dmini))[k][j][i] * (((pfield->mesh)->dmini))[k][j][i]);
		svot = ww*fv3 + vm*fv2;
		gv = pdict->cb1*pfield->q16[n][k][j][i] * svot;
		s11 = gradc[1][k][j][i];
		s22 = gradc[5][k][j][i];
		s33 = gradc[9][k][j][i];
		s12 = 0.50*(gradc[2][k][j][i] + gradc[4][k][j][i]);
		s13 = 0.50*(gradc[3][k][j][i] + gradc[7][k][j][i]);
		s23 = 0.50*(gradc[6][k][j][i] + gradc[8][k][j][i]);
		ss2 = 4.0*(s12*s12 + s13*s13 + s23*s23) + 2.0*(s11*s11 + s22*s22 + s33*s33);
		ss = sqrt(ss2);
		rp1 = pdict->rpm;
		rp2 = 0.0;
		rp3 = 0.0;
		w12 = 0.50*(gradc[2][k][j][i] - gradc[4][k][j][i]) - rp3;
		w13 = 0.50*(gradc[3][k][j][i] - gradc[7][k][j][i]) + rp2;
		w23 = 0.50*(gradc[6][k][j][i] - gradc[8][k][j][i]) - rp1;
		w21 = -w12;
		w31 = -w13;
		w32 = -w23;
		ww2 = 4.0*(w12*w12 + w13*w13 + w23*w23);
		ww = sqrt(ww2);
		fr1r1 = ss / ww;
		dd = max(ss2, 0.09*tur[k][j][i] * tur[k][j][i]);
		fr1r2 = 2.0 / dd / sqrt(dd) / ww*(gradcs[1][k][j][i] * (w12*s12 + w13*s13)
			+ (0.50*(gradcs[2][k][j][i] + gradcs[4][k][j][i]) - s13*rp1)*(w12*s22 + w13*s23 + w21*s11 + w23*s13)
			+ (0.50*(gradcs[3][k][j][i] + gradcs[7][k][j][i]) + s12*rp1)*(w12*s23 + w13*s33 + w31*s11 + w32*s12) + (gradcs[5][k][j][i] - 2.0*s23*rp1)*(w21*s12 + w23*s23)
			+ (0.50*(gradcs[6][k][j][i] + gradcs[8][k][j][i]) + (s22 - s33)*rp1)*(w21*s13 + w23*s33 + w31*s12 + w32*s22) + (gradcs[9][k][j][i] + 2.0*s23*rp1)*(w31*s13 + w32*s23));
		fr1 = (1.0 + pdict->cr1)*2.0*fr1r1 / (1.0 + fr1r1)*(1.0 - pdict->cr3*atan(pdict->cr2*fr1r2)) - pdict->cr1;
		gv = gv*fr1;
		ra = vm / svot;
		ga = ra + pdict->cw2*(pow(ra, 6) - ra);
		fw = pow((pow(ga, -6) + pow(pdict->cw3, -6)) / (1.0 + pow(pdict->cw3, -6)), (-1.0 / 6.0));
		yv = pdict->cw1*fw*pfield->q16[n][k][j][i] * vm*pdict->kap*pdict->kap;
		qv6[k][j][i] = qv6[k][j][i] + (gv - yv + dv)*pgrid->vv[k][j][i];

	}

	Free(tur);
	Free(pwx);
	Free(pwy);
	Free(pwz);
}
void CDiffusion::ComputeViscousData_SA()
{
	double flu6, q16av, tav, cvl;

	double*** tur;
	Malloc(tur,nz + 2,ny + 2,nx + 2);


	for (int i = 0; i<nz + 2; i++)
	for (int j = 0; j<ny + 2; j++)
	for (int k = 0; k<nx + 2; k++)
	{
		tur[i][j][k] = pfield->q16[n][i][j][k] / pfield->q11[n][i][j][k];
	}

	ResetGradient();	
	gradsfaceI(pgrid->s2x, pgrid->s3x, pgrid->s1x, tur, gradfi, 13);
	gradsfaceI(pgrid->s2y, pgrid->s3y, pgrid->s1y, tur, gradfi, 14);
	gradsfaceI(pgrid->s2z, pgrid->s3z, pgrid->s1z, tur, gradfi, 15);

	gradsfaceJ(pgrid->s2x, pgrid->s3x, pgrid->s1x, tur, gradfj, 13);
	gradsfaceJ(pgrid->s2y, pgrid->s3y, pgrid->s1y, tur, gradfj, 14);
	gradsfaceJ(pgrid->s2z, pgrid->s3z, pgrid->s1z, tur, gradfj, 15);

	gradsfaceK(pgrid->s2x, pgrid->s3x, pgrid->s1x, tur, gradfk, 13);
	gradsfaceK(pgrid->s2y, pgrid->s3y, pgrid->s1y, tur, gradfk, 14);
	gradsfaceK(pgrid->s2z, pgrid->s3z, pgrid->s1z, tur, gradfk, 15);

	for (int i = 0; i<nz + 2; i++)
	for (int j = 0; j<ny + 2; j++)
	for (int k = 0; k<nx + 2; k++)
	{
		qv6[i][j][k] = 0;
	}

	for (int i = 1; i<nz + 1; i++)
	for (int j = 1; j<ny + 1; j++)
	for (int k = 1; k<nx + 2; k++)
	{
		tav = 0.50*(pfield->t[i][j][k - 1] + pfield->t[i][j][k]);
		cvl = pdict->cvl0*pow((tav / pdict->t0), 1.5)*(pdict->t0 + pdict->ts) / (tav + pdict->ts);
		q16av = 0.50*(pfield->q16[n][i][j][k - 1] + pfield->q16[n][i][j][k]);
		flu6 = (cvl + q16av)*(gradfi[13][i][j][k] * pgrid->s2x[i][j][k] + gradfi[14][i][j][k] * pgrid->s2y[i][j][k] + gradfi[15][i][j][k] * pgrid->s2z[i][j][k]) / pdict->sigmav;
		qv6[i][j][k] = qv6[i][j][k] + flu6;
		qv6[i][j][k - 1] = qv6[i][j][k - 1] - flu6;
	}

	for (int i = 1; i<nz + 1; i++)
	for (int j = 1; j<ny + 2; j++)
	for (int k = 1; k<nx + 1; k++)
	{
		tav = 0.50*(pfield->t[i][j - 1][k] + pfield->t[i][j][k]);
		cvl = pdict->cvl0*pow((tav / pdict->t0), 1.5)*(pdict->t0 + pdict->ts) / (tav + pdict->ts);
		q16av = 0.50*(pfield->q16[n][i][j - 1][k] + pfield->q16[n][i][j][k]);
		flu6 = (cvl + q16av)*(gradfj[13][i][j][k] * pgrid->s3x[i][j][k] + gradfj[14][i][j][k] * pgrid->s3y[i][j][k] + gradfj[15][i][j][k] * pgrid->s3z[i][j][k]) / pdict->sigmav;
		qv6[i][j][k] = qv6[i][j][k] + flu6;
		qv6[i][j - 1][k] = qv6[i][j - 1][k] - flu6;
	}

	for (int i = 1; i<nz + 2; i++)
	for (int j = 1; j<ny + 1; j++)
	for (int k = 1; k<nx + 1; k++)
	{
		tav = 0.50*(pfield->t[i - 1][j][k] + pfield->t[i][j][k]);
		cvl = pdict->cvl0*pow((tav / pdict->t0), 1.5)*(pdict->t0 + pdict->ts) / (tav + pdict->ts);
		q16av = 0.50*(pfield->q16[n][i - 1][j][k] + pfield->q16[n][i][j][k]);
		flu6 = (cvl + q16av)*(gradfk[13][i][j][k] * pgrid->s1x[i][j][k] + gradfk[14][i][j][k] * pgrid->s1y[i][j][k] + gradfk[15][i][j][k] * pgrid->s1z[i][j][k]) / pdict->sigmav;
		qv6[i][j][k] = qv6[i][j][k] + flu6;
		qv6[i - 1][j][k] = qv6[i - 1][j][k] - flu6;
	}
	Free(tur);

	SAsource();
}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
