/*---------------------------------------------------------------------------*\


\*---------------------------------------------------------------------------*/

#include "RungeKutta.h"
//***************************************************************************


//***************************************************************************
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
void CRungeKutta::ResetP2h()
{
	for (int i = 1; i<nz + 1; i++)   //						
		for (int j = 1; j<ny + 1; j++)
			for (int k = 1; k<nx + 1; k++)
			{
				qp1[i][j][k] = 0.0;
				qp2[i][j][k] = 0.0;
				qp3[i][j][k] = 0.0;
				qp4[i][j][k] = 0.0;
				qp5[i][j][k] = 0.0;
			}
}
void CRungeKutta::InitPressure()
{
	double qq2, cvl, sir, cor;
	double vx, vy, vz, y1, z1, rr, dim;
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
			{
				dim = pfield->q11[n][k][j][i];
				pfield->pvx[k][j][i] = pfield->q12[n][k][j][i] / dim;
				pfield->pvy[k][j][i] = pfield->q13[n][k][j][i] / dim;
				pfield->pvz[k][j][i] = pfield->q14[n][k][j][i] / dim;
				vx = pfield->pvx[k][j][i];
				vy = pfield->pvy[k][j][i];
				vz = pfield->pvz[k][j][i];
				y1 = pgrid->yy0[k][j][i];
				z1 = pgrid->zz0[k][j][i];
				rr = sqrt(y1*y1 + z1*z1);
				sir = z1 / rr;
				cor = y1 / rr;
				pfield->vth[k][j][i] = vz*cor - vy*sir;    //\D6\DC\CF\F2\CBٶ\C8
				pfield->vre[k][j][i] = vz*sir + vy*cor;    //\BE\B6\CF\F2\CBٶ\C8
				qq2 = vx*vx + vy*vy + vz*vz;
				pfield->p[k][j][i] = 0.4*(pfield->q15[n][k][j][i] - 0.5*dim*qq2);
				pfield->t[k][j][i] = pfield->p[k][j][i] / (dim*pdict->rg);
				cvl = pdict->cvl0*pow(pfield->t[k][j][i] / pdict->t0, 1.5)*(pdict->t0 + pdict->ts) / (pfield->t[k][j][i] + pdict->ts);
				pfield->q16[n][k][j][i] = max(pfield->q16[n][k][j][i], pow(10.0, -4)*cvl);
			}
}
void CRungeKutta::InitPreField()
{
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
			{
				q01[k][j][i] = pfield->q11[n][k][j][i];
				q02[k][j][i] = pfield->q12[n][k][j][i];
				q03[k][j][i] = pfield->q13[n][k][j][i];
				q04[k][j][i] = pfield->q14[n][k][j][i];
				q05[k][j][i] = pfield->q15[n][k][j][i];
			}
}
void CRungeKutta::InitPreField_SA()
{
	for (int i = 1; i<nz + 1; i++)
	for (int j = 1; j<ny + 1; j++)
	for (int k = 1; k<nx + 1; k++)
	{
		q06[i][j][k] = pfield->q16[n][i][j][k];
	}
}
void CRungeKutta::RK_I()
{
	InitPressure();
	ExchangeBoundaryCondition();
	step_artificial();
	ComputeConvectiveData();
	ComputeViscousData();
}
void CRungeKutta::RK_C_I()
{
	InitPressure();
	ExchangeBoundaryCondition();
	step_artificial_c();
	ComputeConvectiveData();
	ComputeViscousData();
}
void CRungeKutta::RK_II()
{
	InitPressure();
	ExchangeBoundaryCondition();
	ComputeConvectiveData();
}
void CRungeKutta::RK_III()
{
	InitPressure();
	ExchangeBoundaryCondition();
	ComputeConvectiveData();
}
void CRungeKutta::RungeKutta_I(int nng)
{
	SetMultigridRange(nng);
	InitPreField();
	ta = 1.0;
	timl = 0.6*pdict->cfl;
	RK_I( );
	UpdateFieldValue(ta, timl, 1);
	/*for (int k = 1; k<nz + 1; k++)
		for (int j = 1; j<ny + 1; j++)
			for (int i = 1; i<nx + 1; i++)
			{
				f1 << setw(15) << setprecision(6) << py1[k][j][i] << endl;
				f1 << setw(15) << setprecision(6) << py5[k][j][i] << endl;
			}*/
	



	timl = 0.6*pdict->cfl;
	RK_II();
	UpdateFieldValue(ta, timl, 0);
	




	timl = pdict->cfl;
	RK_III();
	UpdateFieldValue(ta, timl, 1);
}
void CRungeKutta::RungeKutta_II( int nng)
{
	SetMultigridRange(nng);
	InitPreField();
	ta = 1.0;
	timl = 0.6*pdict->cfl;
	RK_C_I();
	UpdateFieldValue(ta, timl, 1);

	timl = 0.6*pdict->cfl;
	RK_II();
	UpdateFieldValue(ta, timl, 0);

	timl = pdict->cfl;
	RK_III();
	UpdateFieldValue(ta, timl, 1);
}
void CRungeKutta::RungeKutta_SA()
{
	timedrivationCmputeSA();
	InitPreField_SA();
	ta = 1.5;
	timl = 0.6*pdict->cfl;
	InitPressure();
	ExchangeBoundaryCondition();
	step_artificial_SA();
	ComputeConvectiveData_SA();
	ComputeViscousData_SA();
	UpdateFieldValue_SA( ta,timl, 1);

	timl = 0.6*pdict->cfl;
	InitPressure();
	ExchangeBoundaryCondition();
	ComputeConvectiveData_SA();
	UpdateFieldValue_SA( ta,timl, 0);

	timl = pdict->cfl;
	InitPressure();
	ExchangeBoundaryCondition();
	ComputeConvectiveData_SA();
	UpdateFieldValue_SA( ta,timl, 1);

	InitPreField_SA();
	timl = 0.125*pdict->cfl;
	InitPressure();
	ExchangeBoundaryCondition();
	step_artificial_SA();
	ComputeConvectiveData_SA();
	ComputeViscousData_SA();
	UpdateFieldValue_SA(ta,timl,  0);
}

void CRungeKutta::residual(double &rmsm)
{
	rms[n] = 0;
	for (int i = 1; i<nz + 1; i++)
	for (int j = 1; j<ny + 1; j++)
	for (int k = 1; k<nx + 1; k++)
	{
		rms[n] = rms[n] + pow((pfield->q11[n][i][j][k] - q01[i][j][k])*time[i][j][k] / pgrid->vv[i][j][k] / pdict->cfl, 2);
	}
	rms[n] = 0.5*log10(rms[n] / double(nx*ny*nz));
	if (rms[n]>rmsm)
		rmsm = rms[n];
}
void CRungeKutta::RestrictionOp(int ign)
{
	int i1, j1, k1;
	double val1, val2, val3, val4, val5, val6, val7, val8, vall;

	for (int i = 1; i<nz + 1; i++)
	{
		i1 = 2 * i;
		for (int j = 1; j<ny + 1; j++)
		{
			j1 = 2 * j;
			for (int k = 1; k<nx + 1; k++)
			{
				k1 = 2 * k;
				val1 = pgrid->vvn[ign + 1][i1][j1][k1];
				val2 = pgrid->vvn[ign + 1][i1][j1][k1 - 1];
				val3 = pgrid->vvn[ign + 1][i1 - 1][j1][k1];
				val4 = pgrid->vvn[ign + 1][i1 - 1][j1][k1 - 1];
				val5 = pgrid->vvn[ign + 1][i1][j1 - 1][k1];
				val6 = pgrid->vvn[ign + 1][i1][j1 - 1][k1 - 1];
				val7 = pgrid->vvn[ign + 1][i1 - 1][j1 - 1][k1];
				val8 = pgrid->vvn[ign + 1][i1 - 1][j1 - 1][k1 - 1];
				vall = val1 + val2 + val3 + val4 + val5 + val6 + val7 + val8;
				for (int m = 1; m<pdict->nt + 1; m++)
				{
					pfield->q11[m][i][j][k] = (pfield->q11[m][i1][j1][k1] * val1 + pfield->q11[m][i1][j1][k1 - 1] * val2 + pfield->q11[m][i1 - 1][j1][k1] * val3 + pfield->q11[m][i1 - 1][j1][k1 - 1] * val4 +
						pfield->q11[m][i1][j1 - 1][k1] * val5 + pfield->q11[m][i1][j1 - 1][k1 - 1] * val6 + pfield->q11[m][i1 - 1][j1 - 1][k1] * val7 + pfield->q11[m][i1 - 1][j1 - 1][k1 - 1] * val8) / vall;

					pfield->q12[m][i][j][k] = (pfield->q12[m][i1][j1][k1] * val1 + pfield->q12[m][i1][j1][k1 - 1] * val2 + pfield->q12[m][i1 - 1][j1][k1] * val3 + pfield->q12[m][i1 - 1][j1][k1 - 1] * val4 +
						pfield->q12[m][i1][j1 - 1][k1] * val5 + pfield->q12[m][i1][j1 - 1][k1 - 1] * val6 + pfield->q12[m][i1 - 1][j1 - 1][k1] * val7 + pfield->q12[m][i1 - 1][j1 - 1][k1 - 1] * val8) / vall;

					pfield->q13[m][i][j][k] = (pfield->q13[m][i1][j1][k1] * val1 + pfield->q13[m][i1][j1][k1 - 1] * val2 + pfield->q13[m][i1 - 1][j1][k1] * val3 + pfield->q13[m][i1 - 1][j1][k1 - 1] * val4 +
						pfield->q13[m][i1][j1 - 1][k1] * val5 + pfield->q13[m][i1][j1 - 1][k1 - 1] * val6 + pfield->q13[m][i1 - 1][j1 - 1][k1] * val7 + pfield->q13[m][i1 - 1][j1 - 1][k1 - 1] * val8) / vall;

					pfield->q14[m][i][j][k] = (pfield->q14[m][i1][j1][k1] * val1 + pfield->q14[m][i1][j1][k1 - 1] * val2 + pfield->q14[m][i1 - 1][j1][k1] * val3 + pfield->q14[m][i1 - 1][j1][k1 - 1] * val4 +
						pfield->q14[m][i1][j1 - 1][k1] * val5 + pfield->q14[m][i1][j1 - 1][k1 - 1] * val6 + pfield->q14[m][i1 - 1][j1 - 1][k1] * val7 + pfield->q14[m][i1 - 1][j1 - 1][k1 - 1] * val8) / vall;

					pfield->q15[m][i][j][k] = (pfield->q15[m][i1][j1][k1] * val1 + pfield->q15[m][i1][j1][k1 - 1] * val2 + pfield->q15[m][i1 - 1][j1][k1] * val3 + pfield->q15[m][i1 - 1][j1][k1 - 1] * val4 +
						pfield->q15[m][i1][j1 - 1][k1] * val5 + pfield->q15[m][i1][j1 - 1][k1 - 1] * val6 + pfield->q15[m][i1 - 1][j1 - 1][k1] * val7 + pfield->q15[m][i1 - 1][j1 - 1][k1 - 1] * val8) / vall;

					pfield->q16[m][i][j][k] = (pfield->q16[m][i1][j1][k1] * val1 + pfield->q16[m][i1][j1][k1 - 1] * val2 + pfield->q16[m][i1 - 1][j1][k1] * val3 + pfield->q16[m][i1 - 1][j1][k1 - 1] * val4 +
						pfield->q16[m][i1][j1 - 1][k1] * val5 + pfield->q16[m][i1][j1 - 1][k1 - 1] * val6 + pfield->q16[m][i1 - 1][j1 - 1][k1] * val7 + pfield->q16[m][i1 - 1][j1 - 1][k1 - 1] * val8) / vall;
				}

				rr1[i][j][k] = rr1[i1][j1][k1] + rr1[i1][j1][k1 - 1] + rr1[i1 - 1][j1][k1] + rr1[i1 - 1][j1][k1 - 1] +
					rr1[i1][j1 - 1][k1] + rr1[i1 - 1][j1 - 1][k1] + rr1[i1][j1 - 1][k1 - 1] + rr1[i1 - 1][j1 - 1][k1 - 1];

				rr2[i][j][k] = rr2[i1][j1][k1] + rr2[i1 - 1][j1][k1] + rr2[i1][j1][k1 - 1] + rr2[i1 - 1][j1][k1 - 1] +
					rr2[i1][j1 - 1][k1] + rr2[i1 - 1][j1 - 1][k1] + rr2[i1][j1 - 1][k1 - 1] + rr2[i1 - 1][j1 - 1][k1 - 1];

				rr3[i][j][k] = rr3[i1][j1][k1] + rr3[i1 - 1][j1][k1] + rr3[i1][j1][k1 - 1] + rr3[i1 - 1][j1][k1 - 1] +
					rr3[i1][j1 - 1][k1] + rr3[i1 - 1][j1 - 1][k1] + rr3[i1][j1 - 1][k1 - 1] + rr3[i1 - 1][j1 - 1][k1 - 1];

				rr4[i][j][k] = rr4[i1][j1][k1] + rr4[i1 - 1][j1][k1] + rr4[i1][j1][k1 - 1] + rr4[i1 - 1][j1][k1 - 1] +
					rr4[i1][j1 - 1][k1] + rr4[i1 - 1][j1 - 1][k1] + rr4[i1][j1 - 1][k1 - 1] + rr4[i1 - 1][j1 - 1][k1 - 1];

				rr5[i][j][k] = rr5[i1][j1][k1] + rr5[i1 - 1][j1][k1] + rr5[i1][j1][k1 - 1] + rr5[i1 - 1][j1][k1 - 1] +
					rr5[i1][j1 - 1][k1] + rr5[i1 - 1][j1 - 1][k1] + rr5[i1][j1 - 1][k1 - 1] + rr5[i1 - 1][j1 - 1][k1 - 1];
			}

		}

	}
}
void CRungeKutta::P2hSetting()
{
	timedrivationCmpute();
	RK_C_I();
	for (int i = 1; i<nz + 1; i++)
		for (int j = 1; j<ny + 1; j++)
			for (int k = 1; k<nx + 1; k++)
			{
				qp1[i][j][k] = rr1[i][j][k] + qc1[i][j][k] - av1[i][j][k] + ts1[i][j][k] * pgrid->vv[i][j][k];
				qp2[i][j][k] = rr2[i][j][k] + qc2[i][j][k] - av2[i][j][k] + ts2[i][j][k] * pgrid->vv[i][j][k] - qv2[i][j][k];
				qp3[i][j][k] = rr3[i][j][k] + qc3[i][j][k] - av3[i][j][k] + ts3[i][j][k] * pgrid->vv[i][j][k] - qv3[i][j][k];
				qp4[i][j][k] = rr4[i][j][k] + qc4[i][j][k] - av4[i][j][k] + ts4[i][j][k] * pgrid->vv[i][j][k] - qv4[i][j][k];
				qp5[i][j][k] = rr5[i][j][k] + qc5[i][j][k] - av5[i][j][k] + ts5[i][j][k] * pgrid->vv[i][j][k] - qv5[i][j][k];
			}
}

void CRungeKutta::UpwardCorrection(int ign,int nng)
{
	int ig, i1, j1, k1, iw, ie, jn, js, kf, kb;
	double dpy, dpiw, dpie, dpjn, dpjs, dpkf, dpkb;
	int nnx, nny, nnz;
	for (ig = ign; ig<nng; ig++)
	{
		nnx = (pfield->mesh->nnx)[ig];
		nny = (pfield->mesh->nny)[ig];
		nnz = (pfield->mesh->nnz)[ig];
		for (int k = 1; k<nnz + 1; k++)
		{
			/*		i1 = 2 * i - 1;
			iw = max(1, i - 1);
			ie = min(nx, i + 1);*/
			k1 = 2 * k - 1;
			kf = max(1, k - 1);
			kb = min(nnz, k + 1);

			for (int j = 1; j<nny + 1; j++)
			{
				j1 = 2 * j - 1;
				js = max(1, j - 1);
				jn = min(nny, j + 1);
				for (int i = 1; i<nnx + 1; i++)
				{
					/*		k1 = 2 * k - 1;
					kf = max(1, k - 1);
					kb = min(nz, k + 1);*/
					i1 = 2 * i - 1;
					iw = max(1, i - 1);
					ie = min(nnx, i + 1);

					dpy = py1[k][j][i];
					dpiw = dp(py1[k][j][iw], dpy, py1[k][j][ie]);
					dpie = dp(py1[k][j][ie], dpy, py1[k][j][iw]);
					dpjs = dp(py1[k][js][i], dpy, py1[k][jn][i]);
					dpjn = dp(py1[k][jn][i], dpy, py1[k][js][i]);
					dpkf = dp(py1[kf][j][i], dpy, py1[kb][j][i]);
					dpkb = dp(py1[kb][j][i], dpy, py1[kf][j][i]);
					q01[k1][j1][i1] = dpy + dpiw + dpjs + dpkf;
					q01[k1][j1][i1 + 1] = dpy + dpie + dpjs + dpkf;
					q01[k1 + 1][j1][i1] = dpy + dpiw + dpjs + dpkb;
					q01[k1 + 1][j1][i1 + 1] = dpy + dpie + dpjs + dpkb;
					q01[k1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkf;
					q01[k1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkf;
					q01[k1 + 1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkb;
					q01[k1 + 1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkb;

					dpy = py2[k][j][i];
					dpiw = dp(py2[k][j][iw], dpy, py2[k][j][ie]);
					dpie = dp(py2[k][j][ie], dpy, py2[k][j][iw]);
					dpjs = dp(py2[k][js][i], dpy, py2[k][jn][i]);
					dpjn = dp(py2[k][jn][i], dpy, py2[k][js][i]);
					dpkf = dp(py2[kf][j][i], dpy, py2[kb][j][i]);
					dpkb = dp(py2[kb][j][i], dpy, py2[kf][j][i]);
					q02[k1][j1][i1] = dpy + dpiw + dpjs + dpkf;
					q02[k1][j1][i1 + 1] = dpy + dpie + dpjs + dpkf;
					q02[k1 + 1][j1][i1] = dpy + dpiw + dpjs + dpkb;
					q02[k1 + 1][j1][i1 + 1] = dpy + dpie + dpjs + dpkb;
					q02[k1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkf;
					q02[k1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkf;
					q02[k1 + 1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkb;
					q02[k1 + 1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkb;

					dpy = py3[k][j][i];
					dpiw = dp(py3[k][j][iw], dpy, py3[k][j][ie]);
					dpie = dp(py3[k][j][ie], dpy, py3[k][j][iw]);
					dpjs = dp(py3[k][js][i], dpy, py3[k][jn][i]);
					dpjn = dp(py3[k][jn][i], dpy, py3[k][js][i]);
					dpkf = dp(py3[kf][j][i], dpy, py3[kb][j][i]);
					dpkb = dp(py3[kb][j][i], dpy, py3[kf][j][i]);
					q03[k1][j1][i1] = dpy + dpiw + dpjs + dpkf;
					q03[k1][j1][i1 + 1] = dpy + dpie + dpjs + dpkf;
					q03[k1 + 1][j1][i1] = dpy + dpiw + dpjs + dpkb;
					q03[k1 + 1][j1][i1 + 1] = dpy + dpie + dpjs + dpkb;
					q03[k1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkf;
					q03[k1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkf;
					q03[k1 + 1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkb;
					q03[k1 + 1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkb;

					dpy = py4[k][j][i];
					dpiw = dp(py4[k][j][iw], dpy, py4[k][j][ie]);
					dpie = dp(py4[k][j][ie], dpy, py4[k][j][iw]);
					dpjs = dp(py4[k][js][i], dpy, py4[k][jn][i]);
					dpjn = dp(py4[k][jn][i], dpy, py4[k][js][i]);
					dpkf = dp(py4[kf][j][i], dpy, py4[kb][j][i]);
					dpkb = dp(py4[kb][j][i], dpy, py4[kf][j][i]);
					q04[k1][j1][i1] = dpy + dpiw + dpjs + dpkf;
					q04[k1][j1][i1 + 1] = dpy + dpie + dpjs + dpkf;
					q04[k1 + 1][j1][i1] = dpy + dpiw + dpjs + dpkb;
					q04[k1 + 1][j1][i1 + 1] = dpy + dpie + dpjs + dpkb;
					q04[k1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkf;
					q04[k1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkf;
					q04[k1 + 1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkb;
					q04[k1 + 1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkb;

					dpy = py5[k][j][i];
					dpiw = dp(py5[k][j][iw], dpy, py5[k][j][ie]);
					dpie = dp(py5[k][j][ie], dpy, py5[k][j][iw]);
					dpjs = dp(py5[k][js][i], dpy, py5[k][jn][i]);
					dpjn = dp(py5[k][jn][i], dpy, py5[k][js][i]);
					dpkf = dp(py5[kf][j][i], dpy, py5[kb][j][i]);
					dpkb = dp(py5[kb][j][i], dpy, py5[kf][j][i]);
					q05[k1][j1][i1] = dpy + dpiw + dpjs + dpkf;
					q05[k1][j1][i1 + 1] = dpy + dpie + dpjs + dpkf;
					q05[k1 + 1][j1][i1] = dpy + dpiw + dpjs + dpkb;
					q05[k1 + 1][j1][i1 + 1] = dpy + dpie + dpjs + dpkb;
					q05[k1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkf;
					q05[k1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkf;
					q05[k1 + 1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkb;
					q05[k1 + 1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkb;
				}
			}
		}

		nnx = (pfield->mesh->nnx)[ig+1];
		nny = (pfield->mesh->nny)[ig+1];
		nnz = (pfield->mesh->nnz)[ig+1];
		/*			py1=q01;
		py2=q02;
		py3=q03;
		py4=q04;
		py5=q05;
		*/
		for (int i = 1; i<nnz + 1; i++)
			for (int j = 1; j<nny + 1; j++)
				for (int k = 1; k<nnx + 1; k++)
				{
					py1[i][j][k] = q01[i][j][k];
					py2[i][j][k] = q02[i][j][k];
					py3[i][j][k] = q03[i][j][k];
					py4[i][j][k] = q04[i][j][k];
					py5[i][j][k] = q05[i][j][k];

				}
	}

	for (int i = 1; i<nnz + 1; i++)
	for (int j = 1; j<nny + 1; j++)
	for (int k = 1; k<nnx + 1; k++)
	{
			q31[n][i][j][k] = q31[n][i][j][k] + py1[i][j][k];
			q32[n][i][j][k] = q32[n][i][j][k] + py2[i][j][k];
			q33[n][i][j][k] = q33[n][i][j][k] + py3[i][j][k];
			q34[n][i][j][k] = q34[n][i][j][k] + py4[i][j][k];
			q35[n][i][j][k] = q35[n][i][j][k] + py5[i][j][k];
	}
}

void CRungeKutta::update(int nng)
{
	int i1, j1, k1, iw, ie, jn, js, kf, kb;
	double dpy, dpiw, dpie, dpjn, dpjs, dpkf, dpkb;

	SetMultigridRange(nng);
	for (int nn = 1; nn<pdict->nt + 1; nn++)
	{
		for (int i = 1; i<nz + 1; i++)
			for (int j = 1; j<ny + 1; j++)
				for (int k = 1; k<nx + 1; k++)
				{
					py1[i][j][k] = pfield->q11[nn][i][j][k];
					py2[i][j][k] = pfield->q12[nn][i][j][k];
					py3[i][j][k] = pfield->q13[nn][i][j][k];
					py4[i][j][k] = pfield->q14[nn][i][j][k];
					py5[i][j][k] = pfield->q15[nn][i][j][k];
					py6[i][j][k] = pfield->q16[nn][i][j][k];
				}

		for (int k = 1; k<nz + 1; k++)
		{
			/*				i1 = 2 * i - 1;
			iw = max(1, i - 1);
			ie = min(nx, i + 1);*/
			k1 = 2 * k - 1;
			kf = max(1, k - 1);
			kb = min(nz, k + 1);
			for (int j = 1; j<ny + 1; j++)
			{
				j1 = 2 * j - 1;
				js = max(1, j - 1);
				jn = min(ny, j + 1);
				for (int i = 1; i<nx + 1; i++)
				{
					/*					k1 = 2 * k - 1;
					kf = max(1, k - 1);
					kb = min(nz, k + 1);
					*/
					i1 = 2 * i - 1;
					iw = max(1, i - 1);
					ie = min(nx, i + 1);

					dpy = py1[k][j][i];
					dpiw = dp(py1[k][j][iw], dpy, py1[k][j][ie]);
					dpie = dp(py1[k][j][ie], dpy, py1[k][j][iw]);
					dpjs = dp(py1[k][js][i], dpy, py1[k][jn][i]);
					dpjn = dp(py1[k][jn][i], dpy, py1[k][js][i]);
					dpkf = dp(py1[kf][j][i], dpy, py1[kb][j][i]);
					dpkb = dp(py1[kb][j][i], dpy, py1[kf][j][i]);

					pfield->q11[nn][k1][j1][i1] = dpy + dpiw + dpjs + dpkf;
					pfield->q11[nn][k1][j1][i1 + 1] = dpy + dpie + dpjs + dpkf;
					pfield->q11[nn][k1 + 1][j1][i1] = dpy + dpiw + dpjs + dpkb;
					pfield->q11[nn][k1 + 1][j1][i1 + 1] = dpy + dpie + dpjs + dpkb;
					pfield->q11[nn][k1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkf;
					pfield->q11[nn][k1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkf;
					pfield->q11[nn][k1 + 1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkb;
					pfield->q11[nn][k1 + 1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkb;

					dpy = py2[k][j][i];
					dpiw = dp(py2[k][j][iw], dpy, py2[k][j][ie]);
					dpie = dp(py2[k][j][ie], dpy, py2[k][j][iw]);
					dpjs = dp(py2[k][js][i], dpy, py2[k][jn][i]);
					dpjn = dp(py2[k][jn][i], dpy, py2[k][js][i]);
					dpkf = dp(py2[kf][j][i], dpy, py2[kb][j][i]);
					dpkb = dp(py2[kb][j][i], dpy, py2[kf][j][i]);

					pfield->q12[nn][k1][j1][i1] = dpy + dpiw + dpjs + dpkf;
					pfield->q12[nn][k1][j1][i1 + 1] = dpy + dpie + dpjs + dpkf;
					pfield->q12[nn][k1 + 1][j1][i1] = dpy + dpiw + dpjs + dpkb;
					pfield->q12[nn][k1 + 1][j1][i1 + 1] = dpy + dpie + dpjs + dpkb;
					pfield->q12[nn][k1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkf;
					pfield->q12[nn][k1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkf;
					pfield->q12[nn][k1 + 1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkb;
					pfield->q12[nn][k1 + 1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkb;

					dpy = py3[k][j][i];
					dpiw = dp(py3[k][j][iw], dpy, py3[k][j][ie]);
					dpie = dp(py3[k][j][ie], dpy, py3[k][j][iw]);
					dpjs = dp(py3[k][js][i], dpy, py3[k][jn][i]);
					dpjn = dp(py3[k][jn][i], dpy, py3[k][js][i]);
					dpkf = dp(py3[kf][j][i], dpy, py3[kb][j][i]);
					dpkb = dp(py3[kb][j][i], dpy, py3[kf][j][i]);

					pfield->q13[nn][k1][j1][i1] = dpy + dpiw + dpjs + dpkf;
					pfield->q13[nn][k1][j1][i1 + 1] = dpy + dpie + dpjs + dpkf;
					pfield->q13[nn][k1 + 1][j1][i1] = dpy + dpiw + dpjs + dpkb;
					pfield->q13[nn][k1 + 1][j1][i1 + 1] = dpy + dpie + dpjs + dpkb;
					pfield->q13[nn][k1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkf;
					pfield->q13[nn][k1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkf;
					pfield->q13[nn][k1 + 1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkb;
					pfield->q13[nn][k1 + 1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkb;

					dpy = py4[k][j][i];
					dpiw = dp(py4[k][j][iw], dpy, py4[k][j][ie]);
					dpie = dp(py4[k][j][ie], dpy, py4[k][j][iw]);
					dpjs = dp(py4[k][js][i], dpy, py4[k][jn][i]);
					dpjn = dp(py4[k][jn][i], dpy, py4[k][js][i]);
					dpkf = dp(py4[kf][j][i], dpy, py4[kb][j][i]);
					dpkb = dp(py4[kb][j][i], dpy, py4[kf][j][i]);
					pfield->q14[nn][k1][j1][i1] = dpy + dpiw + dpjs + dpkf;
					pfield->q14[nn][k1][j1][i1 + 1] = dpy + dpie + dpjs + dpkf;
					pfield->q14[nn][k1 + 1][j1][i1] = dpy + dpiw + dpjs + dpkb;
					pfield->q14[nn][k1 + 1][j1][i1 + 1] = dpy + dpie + dpjs + dpkb;
					pfield->q14[nn][k1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkf;
					pfield->q14[nn][k1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkf;
					pfield->q14[nn][k1 + 1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkb;
					pfield->q14[nn][k1 + 1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkb;

					dpy = py5[k][j][i];
					dpiw = dp(py5[k][j][iw], dpy, py5[k][j][ie]);
					dpie = dp(py5[k][j][ie], dpy, py5[k][j][iw]);
					dpjs = dp(py5[k][js][i], dpy, py5[k][jn][i]);
					dpjn = dp(py5[k][jn][i], dpy, py5[k][js][i]);
					dpkf = dp(py5[kf][j][i], dpy, py5[kb][j][i]);
					dpkb = dp(py5[kb][j][i], dpy, py5[kf][j][i]);

					pfield->q15[nn][k1][j1][i1] = dpy + dpiw + dpjs + dpkf;
					pfield->q15[nn][k1][j1][i1 + 1] = dpy + dpie + dpjs + dpkf;
					pfield->q15[nn][k1 + 1][j1][i1] = dpy + dpiw + dpjs + dpkb;
					pfield->q15[nn][k1 + 1][j1][i1 + 1] = dpy + dpie + dpjs + dpkb;
					pfield->q15[nn][k1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkf;
					pfield->q15[nn][k1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkf;
					pfield->q15[nn][k1 + 1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkb;
					pfield->q15[nn][k1 + 1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkb;

					dpy = py6[k][j][i];
					dpiw = dp(py6[k][j][iw], dpy, py6[k][j][ie]);
					dpie = dp(py6[k][j][ie], dpy, py6[k][j][iw]);
					dpjs = dp(py6[k][js][i], dpy, py6[k][jn][i]);
					dpjn = dp(py6[k][jn][i], dpy, py6[k][js][i]);
					dpkf = dp(py6[kf][j][i], dpy, py6[kb][j][i]);
					dpkb = dp(py6[kb][j][i], dpy, py6[kf][j][i]);

					pfield->q16[nn][k1][j1][i1] = dpy + dpiw + dpjs + dpkf;
					pfield->q16[nn][k1][j1][i1 + 1] = dpy + dpie + dpjs + dpkf;
					pfield->q16[nn][k1 + 1][j1][i1] = dpy + dpiw + dpjs + dpkb;
					pfield->q16[nn][k1 + 1][j1][i1 + 1] = dpy + dpie + dpjs + dpkb;
					pfield->q16[nn][k1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkf;
					pfield->q16[nn][k1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkf;
					pfield->q16[nn][k1 + 1][j1 + 1][i1] = dpy + dpiw + dpjn + dpkb;
					pfield->q16[nn][k1 + 1][j1 + 1][i1 + 1] = dpy + dpie + dpjn + dpkb;
				}
			}
		}
	}
}
void CRungeKutta::FMGCycle(double time_begin,double time_end)
{

	double rmsm, rmsmmax;

	while(V_Cycle())
	{
		
		if(myid==0)LogInfo::Log("The "+Int_to_string(V_stage)+"th stage of  fmg computing start");

		SetMultigridRange(V_stage);
		pgrid->ComputeGeoData(V_stage, pdict->lbb);
		pfield->SetPhysicalBoundary(V_stage, pdict->vxx, pdict->vrr, pdict->vtt, pdict->vee, pdict->ht, pdict->pt, pgrid);

		while (Iteration())
		{
			
			if (myid == 0)LogInfo::Log("    the "+Int_to_string(iter)+"th iteration in stage "+Int_to_string(V_stage)+" begin");

			rmsm = -11;
			
			while (Vitual_timeLoop())
			{
				if (myid == 0)LogInfo::Log("         the "+Int_to_string(n)+"th virtual time step calculation in iteration "+Int_to_string(iter)+" switch on");

				pgrid->SetlayerGeoBoundary(V_stage);
				pfield->GetlayerPhysicalBoundary(n, V_stage);
				ResetP2h();
				timedrivationCmpute();
				RungeKutta_I(V_stage);
				residual(rmsm);
				RK_I();
				InitPreField();
				timl = pdict->cfl*0.125;
				UpdateFieldValue(ta, timl, 0);
				if (V_stage == (pfield->mesh)->ng)
				{
					RungeKutta_SA();
				}
				for (int nn = 1; nn<pdict->nt + 1; nn++)
				for (int k = 1; k<nz + 1; k++)
				for (int j = 1; j<ny + 1; j++)
				for (int i = 1; i<nx + 1; i++)
				{
					q31[nn][k][j][i] = pfield->q11[nn][k][j][i];
					q32[nn][k][j][i] = pfield->q12[nn][k][j][i];
					q33[nn][k][j][i] = pfield->q13[nn][k][j][i];
					q34[nn][k][j][i] = pfield->q14[nn][k][j][i];
					q35[nn][k][j][i] = pfield->q15[nn][k][j][i];
					q36[nn][k][j][i] = pfield->q16[nn][k][j][i];
				}

			
				for (int ign = V_stage - 1; ign >= 1; ign--)
				{
					
						pgrid->SetlayerGeoBoundary(ign);
						SetMultigridRange(ign);
						RestrictionOp(ign);
						pfield->GetlayerPhysicalBoundary(n, ign);
						P2hSetting();
						RungeKutta_II(  ign);
						
						UpwardCorrection(ign, V_stage);
						RK_C_I();
						

						InitPreField();
						timl = pdict->cfl*0.1;
						UpdateFieldValue(ta, timl, 0);
						
						UpwardCorrection(ign, V_stage);
					
				}
				pgrid->SetlayerGeoBoundary(V_stage);
				SetMultigridRange( V_stage);

				for (int i = 1; i<pdict->nt + 1; i++)
				for (int j = 1; j<nz + 1; j++)
				for (int k = 1; k<ny + 1; k++)
				for (int nn = 1; nn<nx + 1; nn++)
				{
						pfield->q11[i][j][k][nn] = q31[i][j][k][nn];
						pfield->q12[i][j][k][nn] = q32[i][j][k][nn];
						pfield->q13[i][j][k][nn] = q33[i][j][k][nn];
						pfield->q14[i][j][k][nn] = q34[i][j][k][nn];
						pfield->q15[i][j][k][nn] = q35[i][j][k][nn];
						pfield->q16[i][j][k][nn] = q36[i][j][k][nn];
				}
			/*	double s = 0.0;
				for (int i = 1; i<pdict->nt + 1; i++)
					for (int j = 1; j<nz + 1; j++)
						for (int k = 1; k<ny + 1; k++)
							for (int n = 1; n<nx + 1; n++)
							{
								s += pfield->q11[i][j][k][n];
							}
				f1 << setw(15) << setprecision(10) << s << endl;*/
				if (myid == 0)LogInfo::Log("         the "+Int_to_string(n)+"th virtual time step calculation in iteration "+Int_to_string(iter)+" switch off");

			}

			time_end = MPI_Wtime();
	
			//输出残差文件
			Write_itertime( iteration, rmsm, time_end, time_begin);

			MPI_Allreduce(&rmsm, &rmsmmax, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

			if (rmsmmax == pdict->rmsm0)
				goto label;
			flow(iteration);

			if (myid == 0)LogInfo::Log("    the "+Int_to_string(iter)+"th iteration in stage "+Int_to_string(V_stage)+" is over");

		}
		if (V_stage < pdict->ng)
		{
			update(V_stage);
			if(myid==0)LogInfo::Log("The "+Int_to_string(V_stage)+"th stage of  fmg computing ends");
		}
		else
		{
			if(myid==0)LogInfo::Log("The "+Int_to_string(V_stage)+"th stage of  fmg computing ends");
			if(myid==0)LogInfo::Log("simulation begin output files");			
			output();
		}
	}
label:;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //



