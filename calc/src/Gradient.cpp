/*---------------------------------------------------------------------------*\


\*---------------------------------------------------------------------------*/


#include "Gradient.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void CGradient::ResetGradient()
{
	for (int i = 1; i<16; i++)
	for (int l = 0; l<nz + 2; l++)
	for (int k = 0; k<ny + 2; k++)
	for (int j = 1; j<nx + 2; j++)
	{
		gradfi[i][l][k][j] = 0;
	}
	//Memset(gradfi[0][0][0], 0, 16*(nz + 2)*(ny + 2)*(nx + 2)*sizeof(double));
	for (int i = 1; i<16; i++)
	for (int l = 0; l<nz + 2; l++)
	for (int k = 1; k<ny + 2; k++)
	for (int j = 0; j<nx + 2; j++)
	{
		gradfj[i][l][k][j] = 0;
	}
	//Memset(gradfj[0][0][0], 0, 16*(nz + 2)*(ny + 2)*(nx + 2)*sizeof(double));
	for (int i = 1; i<16; i++)
	for (int l = 1; l<nz + 2; l++)
	for (int k = 0; k<ny + 2; k++)
	for (int j = 0; j<nx + 2; j++)
	{
		gradfk[i][l][k][j] = 0;
	}
	//Memset(gradfk[0][0][0], 0, 16*(nz + 2)*(ny + 2)*(nx + 2)*sizeof(double));
}
void CGradient::gradsfaceI(double***& si, double***& sj,double***& sk,double***& q,double****& Idqd, int m)
{
	double sx, fg, qav, rvol, sx1, sx2, qav1, qav2;
	for (int k = 0; k<nz + 2; k++)
	for (int j = 0; j<ny + 2; j++)
	for (int i = 1; i<nx + 2; i++)
	{
		Idqd[m][k][j][i] = 0;
	}
	//Memset(Idqd[m][0][0], 0, (nz + 2)*(ny + 2)*(nx + 2)*sizeof(double));
	//

	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	{
		for (int i = 1; i<nx + 1; i++)
		{
			sx = 0.50*(si[k][j][i] + si[k][j][i + 1]);
			fg = q[k][j][i] * sx;
			Idqd[m][k][j][i] = Idqd[m][k][j][i] - fg;
			Idqd[m][k][j][i + 1] = Idqd[m][k][j][i + 1] + fg;
		}
		sx = si[k][j][1];
		fg = 0.50*(q[k][j][0] + q[k][j][1])*sx;
		Idqd[m][k][j][1] = Idqd[m][k][j][1] + fg;
		sx = si[k][j][nx + 1];
		fg = 0.50*(q[k][j][nx] + q[k][j][nx + 1])*sx;
		Idqd[m][k][j][nx + 1] = Idqd[m][k][j][nx + 1] - fg;
	}
	//
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 2; j++)
	{
		for (int i = 2; i<nx + 1; i++)
		{
			sx1 = 0.50*sj[k][j][i];
			qav1 = 0.50*(q[k][j][i] + q[k][j - 1][i]);
			sx2 = 0.50*sj[k][j][i - 1];
			qav2 = 0.50*(q[k][j][i - 1] + q[k][j - 1][i - 1]);
			fg = sx1*qav1 + sx2*qav2;
			Idqd[m][k][j][i] = Idqd[m][k][j][i] + fg;
			Idqd[m][k][j - 1][i] = Idqd[m][k][j - 1][i] - fg;
		}
		sx = sj[k][j][1];
		qav = 0.50*(q[k][j][1] + q[k][j - 1][1]);
		fg = qav*sx;
		Idqd[m][k][j][1] = Idqd[m][k][j][1] + fg;
		Idqd[m][k][j - 1][1] = Idqd[m][k][j - 1][1] - fg;
		sx = sj[k][j][nx];
		qav = 0.50*(q[k][j][nx] + q[k][j - 1][nx]);
		fg = qav*sx;
		Idqd[m][k][j][nx + 1] = Idqd[m][k][j][nx + 1] + fg;
		Idqd[m][k][j - 1][nx + 1] = Idqd[m][k][j - 1][nx + 1] - fg;
	}
	//
	for (int k = 1; k<nz + 2; k++)
	for (int j = 1; j<ny + 1; j++)
	{
		for (int i = 2; i<nx + 1; i++)
		{
			sx1 = 0.50*sk[k][j][i];
			qav1 = 0.50*(q[k][j][i] + q[k - 1][j][i]);
			sx2 = 0.50*sk[k][j][i - 1];
			qav2 = 0.50*(q[k][j][i - 1] + q[k - 1][j][i - 1]);
			fg = sx1*qav1 + sx2*qav2;
			Idqd[m][k][j][i] = Idqd[m][k][j][i] + fg;
			Idqd[m][k - 1][j][i] = Idqd[m][k - 1][j][i] - fg;
		}
		sx = sk[k][j][1];
		qav = 0.50*(q[k][j][1] + q[k - 1][j][1]);
		fg = qav*sx;
		Idqd[m][k][j][1] = Idqd[m][k][j][1] + fg;
		Idqd[m][k - 1][j][1] = Idqd[m][k - 1][j][1] - fg;
		sx = sk[k][j][nx];
		qav = 0.50*(q[k][j][nx] + q[k - 1][j][nx]);
		fg = qav*sx;
		Idqd[m][k][j][nx + 1] = Idqd[m][k][j][nx + 1] + fg;
		Idqd[m][k - 1][j][nx + 1] = Idqd[m][k - 1][j][nx + 1] - fg;
	}

	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	{
		for (int i = 2; i<nx + 1; i++)
		{
			rvol = 2.0 / (pgrid->vv[k][j][i] + pgrid->vv[k][j][i - 1]);
			Idqd[m][k][j][i] = Idqd[m][k][j][i] * rvol;
		}
		rvol = 1.0 / pgrid->vv[k][j][1];
		Idqd[m][k][j][1] = Idqd[m][k][j][1] * rvol;
		rvol = 1.0 / pgrid->vv[k][j][nx];
		Idqd[m][k][j][nx + 1] = Idqd[m][k][j][nx + 1] * rvol;
	}
}

void CGradient::gradsfaceJ(double***& si, double***& sj,double***& sk,double***& q, double****& Jdqd, int m)
{
	double sx, fg, qav, rvol, sx1, sx2, qav1, qav2,fg1,rvol1;
	for (int k = 0; k<nz + 2; k++)
	for (int j = 1; j<ny + 2; j++)
	for (int i = 0; i<nx + 2; i++)
	{
		Jdqd[m][k][j][i] = 0;
	}
	for (int k = 1; k<nz + 1; k++)
	for (int i = 1; i<nx + 1; i++)
	{
		sx = sj[k][1][i];
		fg = 0.50*(q[k][0][i] + q[k][1][i])*sx;
		sx1 = sj[k][ny + 1][i];
		fg1 = 0.50*(q[k][ny][i] + q[k][ny + 1][i])*sx1;
		Jdqd[m][k][1][i] = Jdqd[m][k][1][i] + fg;
		Jdqd[m][k][ny + 1][i] = Jdqd[m][k][ny + 1][i] - fg1;
	}
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
        for (int i = 1; i<nx + 1; i++)
        {
                sx = 0.50*(sj[k][j][i] + sj[k][j + 1][i]);
                fg = q[k][j][i] * sx;
                Jdqd[m][k][j][i] = Jdqd[m][k][j][i] - fg;
                Jdqd[m][k][j + 1][i] = Jdqd[m][k][j + 1][i] + fg;
        }


	for (int k = 1; k<nz + 1; k++)
	for (int i = 1; i<nx + 2; i++)
	{
		sx = si[k][1][i];
		qav = 0.50*(q[k][1][i] + q[k][1][i-1]);
		fg = qav*sx;
		sx1 = si[k][ny][i];
		qav1 = 0.50*(q[k][ny][i] + q[k][ny][i-1]);
		fg1 = qav1*sx1;
		Jdqd[m][k][1][i] = Jdqd[m][k][1][i] + fg;
                Jdqd[m][k][1][i-1] = Jdqd[m][k][1][i-1] - fg;
		Jdqd[m][k][ny + 1][i] = Jdqd[m][k][ny + 1][i] + fg1;
		Jdqd[m][k][ny + 1][i-1] = Jdqd[m][k][ny + 1][i-1] - fg1;
	}
	for (int k = 1; k<nz + 1; k++)
	for (int j = 2; j<ny + 1; j++) 
        for (int i = 1; i<nx + 2; i++)
        {
                sx1 = 0.50*si[k][j][i];
                qav1 = 0.50*(q[k][j][i] + q[k][j][i-1]);
                sx2 = 0.50*si[k][j - 1][i];
                qav2 = 0.50*(q[k][j - 1][i] + q[k][j - 1][i-1]);
                fg = sx1*qav1 + sx2*qav2;
                Jdqd[m][k][j][i] = Jdqd[m][k][j][i] + fg;
                Jdqd[m][k][j][i-1] = Jdqd[m][k][j][i-1] - fg;
        }


	for (int k = 1; k<nz + 2; k++)
	for (int i = 1; i<nx + 1; i++)
	{
		sx = sk[k][1][i];
		qav = 0.50*(q[k][1][i] + q[k-1][1][i]);
		fg = qav*sx;
		sx1 = sk[k][ny][i];
		qav1 = 0.50*(q[k][ny][i] + q[k-1][ny][i]);
		fg1 = qav1*sx1;
		Jdqd[m][k][1][i] = Jdqd[m][k][1][i] + fg;
                Jdqd[m][k-1][1][i] = Jdqd[m][k-1][1][i] - fg;
		Jdqd[m][k][ny + 1][i] = Jdqd[m][k][ny + 1][i] + fg1;
		Jdqd[m][k-1][ny + 1][i] = Jdqd[m][k-1][ny + 1][i] - fg1;
	}
	for (int k = 1; k<nz + 2; k++)
	for (int j = 2; j<ny + 1; j++)
        for (int i = 1; i<nx + 1; i++)
        {
           sx1 = 0.50*sk[k][j][i];
           qav1 = 0.50*(q[k][j][i] + q[k-1][j][i]);
           sx2 = 0.50*sk[k][j - 1][i];
           qav2 = 0.50*(q[k][j - 1][i] + q[k-1][j - 1][i]);
           fg = sx1*qav1 + sx2*qav2;
           Jdqd[m][k][j][i] = Jdqd[m][k][j][i] + fg;
           Jdqd[m][k-1][j][i] = Jdqd[m][k-1][j][i] - fg;
          }

	for (int k = 1; k<nz + 1; k++)
	for (int i = 1; i<nx + 1; i++)
	{
		rvol = 1.0 / pgrid->vv[k][1][i];
		rvol1 = 1.0 / pgrid->vv[k][ny][i];
		Jdqd[m][k][1][i] = Jdqd[m][k][1][i] * rvol;
		Jdqd[m][k][ny + 1][i] = Jdqd[m][k][ny + 1][i] * rvol1;
	}
	for (int k = 1; k<nz + 1; k++)
	for (int j = 2; j<ny + 1; j++)
        for (int i = 1; i<nx + 1; i++)
        {
            	rvol = 2.0 / (pgrid->vv[k][j][i] + pgrid->vv[k][j - 1][i]);
                Jdqd[m][k][j][i] = Jdqd[m][k][j][i] * rvol;
        }
}

void CGradient::gradsfaceK(double***& si, double***& sj,double***& sk,double***& q, double****& Kdqd, int m)
{
	double sx, fg, qav, rvol, sx1, sx2, qav1, qav2,fg1,rvol1;
	for (int k = 1; k<nz + 2; k++)
	for (int j = 0; j<ny + 2; j++)
	for (int i = 0; i<nx + 2; i++)
	{
		Kdqd[m][k][j][i] = 0;
	}
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
	{

		sx = sk[1][j][i];
		fg = 0.50*(q[0][j][i] + q[1][j][i])*sx;
		sx1 = sk[nz+1][j][i];
		fg1 = 0.50*(q[nz][j][i] + q[nz+1][j][i])*sx1;	
		Kdqd[m][1][j][i] = Kdqd[m][1][j][i] + fg;
		Kdqd[m][nz+1][j][i] = Kdqd[m][nz+1][j][i] - fg1;
	}
	for(int k=1;k<nz+1;k++)
        for(int j=1;j<ny + 1; j++)
        for (int i = 1; i<nx + 1; i++)
        {
        	sx = 0.50*(sk[k][j][i] + sk[k+1][j][i]);
        	fg = q[k][j][i] * sx;
		Kdqd[m][k][j][i] = Kdqd[m][k][j][i] - fg;
		Kdqd[m][k+1][j][i] = Kdqd[m][k+1][j][i] + fg;
        }


	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 2; i++)
	{
		sx = si[1][j][i];
		qav = 0.50*(q[1][j][i] + q[1][j][i-1]);
		fg = qav*sx;
		sx1 = si[nz][j][i];
		qav1 = 0.50*(q[nz][j][i] + q[nz][j][i-1]);
		fg1 = qav1*sx1;
		Kdqd[m][1][j][i] = Kdqd[m][1][j][i] + fg;
                Kdqd[m][1][j][i-1] = Kdqd[m][1][j][i-1] - fg;
		Kdqd[m][nz+1][j][i] = Kdqd[m][nz+1][j][i] + fg1;
		Kdqd[m][nz+1][j][i-1] = Kdqd[m][nz+1][j][i-1] - fg1;
	}
	for(int k=2;k<nz+1;k++)
	for (int j = 1; j<ny + 1; j++)
        for (int i = 1; i<nx + 2; i++)
        {
                sx1 = 0.50*si[k][j][i];
                qav1 = 0.50*(q[k][j][i] + q[k][j][i-1]);
                sx2 = 0.50*si[k-1][j][i];
                qav2 = 0.50*(q[k-1][j][i] + q[k - 1][j][i - 1]);
                fg = sx1*qav1 + sx2*qav2;
                Kdqd[m][k][j][i] = Kdqd[m][k][j][i] + fg;
                Kdqd[m][k][j][i-1] = Kdqd[m][k][j][i-1] - fg;
        }

	for (int j = 1; j<ny + 2; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		sx = sj[1][j][i];
		qav = 0.50*(q[1][j][i] + q[1][j - 1][i]);
		fg = qav*sx;
		sx1 = sj[nz][j][i];
		qav1 = 0.50*(q[nz][j][i] + q[nz][j - 1][i]);
		fg1 = qav1*sx1;			
		Kdqd[m][1][j][i] = Kdqd[m][1][j][i] + fg;
                Kdqd[m][1][j - 1][i] = Kdqd[m][1][j - 1][i] - fg;
		Kdqd[m][nz+1][j][i] = Kdqd[m][nz+1][j][i] + fg1;
		Kdqd[m][nz+1][j - 1][i] = Kdqd[m][nz+1][j - 1][i] - fg1;
	}
	for (int k = 2; k<nz + 1; k++)
	for (int j = 1; j<ny + 2; j++)
        for (int i = 1; i<nx + 1; i++)
        {
            sx1 = 0.50*sj[k][j][i];
            qav1 = 0.50*(q[k][j][i] + q[k][j - 1][i]);
            sx2 = 0.50*sj[k-1][j][i];
            qav2 = 0.50*(q[k-1][j][i] + q[k-1][j - 1][i]);
            fg = sx1*qav1 + sx2*qav2;
            Kdqd[m][k][j][i] = Kdqd[m][k][j][i] + fg;
            Kdqd[m][k][j - 1][i] = Kdqd[m][k][j - 1][i] - fg;
        }

	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		rvol = 1.0 / pgrid->vv[1][j][i];
		rvol1 = 1.0 / pgrid->vv[nz][j][i];
		Kdqd[m][1][j][i] = Kdqd[m][1][j][i] * rvol;
		Kdqd[m][nz+1][j][i] = Kdqd[m][nz+1][j][i] * rvol1;
	}
	for (int k = 2; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
        for (int i = 1; i<nx + 1; i++)
        {
                rvol = 2.0 / (pgrid->vv[k][j][i] + pgrid->vv[k-1][j][i]);
                Kdqd[m][k][j][i] = Kdqd[m][k][j][i] * rvol;
        }
}


void CGradient::gradsface()
{
	ResetGradient();
	gradsfaceI(pgrid->s2x, pgrid->s3x, pgrid->s1x, pfield->pvx, gradfi, 1);
	gradsfaceI(pgrid->s2y, pgrid->s3y, pgrid->s1y, pfield->pvx, gradfi, 2);
	gradsfaceI(pgrid->s2z, pgrid->s3z, pgrid->s1z, pfield->pvx, gradfi, 3);
	gradsfaceI(pgrid->s2x, pgrid->s3x, pgrid->s1x, pfield->pvy, gradfi, 4);
	gradsfaceI(pgrid->s2y, pgrid->s3y, pgrid->s1y, pfield->pvy, gradfi, 5);
	gradsfaceI(pgrid->s2z, pgrid->s3z, pgrid->s1z, pfield->pvy, gradfi, 6);
	gradsfaceI(pgrid->s2x, pgrid->s3x, pgrid->s1x, pfield->pvz, gradfi, 7);
	gradsfaceI(pgrid->s2y, pgrid->s3y, pgrid->s1y, pfield->pvz, gradfi, 8);
	gradsfaceI(pgrid->s2z, pgrid->s3z, pgrid->s1z, pfield->pvz, gradfi, 9);
	gradsfaceI(pgrid->s2x, pgrid->s3x, pgrid->s1x, pfield->t, gradfi, 10);
	gradsfaceI(pgrid->s2y, pgrid->s3y, pgrid->s1y, pfield->t, gradfi, 11);
	gradsfaceI(pgrid->s2z, pgrid->s3z, pgrid->s1z, pfield->t, gradfi, 12);
	gradsfaceJ(pgrid->s2x, pgrid->s3x, pgrid->s1x, pfield->pvx, gradfj, 1);
	gradsfaceJ(pgrid->s2y, pgrid->s3y, pgrid->s1y, pfield->pvx, gradfj, 2);
	gradsfaceJ(pgrid->s2z, pgrid->s3z, pgrid->s1z, pfield->pvx, gradfj, 3);
	gradsfaceJ(pgrid->s2x, pgrid->s3x, pgrid->s1x, pfield->pvy, gradfj, 4);
	gradsfaceJ(pgrid->s2y, pgrid->s3y, pgrid->s1y, pfield->pvy, gradfj, 5);
	gradsfaceJ(pgrid->s2z, pgrid->s3z, pgrid->s1z, pfield->pvy, gradfj, 6);
	gradsfaceJ(pgrid->s2x, pgrid->s3x, pgrid->s1x, pfield->pvz, gradfj, 7);
	gradsfaceJ(pgrid->s2y, pgrid->s3y, pgrid->s1y, pfield->pvz, gradfj, 8);
	gradsfaceJ(pgrid->s2z, pgrid->s3z, pgrid->s1z, pfield->pvz, gradfj, 9);
	gradsfaceJ(pgrid->s2x, pgrid->s3x, pgrid->s1x, pfield->t, gradfj, 10);
	gradsfaceJ(pgrid->s2y, pgrid->s3y, pgrid->s1y, pfield->t, gradfj, 11);
	gradsfaceJ(pgrid->s2z, pgrid->s3z, pgrid->s1z, pfield->t, gradfj, 12);
	gradsfaceK(pgrid->s2x, pgrid->s3x, pgrid->s1x, pfield->pvx, gradfk, 1);
	gradsfaceK(pgrid->s2y, pgrid->s3y, pgrid->s1y, pfield->pvx, gradfk, 2);
	gradsfaceK(pgrid->s2z, pgrid->s3z, pgrid->s1z, pfield->pvx, gradfk, 3);
	gradsfaceK(pgrid->s2x, pgrid->s3x, pgrid->s1x, pfield->pvy, gradfk, 4);
	gradsfaceK(pgrid->s2y, pgrid->s3y, pgrid->s1y, pfield->pvy, gradfk, 5);
	gradsfaceK(pgrid->s2z, pgrid->s3z, pgrid->s1z, pfield->pvy, gradfk, 6);
	gradsfaceK(pgrid->s2x, pgrid->s3x, pgrid->s1x, pfield->pvz, gradfk, 7);
	gradsfaceK(pgrid->s2y, pgrid->s3y, pgrid->s1y, pfield->pvz, gradfk, 8);
	gradsfaceK(pgrid->s2z, pgrid->s3z, pgrid->s1z, pfield->pvz, gradfk, 9);
	gradsfaceK(pgrid->s2x, pgrid->s3x, pgrid->s1x, pfield->t, gradfk, 10);
	gradsfaceK(pgrid->s2y, pgrid->s3y, pgrid->s1y, pfield->t, gradfk, 11);
	gradsfaceK(pgrid->s2z, pgrid->s3z, pgrid->s1z, pfield->t, gradfk, 12);
}

void CGradient::gradscentre(int direction, double***& q, double****& dqd, int m)
{
	double si, flu;
	for (int k = 0; k<nz + 2; k++)
	for (int j = 0; j<ny + 2; j++)
	for (int i = 0; i<nx + 2; i++)
	{
		dqd[m][k][j][i] = 0;
	}

	//********x
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 2; i++)
	{
		if (direction == 1)
			si = pgrid->s2x[k][j][i];
		else if (direction == 2)
			si = pgrid->s2y[k][j][i];
		else if (direction == 3)
			si = pgrid->s2z[k][j][i];
		flu = 0.50*(q[k][j][i] + q[k][j][i - 1])*si;
		dqd[m][k][j][i] = dqd[m][k][j][i] + flu;
		dqd[m][k][j][i - 1] = dqd[m][k][j][i - 1] - flu;
	}
	//*********y
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 2; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		if (direction == 1)
			si = pgrid->s3x[k][j][i];
		else if (direction == 2)
			si = pgrid->s3y[k][j][i];
		else if (direction == 3)
			si = pgrid->s3z[k][j][i];

		flu = 0.50*(q[k][j - 1][i] + q[k][j][i])*si;
		dqd[m][k][j][i] = dqd[m][k][j][i] + flu;
		dqd[m][k][j - 1][i] = dqd[m][k][j - 1][i] - flu;
	}
	//*********z
	for (int k = 1; k<nz + 2; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		if (direction == 1)
			si = pgrid->s1x[k][j][i];
		else if (direction == 2)
			si = pgrid->s1y[k][j][i];
		else if (direction == 3)
			si = pgrid->s1z[k][j][i];
		flu = 0.50*(q[k][j][i] + q[k - 1][j][i])*si;
		dqd[m][k][j][i] = dqd[m][k][j][i] + flu;
		dqd[m][k - 1][j][i] = dqd[m][k - 1][j][i] - flu;
	}
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		dqd[m][k][j][i] = dqd[m][k][j][i] / pgrid->vv[k][j][i];
	}
}
void CGradient::dsdt(double****&q, double****& s, int m,double rpm)
{
	double vf, rf, qq1, flu;
	double vx, vy, vz;
	for (int k = 0; k<nz + 2; k++)
	for (int j = 0; j<ny + 2; j++)
	for (int i = 0; i<nx + 2; i++)
	{
		s[m][k][j][i] = 0;
	}
	//******************x********************
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 2; i++)
	{
		vx = 0.50*(pfield->pvx[k][j][i] + pfield->pvx[k][j][i - 1]);
		vy = 0.50*(pfield->pvy[k][j][i] + pfield->pvy[k][j][i - 1]);
		vz = 0.50*(pfield->pvz[k][j][i] + pfield->pvz[k][j][i - 1]);
		vf = vx*pgrid->s2x[k][j][i] + vy*pgrid->s2y[k][j][i] + vz*pgrid->s2z[k][j][i];
		rf = rpm*pgrid->zz02[k][j][i] * pgrid->s2y[k][j][i] - rpm*pgrid->yy02[k][j][i] * pgrid->s2z[k][j][i];
		vf = vf + rf;
		qq1 = 0.50*(q[m][k][j][i] + q[m][k][j][i - 1]);
		flu = qq1*vf;
		s[m][k][j][i] = s[m][k][j][i] + flu;
		s[m][k][j][i - 1] = s[m][k][j][i - 1] - flu;
	}
	//******************y********************
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 2; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		vx = 0.50*(pfield->pvx[k][j][i] + pfield->pvx[k][j - 1][i]);
		vy = 0.50*(pfield->pvy[k][j][i] + pfield->pvy[k][j - 1][i]);
		vz = 0.50*(pfield->pvz[k][j][i] + pfield->pvz[k][j - 1][i]);
		vf = vx*pgrid->s3x[k][j][i] + vy*pgrid->s3y[k][j][i] + vz*pgrid->s3z[k][j][i];
		rf = rpm*pgrid->zz03[k][j][i] * pgrid->s3y[k][j][i] - rpm*pgrid->yy03[k][j][i] * pgrid->s3z[k][j][i];
		vf = vf + rf;
		qq1 = 0.50*(q[m][k][j][i] + q[m][k][j - 1][i]);
		flu = qq1*vf;
		s[m][k][j][i] = s[m][k][j][i] + flu;
		s[m][k][j - 1][i] = s[m][k][j - 1][i] - flu;
	}
	//*****************z********************
	for (int k = 1; k<nz + 2; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		vx = 0.50*(pfield->pvx[k][j][i] + pfield->pvx[k - 1][j][i]);
		vy = 0.50*(pfield->pvy[k][j][i] + pfield->pvy[k - 1][j][i]);
		vz = 0.50*(pfield->pvz[k][j][i] + pfield->pvz[k - 1][j][i]);
		vf = vx*pgrid->s1x[k][j][i] + vy*pgrid->s1y[k][j][i] + vz*pgrid->s1z[k][j][i];
		rf = rpm*pgrid->zz01[k][j][i] * pgrid->s1y[k][j][i] - rpm*pgrid->yy01[k][j][i] * pgrid->s1z[k][j][i];
		vf = vf + rf;
		qq1 = 0.50*(q[m][k][j][i] + q[m][k - 1][j][i]);
		flu = qq1*vf;
		s[m][k][j][i] = s[m][k][j][i] + flu;
		s[m][k - 1][j][i] = s[m][k - 1][j][i] - flu;
	}
	for (int k = 1; k<nz + 1; k++)
	for (int j = 1; j<ny + 1; j++)
	for (int i = 1; i<nx + 1; i++)
	{
		s[m][k][j][i] = s[m][k][j][i] / pgrid->vv[k][j][i];
	}
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
