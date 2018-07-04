/*---------------------------------------------------------------------------*\


\*---------------------------------------------------------------------------*/

#ifndef Mesh_H
#define Mesh_H

#include "Memory.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //




/*---------------------------------------------------------------------------*\
                           Class Mesh Declaration
\*---------------------------------------------------------------------------*/


class CMesh:
public CMemory<double>
{
private:
	string file1;


public:


	int nxm1,nym1,nzm1,ng;

	double*** dmini;
	double*** xf, ***yf, ***zf;
	int* nnx, *nny, *nnz, *nib, *nit, *njb, *njt;
	double*** x,***y,***z;



public:

    // Constructors

        CMesh(string s1,int Gridnum)
	:file1(s1),nxm1(0),nym1(0),nzm1(0),ng(Gridnum)
        {

		nnx=new int[ng+1];
		nny=new int[ng+1];
		nnz=new int[ng+1];
		nib=new int[ng+1];
		nit=new int[ng+1];
		njb=new int[ng+1];
		njt=new int[ng+1];		
	}

	~CMesh()
	{
		Free(dmini);
		Free(xf);
		Free(yf);
		Free(zf);
		Free(x);
		Free(y);
		Free(z);
		delete[] nnx;
		delete[] nny;
		delete[] nnz;
		delete[] nib;
		delete[] nit;
		delete[] njb;
		delete[] njt;
	}
	bool getOrgData();
	void GenerateMesh(int myid,int lbb);
	bool SetMultiGrid(int ibm,int jbm,int itm,int jtm);
	double*** min_distance();
	bool SetlayerMesh(int LayerRange);

private:
	bool initvector(int x,int y,int z);

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
