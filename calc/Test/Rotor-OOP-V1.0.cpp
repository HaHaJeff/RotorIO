#include"headfiles.h"

int main(int argc,char** argv)
{
	MPI_Init(&argc, &argv);
	int myid, numprocs;
	int nng0 = 1;
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);

	if(myid==0)LogInfo::Log("Program begin");

	double time_begin= MPI_Wtime();
	double time_end;
	
	if(myid==0)LogInfo::Log("Start reading configuration");
	
	//读入配置参数
	CDictionary *pdict = new CDictionary("../input/ini3.dat");
	pdict->readConfig();

	if(myid==0)LogInfo::Log("Configuration reading ends");

	if(myid==0)LogInfo::Log("Start fetch mesh data");

	//读入网格文件
	CMesh *mesh = new CMesh("../input/grid.dat",pdict->ng);
	mesh->getOrgData();
	mesh->GenerateMesh(myid,pdict->lbb);
	mesh->SetMultiGrid(pdict->ibm, pdict->jbm, pdict->itm, pdict->jtm);
	mesh->SetlayerMesh(nng0);
	if(myid==0)LogInfo::Log("Compute minimum distance");
	mesh->min_distance();

	if(myid==0)LogInfo::Log("Mesh preprocessing ends");


	if(myid==0)LogInfo::Log("Initialize field value from physics");

	CMultiGrid *pgrid = new CMultiGrid(mesh, pdict->ng);
	CGeoField *pfield = new CGeoField(mesh, pdict->nt);
	pfield->SetPhysicalBoundary(nng0, pdict->vxx, pdict->vrr, pdict->vtt, pdict->vee, pdict->ht, pdict->pt);
	pfield->FieldInitialization(nng0, pdict->ma, pdict->cvl0);

	if(myid==0)LogInfo::Log("Filed initialization ends");

	if(myid==0)LogInfo::Log("Full Multi-Grid Cycle start");

	CRungeKutta *kutta = new CRungeKutta(pgrid, pdict, pfield, numprocs,myid);
	kutta->SetTimespectrum();
	kutta->FMGCycle(time_begin,time_end);

	MPI_Finalize();
	if(myid==0)LogInfo::Log("Program finished");
	return 0;
}
