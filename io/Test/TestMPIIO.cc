#include <strings.h> 
#include <string.h>
#include <stdio.h> 
#include "IOStrategy/MPIIO.h"

void printarr(double*** data, int dims[], char *str);
double*** allocarr(int dims[]);
void freearr(double**** data);
void initarr(double**** arr, int dims[],  int rank);

int main(int argc, char** argv) {

	int rank = 0, size = 0;
	MPI_Comm comm = MPI_COMM_WORLD;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	MPI_Comm_split(MPI_COMM_WORLD, rank / 1, rank % 1, &comm);
	//get newrank and newsize
	int newrank = 0, newsize = 0;
	MPI_Comm_rank(comm, &newrank);
	MPI_Comm_size(comm, &newsize);


	char filename[5] = {};
	snprintf(filename, sizeof(int), "%d", rank);

	
	int dims[3] = {4, 4, 4};

	double ***arr = allocarr(dims);
	
	initarr(&arr, dims, rank);

	Data_3D data(&arr[0][0][0], dims[0], dims[1], dims[2], filename);

	MPIIO mpiio(comm);

	MPIIOStrategy* strategy = mpiio.GetIOStrategy( static_cast<TYPE>(0) );

	char buf[32] = {};
	int len = 32;

//	MPIIOStrategyOneFilePerProcessAllWrite io(comm);

/*
	MPI_File fh;
	MPI_File_open(comm, "3", MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &fh);
	int arr1[] = {1,2,3,4};
	MPI_File_write_all(fh, &arr1, 3, MPI_INT, NULL);
	if (rank == 0) {
		perror("raw open");
	}
*/

	int ret = 0;
	ret = strategy->Open(filename);
//	ret = io.Open(filename);
	int count = 0;
//	count = io.Write(data);
	count = strategy->Write(data);

	strategy->Close();
	delete(strategy);
	strategy = NULL;
	freearr(&arr);

//	io.Close();
	MPI_Finalize();

	return 0;
}


void printfarr(double*** data, int dims[], int rank) {
	fprintf(stderr, "-- rank: %d --\n", rank);
	for (int i = 0; i < dims[0]; i++) {
		for (int j = 0; j < dims[1]; j++) {
			for (int k = 0; k < dims[2]; k++) {
				fprintf(stderr, "%5lf ", data[i][j][k]);
			}
			fprintf(stderr, "%s" ,"\n");
		}
		fprintf(stderr, "%s" ,"\n");
	}
}

double*** allocarr(int dims[]) {
	double*** arr;

	arr 		= (double***)malloc(dims[0]*sizeof(double**));
	arr[0] 		= (double**)malloc(dims[0]*dims[1]*sizeof(double*));
	arr[0][0]   = (double*)malloc(dims[0]*dims[1]*dims[2]*sizeof(double));

	for(int i = 1; i < dims[0]; i++) {
		arr[i] = arr[i-1] + dims[1];
		arr[i][0] = arr[i-1][0] + dims[1]*dims[2];
	}

	for (int i = 0; i < dims[0]; i++) {
		for (int j = 1; j < dims[1]; j++) {
			arr[i][j] = arr[i][j-1] + dims[2];
		}
	}

	bzero(arr[0][0], dims[0]*dims[1]*dims[2]*sizeof(double));
	
	return arr;
}

void freearr(double**** arr) {
	double*** pArr = *arr;
	if (pArr[0][0]) {
		free(pArr[0][0]);
	}
	if (pArr[0]) {
		free(pArr[0]);
	}
	if (pArr) {
		free(pArr);
	}

	*arr = NULL;

}

void initarr(double**** arr, int dims[], int rank) {
	int start = dims[0]*dims[1]*dims[2]*rank;

	for (int i = 0; i < dims[0]; i++) {
		for (int j = 0; j < dims[1]; j++) {
			for (int k = 0; k < dims[2]; k++) {
				(*arr)[i][j][k] = i*dims[1]*dims[2] + j*dims[2]+ k + start;
			}
		}
	}
}

