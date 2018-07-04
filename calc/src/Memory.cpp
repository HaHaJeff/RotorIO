/*---------------------------------------------------------------------------*\


\*---------------------------------------------------------------------------*/
#ifndef Memory_CPP
#define Memory_CPP

#include "Memory.h"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

template<class T>
bool CMemory<T>::Malloc(T* &pArr,int size1)
{
	pArr = (T *)malloc(sizeof(T ) * size1 );
 	if(pArr == NULL)
 	{	
		cout<<" memory allocation error"<<endl;
		return -1;	
 	}
	memset(pArr, 0, size1*sizeof(T));
	return true;
}

template<class T>
bool CMemory<T>::Malloc(T** &pArr,int size1,int size2)
{
	pArr = (T **)malloc(sizeof(T *) *size1); 
	pArr[0] = (T *)malloc(sizeof(T ) * size1 * size2);

 	if(pArr == NULL)
 	{	
		cout<<" memory allocation error"<<endl;
		return -1;	
 	}

	for(int n=1; n<size1; n++)  
	{  
		pArr[n] = pArr[n-1] + size2; 
	}

	
	memset(pArr[0], 0, size1*size2*sizeof(T));

}

template<class T>
bool CMemory<T>::Malloc(T*** &pArr,int size1,int size2,int size3)
{
	pArr = (T ***)malloc(sizeof(T **) *size1); 
	pArr[0] = (T **)malloc(sizeof(T *) * size1 * size2);
	pArr[0][0]= (T *)malloc(sizeof(T ) * size1 * size2* size3);

 	if(pArr == NULL)
 	{	
		cout<<" memory allocation error"<<endl;
		return -1;	
 	}

	for(int n=1; n<size1; n++)  
	{  
		pArr[n] = pArr[n-1] + size2;
		pArr[n][0]= pArr[n-1][0]+size2*size3; 
	}

	for(int i=0; i<size1; i++) 
	for(int j=1; j<size2; j++) 
	{  
		pArr[i][j] = pArr[i][j-1] + size3;  
	}	
	memset(pArr[0][0], 0, size1*size2*size3*sizeof(T));

	return true;
}

template<class T>
bool CMemory<T>::Malloc(T**** &pArr,int size1,int size2,int size3,int size4)
{
	pArr = (T ****)malloc(sizeof(T ***) *size1); 
	pArr[0] = (T ***)malloc(sizeof(T **) * size1 * size2);
	pArr[0][0]= (T **)malloc(sizeof(T *) * size1 * size2* size3);
	pArr[0][0][0]= (T *)malloc(sizeof(T ) * size1 * size2* size3*size4);

	for(int n=1; n<size1; n++)  
	{  
		pArr[n] = pArr[n-1] + size2;
		pArr[n][0]= pArr[n-1][0]+size2*size3; 
		pArr[n][0][0]= pArr[n-1][0][0]+size2*size3*size4;
	}

	for(int n1=0; n1<size1; n1++) 
	for(int n2=1; n2<size2; n2++) 
	{  
		pArr[n1][n2] = pArr[n1][n2-1] + size3;
		pArr[n1][n2][0]= pArr[n1][n2-1][0]+size3*size4;  
	}

	for(int i=0; i<size1; i++) 
	for(int j=0; j<size2; j++) 
	for(int k=1; k<size3; k++)
	{  
		pArr[i][j][k] = pArr[i][j][k-1] + size4;  
	}

    if(pArr == 0)
    {
        return -1;
    }
	
    memset(pArr[0][0][0], 0, size1*size2*size3*size4*sizeof(T));
    return true;

}

template<class T>
bool CMemory<T>::Free(T* &pArr)
{	if(pArr)
	free(pArr);
	pArr=NULL;
	return true;
}

template<class T>
bool CMemory<T>::Free(T** &pArr)
{	
	if(pArr[0]) 
	free(pArr[0]);
	if(pArr)
	free(pArr);

	pArr=NULL;
	return true;
}

template<class T>
bool CMemory<T>::Free(T*** &pArr)
{
	if(pArr[0][0])
	free(pArr[0][0]);
	if(pArr[0])
	free(pArr[0]);
	if(pArr)
	free(pArr); 
	

	pArr=NULL;
	return true;	
}

template<class T>
bool CMemory<T>::Free(T**** &pArr)
{

	if(pArr[0][0][0])
	free(pArr[0][0][0]);	
	if(pArr[0][0])
	free(pArr[0][0]);
	if(pArr[0])
	free(pArr[0]);
	if(pArr)
	free(pArr); 
		
	pArr=NULL;
	return true;
}
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


#endif
