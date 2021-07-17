//Project v1
#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

void computeass(int nrighe, int nproc, int* vector)
{
	int div;
	int resto;
	div= nrighe/nproc;
	resto = nrighe%nproc;
	for(int i=0; i<nproc;i++)
	{
		vector[i]=div;
	}
	int j=0;
	while(resto>0)
	{
		vector[j]=vector[j]+1;
		resto--;
		j++;
	}
}

void printm(int** mat, int count)
{
	for(int i=0; i<count;i++)
		{
			for(int j =0 ; j< count; j++)
			{
				if(mat[i][j]==2147483647)
				{
					printf(" i ");
				}
				else
					printf(" %d " , mat[i][j]);
			}
			printf("\n");
		}
}

int ** Readtxt(int * k, char* argv[])
{
	int kk=0;
	FILE *arc;
	
    arc = fopen(argv[1],"r");
    if (arc == NULL)
        exit(1);
	char key;
	int count1=1;
	
	
	for (key = getc(arc); key!= '\n'; key = getc(arc)) //bisogna modificare qui se necessario
		if(key!=' ')
		{
        	//c[kk] = key;
			kk++;
		}else
		{
			kk=0;
			count1++;
		}
	
	count1--; //Erase this line if the matrix is of the form value-space-value-\n
	int **mat = (int **)malloc(count1 * sizeof(int*));
	for(int i = 0; i < count1; i++) mat[i] = (int *)malloc(count1 * sizeof(int));
	fclose(arc);
	
	arc = fopen(argv[1],"r");
	
	for(int i = 0; i< count1; i++)
	{
		for(int j=0; j<count1;j++)
		{
			if(fscanf(arc, "%d", &mat[i][j]) != 1)
			{
				fclose(arc);
				exit(1);
			}
				
		}
	}
	fclose(arc);
	*k = count1;
	return mat;
	
}
void scatter_val(int rank, int k, int * vector, int nofp, MPI_Status status)
{
	if(rank==0)
	{
		MPI_Send(vector,1,MPI_INT,(rank+1)%nofp,0,MPI_COMM_WORLD);
	}
	else if(rank==nofp-1)
	{
		MPI_Recv(vector,1,MPI_INT,(rank-1)%nofp,0,MPI_COMM_WORLD, &status);
	}
	else
	{
		MPI_Recv(vector,1,MPI_INT,(rank-1)%nofp,0,MPI_COMM_WORLD, &status);
		MPI_Send(vector,1,MPI_INT,(rank+1)%nofp,0,MPI_COMM_WORLD);
		
	}
}




void scattermatrix(int rank, int root, int nofp, int** matrix, int **semimatrix, int ndirighe, int* arr, MPI_Status status)
{
	if(rank==0)
	{
		int send=0;
		int* savedarr= (int*)malloc(sizeof(int*)*ndirighe); 

		//count how many send i have to do
		for(int i=0; i<nofp; i++)
			send=send+arr[i];

		for(int i = ndirighe-1; i>=0 ; i--)
		{
			for(int j = 0 ; j<ndirighe;j++)
			{
				savedarr[j]=matrix[i][j];
			}
			if(send>arr[0])
			{
				MPI_Send(savedarr,ndirighe,MPI_INT,(rank+1)%nofp,0,MPI_COMM_WORLD);
				send--;
			}
		}
		free(savedarr);
	}
	else
	{
		int* savedarr= (int*)malloc(sizeof(int*)*ndirighe);  
		int send=0;
		//counts how many rows I need to propagate
		for(int i=nofp-1; i>rank; i--)
		{
			send=send+arr[i];
		}
		//propagate
		for(int i=0; i<send;i++)
		{
			MPI_Recv(savedarr,ndirighe,MPI_INT,(rank-1)%nofp,0,MPI_COMM_WORLD, &status);
			MPI_Send(savedarr, ndirighe ,MPI_INT, (rank+1)%nofp,0, MPI_COMM_WORLD);
		}
		//save my rows
		for(int i=arr[rank]-1; i>=0;i--)
		{
			MPI_Recv(semimatrix[i],ndirighe,MPI_INT,(rank-1)%nofp,0,MPI_COMM_WORLD, &status);
		}
		free(savedarr);
	}
	
}
void gathermatrix(int rank, int root, int nofp, int** matrix, int **semimatrix, int ndirighe, int* arr, MPI_Status status)
{
	if(rank==0)
	{
		//Build up my submatrix and receive the rest of the matrix
		for(int i=0;i<ndirighe;i++)
		{
			if(i<arr[rank])
			{
				matrix[i]=semimatrix[i];
			}
			else
			{
				MPI_Recv(matrix[i],ndirighe,MPI_INT,(rank-1)%nofp,0,MPI_COMM_WORLD, &status);
				
			}
		}
	}
	else
	{
		int* savedarr=(int*)malloc(sizeof(int*)*ndirighe); 
		int send=0;
		for(int i=1; i<rank; i++)
		{
			send=send+arr[i];//count how many rows i need to propagate
		}
		//propagate all the rows i have to
		for(int i=0;i<send;i++)
		{
			MPI_Recv(savedarr,ndirighe,MPI_INT,(rank-1)%nofp,0,MPI_COMM_WORLD, &status);
			MPI_Send(savedarr,ndirighe,MPI_INT,(rank+1)%nofp,0, MPI_COMM_WORLD);
			
		}
		//finally send my rows
		for(int i=0;i<arr[rank];i++)
		{
			MPI_Send(semimatrix[i], ndirighe ,MPI_INT, (rank+1)%nofp,0, MPI_COMM_WORLD);
		}
		free(savedarr);
	}
}



void floydplusrotation(int ndirighe, int* arr, int nofp, int rank, int **mat, int **semimatrix, int** matfinale, MPI_Status status, int* backarr)
{
	if(rank==0)
	{
		
		int dummy=2147483647;
		int minore=2147483647;
		int ndr=0;
		for(int i=0; i<ndirighe; i++)
		{
			//sending column x column
			for(int mm=0; mm<ndirighe;mm++)
			{
				backarr[mm]=mat[mm][ndr];
			}
			ndr++;			
			MPI_Send(backarr,ndirighe,MPI_INT,(rank+1)%nofp,0,MPI_COMM_WORLD);
			for(int j=0; j<arr[rank];j++)
			{
				dummy=2147483647;
				minore=2147483647;
				//do the computation
				for(int z=0; z<ndirighe;z++)
				{
					if(backarr[z] == 2147483647 || semimatrix[j][z] == 2147483647)
					{
						dummy=2147483647;
					}
					else
					{
						dummy=backarr[z] + semimatrix[j][z];
					}
					if(dummy<minore)
						minore=dummy;	
				}
				//saving the result
				matfinale[j][i]=minore;
				
			}
		}
	}
	else if(rank==nofp-1)
	{
		int dummy=2147483647;
		int minore=2147483647;
		for(int i=0; i<ndirighe; i++)
		{
			MPI_Recv(backarr,ndirighe,MPI_INT,(rank-1)%nofp,0,MPI_COMM_WORLD,&status);
			for(int j=0; j<arr[rank];j++)
			{
				dummy=2147483647;
				minore=2147483647;
				for(int z=0; z<ndirighe;z++)
				{
					if(backarr[z] == 2147483647 || semimatrix[j][z] == 2147483647)
					{
						dummy=2147483647;
					}
					else
					{
						dummy=backarr[z] + semimatrix[j][z];
					}
					if(dummy<minore)
						minore=dummy;	
				}
				matfinale[j][i]=minore;
			}
		}
		
	}
	else
	{
		int dummy=2147483647;
		int minore=2147483647;
		for(int i=0; i<ndirighe; i++)
		{
			MPI_Recv(backarr,ndirighe,MPI_INT,(rank-1)%nofp,0,MPI_COMM_WORLD,&status);
			MPI_Send(backarr,ndirighe,MPI_INT,(rank+1)%nofp,0,MPI_COMM_WORLD);
			for(int j=0; j<arr[rank];j++)
			{
				dummy=2147483647;
				minore=2147483647;
				for(int z=0; z<ndirighe;z++)
				{
					if(backarr[z] == 2147483647 || semimatrix[j][z] == 2147483647)
					{
						dummy=2147483647;
					}
					else
					{
						dummy=backarr[z] + semimatrix[j][z];
						
					}
					if(dummy<minore)
						minore=dummy;	
				}
				
				matfinale[j][i]=minore;
			}
		}
	}
}
void trasf_graf_ad(int **mat,int count)
{
	for(int i =0; i<count;i++)
	{
		for(int j=0;j<count;j++)
		{
			if(mat[i][j] == 0 && i!=j)
			{
				mat[i][j]=2147483647; //inf biggest int
			}
			if(i==j)
				mat[i][j]=0;
		}
	}
}
int main(int argc, char **argv)
{
	int nofp;
  	int rank;
  	MPI_Init(&argc,&argv);
  	MPI_Status status;
  	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
  	MPI_Comm_size(MPI_COMM_WORLD,&nofp);
	int arr[nofp];
	int ndirighe;
	int k;
	int **matrice;
	int **matfinale;
	
	if(rank==0)
	{
		//Reading the matrix
		matrice=Readtxt(&ndirighe, argv);
		//trasf adj graph
		trasf_graf_ad(matrice, ndirighe);
		//compute vector of rows assig
		computeass(ndirighe,nofp,arr);
		//Scatter n of rows
		scatter_val(rank,0,&ndirighe,nofp,status);
		//build up partial matrix
		int **parmat = (int **)malloc(arr[rank] * sizeof(int*));
		for(int i = 0; i < arr[rank]; i++) parmat[i] = (int *)malloc(ndirighe * sizeof(int));
		//buildup final matrix
		int **matfinale = (int **)malloc(arr[rank] * sizeof(int*));
		for(int i = 0; i < arr[rank]; i++) matfinale[i] = (int *)malloc(ndirighe * sizeof(int));
		//Build up backup array
		int* backarr = (int*) malloc(sizeof(int*)*ndirighe);
		//Do Floyd-Warshall	
		for(int zz=1;zz<=ndirighe;zz=zz*2)
		{
			scattermatrix(rank,0,nofp,matrice,parmat,ndirighe,arr,status);
			//Build sub matrix of p0
			for(int i=0;i<arr[rank];i++)
			{
				for(int j=0;j<ndirighe;j++)
				{
					parmat[i][j]=matrice[i][j];

				}
			}
			floydplusrotation(ndirighe, arr, nofp, rank, matrice, parmat, matfinale, status,backarr);
			gathermatrix(rank,0,nofp,matrice,matfinale,ndirighe,arr,status);	
		}
		printm(matrice,ndirighe);
	}
	else
	{
		scatter_val(rank,0,&ndirighe,nofp,status);
		computeass(ndirighe,nofp,arr);
		//build up my submatrix
		int **parmat = (int **)malloc(arr[rank] * sizeof(int*));
		for(int i = 0; i < arr[rank]; i++) parmat[i] = (int *)malloc(ndirighe * sizeof(int));
		//build up my finalmatrix
		int **matfinale = (int **)malloc(arr[rank] * sizeof(int*));
		for(int i = 0; i < arr[rank]; i++) matfinale[i] = (int *)malloc(ndirighe * sizeof(int));
		//build up backup array for computation
		int* backarr = (int*) malloc(sizeof(int*)*ndirighe);
		//do the algorithm
		for(int zz=1;zz<=ndirighe;zz=zz*2)
		{
			scattermatrix(rank,0,nofp,matrice,parmat,ndirighe,arr,status);
			floydplusrotation(ndirighe, arr, nofp, rank, matrice, parmat, matfinale, status,backarr);
			gathermatrix(rank,0,nofp,matrice,matfinale,ndirighe,arr,status);
			
		}
	}
	MPI_Finalize();
}