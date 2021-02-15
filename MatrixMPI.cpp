#include <stdio.h>
#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>

#define NUM_ROWS_A 5//3
#define NUM_COLUMNS_A 5//3
#define NUM_ROWS_B 5//3
#define NUM_COLUMNS_B 5//3
#define MASTER_TO_SLAVE_TAG 1
#define SLAVE_TO_MASTER_TAG 4
void makeAB();
void printArray();
int rank;
int size;
int i, j, k;
double mat_a[NUM_ROWS_A][NUM_COLUMNS_A];
double mat_b[NUM_ROWS_B][NUM_COLUMNS_B];
double mat_result[NUM_ROWS_A][NUM_COLUMNS_B];
double start_time;
double end_time;
int low_bound; //low bound of the number of rows of [A] allocated to a slave
int upper_bound; //upper bound of the number of rows of [A] allocated to a slave
int portion; //portion of the number of rows of [A] allocated to a slave
MPI_Status status;
MPI_Request request; //capture request of a MPI_Isend


void inmulltesteParalel();
MPI_Status mpiStatus;
int* nrOfDimensions;
int nrOfRealDimensions;
int* period;
int reorder = 1;
int* coord;
void creeazaTopologie(int[], int);
MPI_Comm comm;

int main(int argc, char* argv[])
{
	MPI_Init(&argc, &argv);
	//MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	MPI_Comm_size(MPI_COMM_WORLD, &size);

	nrOfRealDimensions = atoi(argv[1]);

	//nrOfDimensions = new int[nrOfRealDimensions];
	nrOfDimensions = (int*)malloc(nrOfRealDimensions * sizeof(*nrOfDimensions));
	/*
	nrOfDimensions[0] = atoi(argv[2]);
	nrOfDimensions[1] = atoi(argv[3]);
	nrOfDimensions[2] = atoi(argv[4]);
	*/

	creeazaTopologie(nrOfDimensions, nrOfRealDimensions);

	inmulltesteParalel();


	MPI_Finalize(); //finalize MPI operations
	return 0;
}

void creeazaTopologie(int nrOfDimensions[], int nrOfRealDimensions)
{
	//coord = new int[nrOfRealDimensions];
	coord = (int*)malloc(nrOfRealDimensions * sizeof(*coord));//new int[nrOfRealDimensions];


	/*
	if (size != totalNrOfProcessors)
	{
	printf("Please run with %d processes.\n", totalNrOfProcessors);fflush(stdout);
	MPI_Abort(MPI_COMM_WORLD, 1);
	}
	*/
	//period = new int[nrOfRealDimensions];
	period = (int*)malloc(nrOfRealDimensions * sizeof(*period));//new int[nrOfRealDimensions];

	if (nrOfRealDimensions == 2)
	{
		nrOfDimensions[0] = 2;
		nrOfDimensions[1] = size / 2;
		period[0] = 1; period[1] = 1;

		coord[0] = coord[1] = 0;
	}
	else
	{
		for (int i = 0; i < nrOfRealDimensions; i++)
		{
			nrOfDimensions[i] = (int)pow((double)size, (double)1 / nrOfRealDimensions); //dim_size[1] = q;
			period[i] = 1; //periods[1]= 1;
		}

		coord[0] = coord[1] = coord[2] = 0;
	}

	int totalNrOfProcessors = nrOfDimensions[0] * nrOfDimensions[1];

	MPI_Cart_create(MPI_COMM_WORLD, nrOfRealDimensions, nrOfDimensions, period, reorder, &comm);

	MPI_Comm_rank(comm, &rank);

	MPI_Cart_coords(comm, rank, nrOfRealDimensions, coord);

	switch (nrOfRealDimensions)
	{
	case 1: printf("\nProcesul %d - coord = (%d)\n", rank, coord[0]);
		break;
	case 2: printf("\nProcesul %d - coord = (%d %d)\n", rank, coord[0], coord[1]);
		break;
	case 3: printf("\nProcesul %d - coord = (%d %d %d)\n", rank, coord[0], coord[1], coord[2]);
		break;
	}
	/*
	if(rank!=0)
	{
		int localCounter = 0;
		MPI_Recv(&localCounter, 1, MPI_INT, (rank-1), 0, comm, &mpiStatus);
		printf("\n\tEu, procesul %d, am primit %d, ", rank, localCounter);
		localCounter++;
		printf("am trimis %d", localCounter); fflush(stdout);
		MPI_Send(&localCounter, 1, MPI_INT, ((rank+1)%size), 0, comm);
	}

	if(rank == 0)
	{
		int localCounter = 1;

		//MPI_Send(void *buf, int count, MPI_Datatype datatype, int dest, int tag, MPI_Comm comm)
		MPI_Send(&localCounter, 1, MPI_INT, (rank+1), 0, comm);

		//MPI_Recv(void *buf, int count, MPI_Datatype datatype, int source, int tag,MPI_Comm comm, MPI_Status *status)
		MPI_Recv(&localCounter, 1, MPI_INT, (size-1), 0, comm, &mpiStatus);


		printf("\n\n\tEu, procesul %d, am primit %d ", rank, localCounter);
		printf("\n\n\tCounter= %d\n", localCounter);
		*/

		/**
		printf("\tTopologia este:\n");fflush(stdout);
		int k = 0;
		for(int i = 0; i<nrOfDimensions[0]; i++)
		{
		printf("\t"); fflush(stdout);
		for(int j = 0; j< nrOfDimensions[1]; j++)
		{
		MPI_Cart_coords(comm, k, nrOfRealDimensions, coord);
		printf("(%d,%d)\t", coord[0], coord[1]);fflush(stdout);
		k++;
		}

		printf("\n");
		}
		**/
		//}
}//.end creeazaTopologie()

void inmulltesteParalel()
{
	/* master initializes work*/
	if (rank == 0) {
		makeAB();
		start_time = MPI_Wtime();
		for (i = 1; i < size; i++) {//for each slave other than the master
			portion = (NUM_ROWS_A / (size - 1)); // calculate portion without master
			low_bound = (i - 1) * portion;
			if (((i + 1) == size) && ((NUM_ROWS_A % (size - 1)) != 0)) {//if rows of [A] cannot be equally divided among slaves
				upper_bound = NUM_ROWS_A; //last slave gets all the remaining rows
			}
			else {
				upper_bound = low_bound + portion; //rows of [A] are equally divisable among slaves
			}
			//send the low bound first without blocking, to the intended slave
			MPI_Isend(&low_bound, 1, MPI_INT, i, MASTER_TO_SLAVE_TAG, comm, &request);
			//next send the upper bound without blocking, to the intended slave
			MPI_Isend(&upper_bound, 1, MPI_INT, i, MASTER_TO_SLAVE_TAG + 1, comm, &request);
			//finally send the allocated row portion of [A] without blocking, to the intended slave
			MPI_Isend(&mat_a[low_bound][0], (upper_bound - low_bound) * NUM_COLUMNS_A, MPI_DOUBLE, i, MASTER_TO_SLAVE_TAG + 2, comm, &request);
		}
	}
	//broadcast [B] to all the slaves
	MPI_Bcast(&mat_b, NUM_ROWS_B * NUM_COLUMNS_B, MPI_DOUBLE, 0, comm);
	/* work done by slaves*/
	if (rank > 0) {
		//receive low bound from the master
		MPI_Recv(&low_bound, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG, comm, &status);
		//next receive upper bound from the master
		MPI_Recv(&upper_bound, 1, MPI_INT, 0, MASTER_TO_SLAVE_TAG + 1, comm, &status);
		//finally receive row portion of [A] to be processed from the master
		MPI_Recv(&mat_a[low_bound][0], (upper_bound - low_bound) * NUM_COLUMNS_A, MPI_DOUBLE, 0, MASTER_TO_SLAVE_TAG + 2, comm, &status);
		for (i = low_bound; i < upper_bound; i++) {//iterate through a given set of rows of [A]
			for (j = 0; j < NUM_COLUMNS_B; j++) {//iterate through columns of [B]
				for (k = 0; k < NUM_ROWS_B; k++) {//iterate through rows of [B]
					mat_result[i][j] += (mat_a[i][k] * mat_b[k][j]);
				}
			}
		}
		//send back the low bound first without blocking, to the master
		MPI_Isend(&low_bound, 1, MPI_INT, 0, SLAVE_TO_MASTER_TAG, comm, &request);
		//send the upper bound next without blocking, to the master
		MPI_Isend(&upper_bound, 1, MPI_INT, 0, SLAVE_TO_MASTER_TAG + 1, comm, &request);
		//finally send the processed portion of data without blocking, to the master
		MPI_Isend(&mat_result[low_bound][0], (upper_bound - low_bound) * NUM_COLUMNS_B, MPI_DOUBLE, 0, SLAVE_TO_MASTER_TAG + 2, comm, &request);
	}
	/* master gathers processed work*/
	if (rank == 0) {
		for (i = 1; i < size; i++) {// untill all slaves have handed back the processed data
			//receive low bound from a slave
			MPI_Recv(&low_bound, 1, MPI_INT, i, SLAVE_TO_MASTER_TAG, comm, &status);
			//receive upper bound from a slave
			MPI_Recv(&upper_bound, 1, MPI_INT, i, SLAVE_TO_MASTER_TAG + 1, comm, &status);
			//receive processed data from a slave
			MPI_Recv(&mat_result[low_bound][0], (upper_bound - low_bound) * NUM_COLUMNS_B, MPI_DOUBLE, i, SLAVE_TO_MASTER_TAG + 2, comm, &status);
		}
		end_time = MPI_Wtime();
		printf("\nRunning Time = %f\n\n", end_time - start_time);
		printArray();
	}
}

void makeAB()
{
	for (i = 0; i < NUM_ROWS_A; i++) {
		for (j = 0; j < NUM_COLUMNS_A; j++) {
			mat_a[i][j] = rand() % 10;//i + j;
		}
	}
	for (i = 0; i < NUM_ROWS_B; i++) {
		for (j = 0; j < NUM_COLUMNS_B; j++) {
			mat_b[i][j] = rand() % 10;//i*j;
		}
	}
}
void printArray()
{
	for (i = 0; i < NUM_ROWS_A; i++) {
		printf("\n");
		for (j = 0; j < NUM_COLUMNS_A; j++)
			printf("%8.2f ", mat_a[i][j]);
	}
	printf("\n\n\n");
	for (i = 0; i < NUM_ROWS_B; i++) {
		printf("\n");
		for (j = 0; j < NUM_COLUMNS_B; j++)
			printf("%8.2f ", mat_b[i][j]);
	}
	printf("\n\n\n");
	for (i = 0; i < NUM_ROWS_A; i++) {
		printf("\n");
		for (j = 0; j < NUM_COLUMNS_B; j++)
			printf("%8.2f ", mat_result[i][j]);
	}
	printf("\n\n");
}

