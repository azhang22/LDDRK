// Ultrasound_test.cpp : Defines the entry point for the console application.
//

//+++++++++++++++++++++++++++filename: main.cpp +++++++++++++++++++++++++++++++//
#include <stdlib.h>
#include <string>
#include <math.h>
#include "mpi.h"
#include "airpore.h"

void main(int argc, char *argv[])
{
	int mpi_rank, mpi_size;
	MPI_Init(&argc, &argv);
	//mpi_rank=MPI::COMM_WORLD.Get_rank();
	MPI_Comm_size(MPI_COMM_WORLD, &mpi_size);
	// mpi_size=MPI::COMM_WORLD.Get_size();
	MPI_Comm_rank(MPI_COMM_WORLD, &mpi_rank);
	printf("size = %d, rank = %d\n", mpi_size, mpi_rank);
	char infile[100] = " ";
	strcpy_s(infile, "input.txt");
	airpore *AirPore;
	printf("%d,%d\n", mpi_rank, 1);
	AirPore = new airpore(infile);
	printf("%d,%d\n", mpi_rank, 2);
	//	goto loop;
	AirPore->get_output(mpi_rank, mpi_size);
	printf("%d,%d\n", mpi_rank, 3);
	//	loop:
	if (mpi_rank == 0)
	{
		AirPore->get_FFT_y1(2);// for the long-distance sound propagatio
		AirPore->get_FFT_y11(2);
		AirPore->get_FFT_y12(2);
		AirPore->get_FFT_y13(2);
		AirPore->get_FFT_y14(2);
		AirPore->get_FFT_y15(2);
		AirPore->get_FFT_y16(2);
	}

	delete AirPore;
	MPI_Finalize();
}