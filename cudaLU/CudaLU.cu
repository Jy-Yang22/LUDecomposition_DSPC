#include <iostream>
#include <iomanip>
#include <cmath>
#include <omp.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>  
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

using namespace std;
#define TILE 3

//Get input at the argument in the properties
bool inputMatrix(int argc, char* argv[], int& n, int& isPrint)
{
	bool correct = true;

	if (argc < 2)
	{
		cout << "Arguments:X Y" << endl;
		cout << "X : Matrix size [N x N]" << endl;
		cout << "Y = 1: print the input/output matrix if X < 10" << endl;
		correct = false;
	}
	else
	{
		//get matrix size in argument pos 1
		n = atoi(argv[1]);
		if (n <= 0)
		{
			cout << "Matrix size must be larger than 0" << endl;
			correct = false;
		}

		//is print the input/output matrix
		if (argc >= 3)
		{
			isPrint = (atoi(argv[2]) == 1 && n <= 9) ? 1 : 0;
		}
		else
		{
			isPrint = 0;
		}
	}
	return correct;
}

//Initialize the value of matrix a[n x n]
void InitializeMatrix(float**& a, int n)
{
	a = new float* [n];
	a[0] = new float[n * n];

	for (int i = 1; i < n; i++)
	{
		a[i] = a[i - 1] + n;
	}

	for (int i = 0; i < n; i++)
	{
		for (int j = 0; j < n; j++)
		{
			if (i == j)
			{
				a[i][j] = (((float)i + 1) * ((float)i + 1)) / (float)2;
			}

			else
			{
				a[i][j] = (((float)i + 1) + ((float)j + 1)) / (float)2;
			}
		}
	}
}

//Initialize the value of matrix L[n x n] for lower triangular matrix
void InitializeLowerMatrix(float**& L, int n) {
	L = new float* [n];
	L[0] = new float[n * n];

	for (int i = 1; i < n; i++)
	{
		L[i] = L[i - 1] + n;
	}

	for (int j = 0; j < n; j++)
	{
		for (int i = 0; i < n; i++)
		{
			if (i == j)
			{
				L[j][i] = 1;
			}
			else
			{
				L[j][i] = 0;
			}
		}
	}
}

//Compute the LU Decomposition for matrix a[n x n] and L[n x n]
__global__ void ComputeLUDecomposition(float* a, float* L, int n)
{
	//extern __shared__ float pivot;
	//Define the variables
	float pivot, gmax, pmax, temp;
	int  pindmax, gindmax, i, j, k;

	int tx = threadIdx.x; //threadIdx.x to access thread index within block
	int ty = threadIdx.y;
	int Row = blockIdx.x * TILE + tx; // blockIdx.x to access block index within grid
	int Col = blockIdx.y * TILE + ty;

	// Synchronize to make sure the sub-matrices are loaded
	// before starting the computation
	__syncthreads(); //a barrier use to prevent data hazards
	if (Row < n && Col < n)
	{
		//Perform rowwise elimination
		for (k = 0; k < n - 1; k++)
		{
			gmax = 0.0;

			//Find the pivot row among rows k, k+1,...n
			//Each thread works on a number of rows to find the local max value pmax
			//Then update this max local value to the global variable gmax
			{
				pmax = 0.0;
				for (i = k; i < n; i++)
				{
					temp = abs(a[i * n + k]);

					if (temp > pmax)
					{
						pmax = temp;
						pindmax = i;
					}
				}

				if (gmax < pmax)
				{
					gmax = pmax;
					gindmax = pindmax;
				}
			}
			//If matrix is singular set the flag & quit
			if (gmax == 0)
			{
				return;
			}

			//Swap rows if necessary
			if (gindmax == k)
			{
				for (j = k; j < n; j++)
				{
					temp = a[gindmax * n + j];
					a[gindmax * n + j] = a[k * n + j];
					a[k * n + j] = temp;
				}
			}

			//Compute the pivot
			pivot = -1.0 / a[k * n + k];

			//Perform row reductions
			for (i = k + 1; i < n; i++)
			{
				temp = pivot * a[i * n + k];
				L[i * n + k] = ((-1.0) * temp);
				for (j = k; j < n; j++)
				{
					a[i * n + j] = a[i * n + j] + temp * a[k * n + j];
				}
			}
		}
	}
	return;
}

//Print lower triangular matrix	
void PrintMatrix(float** a, int n)
{
	for (int i = 0; i < n; i++)
	{
		cout << "Row " << (i + 1) << ":\t";
		for (int j = 0; j < n; j++)
		{
			printf("%.2f\t", a[i][j]);
		}
		cout << endl;
	}
}

//Print upper triangular matrix	
void PrintMatrixU(float** a, int n)
{
	for (int i = 0; i < n; i++)
	{
		cout << "Row " << (i + 1) << ":\t";
		for (int j = 0; j < n; j++)
		{
			if (j < i) {
				a[i][j] = 0;
			}
			printf("%.2f\t", a[i][j]);
		}
		cout << endl;
	}
}

int main(int argc, char* argv[])
{
	int n = 0, isPrintMatrix = 0;
	float** a;
	float** L;
	float* da, * dl, * du; //device pointers
	double runtime;
	bool correct;

	if (inputMatrix(argc, argv, n, isPrintMatrix) == false)
	{
		return 1;
	}

	cout << "Cuda 1 - gpu matrix " << endl;
	cout << "matrix size is " << n << endl;

	runtime = clock() / (double)CLOCKS_PER_SEC;

	//Initialize the value of matrix A[n x n]
	InitializeMatrix(a, n);
	InitializeLowerMatrix(L, n);

	if (isPrintMatrix == 1)
	{
		cout << "CUDA LU decomposition" << endl;
		cout << "=========================" << endl << endl;
		cout << "Generated : " << n << " x " << n << " Matrix" << endl;
		PrintMatrix(a, n);
	}

	//Declare grid size and block size
	int numblock = n / TILE + ((n % TILE) ? 1 : 0);
	dim3 dimGrid(numblock, numblock); //Dimensions of the grid in blocks
	dim3 dimBlock(TILE, TILE);// Dimensions of the block in threads

	//Allocate memory on device
	//cudaMalloc((void**)&da, n * n * sizeof(float));
	cudaMalloc((void**)&dl, n * n * sizeof(float));
	cudaMalloc((void**)&du, n * n * sizeof(float));
	
	//Copy data to the device
	cudaMemcpy(du, a[0], n * n * sizeof(float), cudaMemcpyHostToDevice);
	cudaMemcpy(dl, L[0], n * n * sizeof(float), cudaMemcpyHostToDevice);

	//Compute the LU decomposition for matrix a[n x n]
	//correct = ComputeLUDecomposition(a, L, n);

	//Do the matrix multiplication on the device (GPU)
	ComputeLUDecomposition << < dimGrid, dimBlock >> > (du, dl, n);

	// wait for the gpu to finish
	cudaDeviceSynchronize();

	//Get results from the device
	cudaMemcpy(L[0], dl, n * n * sizeof(float), cudaMemcpyDeviceToHost);
	cudaMemcpy(a[0], du, n * n * sizeof(float), cudaMemcpyDeviceToHost);

	runtime = (clock() / (double)CLOCKS_PER_SEC) - runtime;


	//The eliminated matrix is as below:
	if (isPrintMatrix == 1)
	{
		cout << "\nLower Triangular Matrix:" << endl;
		PrintMatrix(L, n);
		cout << "\nUpper Triangular Matrix:" << endl;
		PrintMatrixU(a, n);
	}

	//print computing time
	cout << "\n\nLU Decomposition take: " << setiosflags(ios::fixed) << setprecision(8) << runtime << " seconds\n\n";
	cout << "Matrix size  = " << n << endl;

	//cudaFree(da);
	cudaFree(dl);
	cudaFree(du);

	delete[] a[0];
	delete[] L[0];
	delete[] a;
	delete[] L;

	return 0;
}