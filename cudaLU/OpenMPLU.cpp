#include <iostream>
#include <iomanip>
#include <cmath>
#include <omp.h>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>  
using namespace std;

//Get input at the argument in the properties
bool inputMatrix(int argc, char* argv[], int& n, int& numThreads, int& isPrint)
{
	bool correct = true;

	if (argc < 3)
	{
		cout << "Arguments:X Y Z" << endl;
		cout << "X : Matrix size [N x N]" << endl;
		cout << "Y : Number of threads" << endl;
		cout << "Z = 1: print the input/output matrix if X < 10" << endl;
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
		if (argc >= 4)
		{
			isPrint = (atoi(argv[3]) == 1 && n <= 9) ? 1 : 0;
		}
		else
		{
			isPrint = 0;
		}

		//get threads number in argument pos 2
		numThreads = atoi(argv[2]);
		if (numThreads <= 0)
		{
			cout << "Number of threads must be larger than 0" << endl;
			correct = false;
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

	//each thread divide the work 1 time the loop until N time (if didnt set the chunk-size then it will run 1 time) 
#pragma omp parallel for schedule(static) 
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

#pragma omp parallel for schedule(static)
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
bool ComputeLUDecomposition(float** a, float** L, int n)
{
	//Define the variables
	float pivot, gmax, pmax, temp;
	int  pindmax, gindmax, i, j, k;

	omp_lock_t lock;

	omp_init_lock(&lock);

	//Perform rowwise elimination
	for (k = 0; k < n - 1; k++)
	{
		gmax = 0.0;

		//Find the pivot row among rows k, k+1,...n
		//Each thread works on a number of rows to find the local max value pmax
		//Then update this max local value to the global variable gmax
#pragma omp parallel shared(a,gmax,gindmax) firstprivate(n,k) private(pmax,temp,pindmax,i,j)
		{
			pmax = 0.0;
#pragma omp for schedule(dynamic) 
			for (i = k; i < n; i++)
			{
				temp = abs(a[i][k]);

				if (temp > pmax)
				{
					pmax = temp;
					pindmax = i;
				}
			}

			omp_set_lock(&lock);
			if (gmax < pmax)
			{
				gmax = pmax;
				gindmax = pindmax;
			}
			omp_unset_lock(&lock);
		}

#pragma omp critical //identifies this section must be executed by a single thread at a time
		{
			//If matrix is singular set the flag & quit
			if (gmax == 0)
			{
				return false;
			}
		}

		//Swap rows if necessary
		if (gindmax == k)
		{
#pragma omp parallel for shared(a) firstprivate(n,k,gindmax) private(j,temp) schedule(dynamic)
			for (j = k; j < n; j++)
			{
				temp = a[gindmax][j];
				a[gindmax][j] = a[k][j];
				a[k][j] = temp;
			}
		}

		//Compute the pivot
		pivot = -1.0 / a[k][k];

		//Perform row reductions
#pragma omp parallel for shared(a,L) firstprivate(pivot,n,k) private(i,j,temp) schedule(dynamic)
		for (i = k + 1; i < n; i++)
		{
			temp = pivot * a[i][k];
			L[i][k] = ((-1.0) * temp);
			for (j = k; j < n; j++)
			{
				a[i][j] = a[i][j] + temp * a[k][j];
			}
		}
	}
	omp_destroy_lock(&lock);

	return true;
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
	int n = 0, numThreads = 0, isPrintMatrix = 0;
	float** a;
	float** L;
	double runtime;
	bool correct;

	if (inputMatrix(argc, argv, n, numThreads, isPrintMatrix) == false)
	{
		return 1;
	}

	//specify number of threads created in parallel region
	omp_set_num_threads(numThreads);

	runtime = omp_get_wtime();

	//Initialize the value of matrix A[n x n]
	InitializeMatrix(a, n);
	InitializeLowerMatrix(L, n);

	if (isPrintMatrix == 1)
	{
		cout << "OpenMP LU decomposition" << endl;
		cout << "=========================" << endl << endl;
		cout << "Generated : " << n << " x " << n << " Matrix" << endl;
		PrintMatrix(a, n);
	}

	//Compute the LU decomposition for matrix a[n x n]
	correct = ComputeLUDecomposition(a, L, n);

	runtime = omp_get_wtime() - runtime;

	if (correct == true)
	{
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
		cout << "Total threads = " << omp_get_max_threads() << endl << endl;
		cout << "Matrix size  = " << n << endl;
	}
	else
	{
		cout << "The matrix is singular" << endl;
	}

	delete[] a[0];
	delete[] L[0];
	delete[] a;
	delete[] L;
	return 0;
}