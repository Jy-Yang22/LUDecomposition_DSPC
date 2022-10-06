//#include <iostream>
//#include <iomanip>
//#include <cmath>
//#include <stdio.h>
//#include <stdlib.h> 
//#include <mpi.h>    
//
//using namespace std;
//
////Get user input of matrix dimension and printing option
//bool GetUserInput(int argc, char* argv[], int& n, int& isPrint, int numProcesses, int rank)
//{
//	bool isOK = true;
//
//	if (argc < 2)
//	{
//		if (rank == 0)
//		{
//			cout << "Arguments:<X> [<Y>]" << endl;
//			cout << "X : Matrix size [X x X]" << endl;
//			cout << "Y = 1: print the input/output matrix if X < 10" << endl;
//		}
//		isOK = false;
//	}
//	else
//	{
//		//get matrix size
//		n = atoi(argv[1]);
//		if (n <= 0)
//		{
//			if (rank == 0)
//			{
//				cout << "Matrix size must be larger than 0" << endl;
//			}
//			isOK = false;
//		}
//		//check if matrix size is multiple of processes
//		if ((n % numProcesses) != 0)
//		{
//			if (rank == 0)
//			{
//				cout << "Matrix size must be multiple of the number of processes" << endl;
//			}
//			isOK = false;
//		}
//
//		//is print the input/output matrix
//		if (argc >= 3)
//		{
//			isPrint = (atoi(argv[2]) == 1 && n <= 9) ? 1 : 0;
//		}
//		else
//		{
//			isPrint = 0;
//		}
//	}
//	return isOK;
//}
//
////Initialize the value of matrix a[n x n]
//void InitializeMatrix(float**& a, int n)
//{
//	a = new float* [n];
//	a[0] = new float[n * n];
//
//	for (int i = 1; i < n; i++)
//	{
//		a[i] = a[i - 1] + n;
//	}
//
//	for (int j = 0; j < n; j++)
//	{
//		for (int i = 0; i < n; i++)
//		{
//			if (i == j)
//			{
//				a[j][i] = (((float)i + 1) * ((float)i + 1)) / (float)2;
//			}
//			else
//			{
//				a[j][i] = (((float)i + 1) + ((float)j + 1)) / (float)2;
//			}
//		}
//	}
//}
//
//void InitializeMatrixL(float**& L, int n) {
//	L = new float* [n];
//	L[0] = new float[n * n];
//
//	for (int i = 1; i < n; i++)
//	{
//		L[i] = L[i - 1] + n;
//	}
//
//	for (int j = 0; j < n; j++)
//	{
//		for (int i = 0; i < n; i++)
//		{
//			if (i == j)
//			{
//				L[j][i] = 1;
//			}
//			else
//			{
//				L[j][i] = 0;
//			}
//		}
//	}
//}
//
////Compute the lu decomposition for matrix a[n x n]
//bool ComputeLUDecomposition(float** a, float** L, int n, int numProcesses, int rank)
//{
//	float pivot, max, temp;
//	int gindmax, i, j, k, lk, master;
//	int nCols = n / numProcesses;
//	float* tmp = new float[n];
//	float** b = new float* [n]; //create local matrix
//	b[0] = new float[n * nCols];
//	MPI_Status status;
//
//	for (int j = 1; j < n; j++)
//	{
//		b[j] = b[j - 1] + n;
//	}
//
//	//process 0 send the data to compute nodes
//	if (rank == 0)
//	{
//		//copy data part of each proccess from matrix a to local matrix b
//		for (k = numProcesses - 1; k >= 0; k--)
//		{
//			for (j = 0; j < nCols; j++)
//			{
//				for (i = 0; i < n; i++)
//				{
//					b[j][i] = a[j * numProcesses + k][i];
//				}
//			}
//
//			//send it to work proceess
//			if (k != 0)
//			{
//				MPI_Send(b[0], n * nCols, MPI_FLOAT, k, 0, MPI_COMM_WORLD);
//			}
//		}
//	}
//	else
//	{
//		MPI_Recv(b[0], n * nCols, MPI_FLOAT, 0, 0, MPI_COMM_WORLD, &status);
//	}
//
//	//Perform rowwise elimination
//	for (k = 0; k < n - 1; k++)
//	{
//		max = 0.0;
//		gindmax = k;
//
//		//master process ID
//		master = k % numProcesses;
//		//local k
//		lk = k / numProcesses;
//
//		//master process find the pivot row
//		//Then broadcast it to all other processes
//		if (rank == master)
//		{
//			//Find the pivot row
//			for (i = k; i < n; i++)
//			{
//				temp = abs(b[lk][i]);
//				if (temp > max)
//				{
//					max = temp;
//					gindmax = i;
//				}
//			}
//		}
//		MPI_Bcast(&max, 1, MPI_FLOAT, master, MPI_COMM_WORLD);
//		MPI_Bcast(&gindmax, 1, MPI_INT, master, MPI_COMM_WORLD);
//
//		//If matrix is singular set the flag & quit
//		if (max == 0) return false;
//		// Swap rows if necessary
//		if (gindmax == k)
//		{
//			//j = lk;
//			if (rank < master)
//			{
//				j++;
//			}
//
//			for (j = lk; j < nCols; j++)
//			{
//				temp = b[j][gindmax];
//				b[j][gindmax] = b[j][k];
//				b[j][k] = temp;
//			}
//		}
//		//Master 
//		if (rank == master)
//		{
//			pivot = -1.0 / b[lk][k];
//			for (i = k + 1; i < n; i++)
//			{
//				tmp[i] = pivot * b[lk][i];
//
//			}
//		}
//		MPI_Bcast(tmp + k + 1, n - k - 1, MPI_FLOAT, master, MPI_COMM_WORLD);
//		// after tmp is broadcast to all processes, add it to L only on pid 0
//		if (rank == 0)
//		{
//			for (i = k + 1; i < n; i++)
//			{
//				L[k][i] = ((-1.0) * tmp[i]);
//			}
//		}
//
//		//Perform row reductions
//		//j = lk;
//		if (rank < master)
//		{
//			j++;
//		}
//
//		for (j = lk; j < nCols; j++)
//		{
//			for (i = k + 1; i < n; i++)
//			{
//				b[j][i] = b[j][i] + tmp[i] * b[j][k];
//			}
//		}
//	}
//
//	//process 0 collects results from the worker processes
//	if (rank == 0)
//	{
//		//copy data part of each proccess from matrix a to local matrix b
//		for (k = 0; k < numProcesses; k++)
//		{
//			//receive data from worker proceess
//			if (k != 0)
//			{
//				MPI_Recv(b[0], n * nCols, MPI_FLOAT, k, 0, MPI_COMM_WORLD, &status);
//			}
//
//			for (j = 0; j < nCols; j++)
//			{
//				for (i = 0; i < n; i++)
//				{
//					a[j * numProcesses + k][i] = b[j][i];
//				}
//			}
//		}
//	}
//	else
//	{
//		//worker processes send the data to process 0
//		MPI_Send(b[0], n * nCols, MPI_FLOAT, 0, 0, MPI_COMM_WORLD);
//	}
//
//	delete[] b[0];
//	delete[] b;
//	delete[] tmp;
//	return true;
//}
//
////Print matrix	
//void PrintMatrix(float** a, int n)
//{
//	for (int i = 0; i < n; i++)
//	{
//		cout << "Row " << (i + 1) << ":\t";
//		for (int j = 0; j < n; j++)
//		{
//			printf("%.2f\t", a[j][i]);
//		}
//		cout << endl;
//	}
//}
//
//void PrintMatrixU(float** a, int n)
//{
//	for (int i = 0; i < n; i++)
//	{
//		cout << "Row " << (i + 1) << ":\t";
//		for (int j = 0; j < n; j++)
//		{
//			if (j < i)
//			{
//				a[j][i] = 0;
//			}
//			printf("%.2f\t", a[j][i]);
//		}
//		cout << endl;
//	}
//}
//
//// Main Program
//int main(int argc, char* argv[])
//{
//	float** a{};
//	float** L{};
//	int	n, isPrintMatrix, size, rank;
//	bool isOK;
//	double begin = 0.0, end = 0.0;
//
//	MPI_Init(&argc, &argv);
//	MPI_Comm_size(MPI_COMM_WORLD, &size);
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//
//	//Master process get program input
//	if (GetUserInput(argc, argv, n, isPrintMatrix, size, rank) == false)
//	{
//		MPI_Finalize();
//		return 1;
//	}
//
//	//Master part
//	if (rank == 0)
//	{
//		//Get start time
//		begin = MPI_Wtime();
//
//		//Initialize the value of matrix A[N][N]
//		InitializeMatrix(a, n);
//		InitializeMatrixL(L, n);
//
//		//Prints the input maxtrix if needed
//		if (isPrintMatrix == 1)
//		{
//			cout << "MPI LU decomposition" << endl;
//			cout << "=========================" << endl << endl;
//			cout << "Generated : " << n << " x " << n << " Matrix" << endl;
//			PrintMatrix(a, n);
//			cout << endl;
//			//cout<< "L initialized to:" << endl;
//			//PrintMatrix(L,n);
//		}
//	}
//
//	//Compute the Gaussian Elimination for matrix a[n x n]
//	isOK = ComputeLUDecomposition(a, L, n, size, rank);
//
//	//Master process gets end time and print results
//	if (rank == 0)
//	{
//		end = MPI_Wtime();
//
//		if (isOK == true)
//		{
//			//Print result matrix
//			if (isPrintMatrix == 1)
//			{
//				cout << "Lower matrix:" << endl;
//				PrintMatrix(L, n);
//				cout << endl;
//				cout << "Upper matrix:" << endl;
//				PrintMatrixU(a, n);
//				cout << endl;
//			}
//			cout << "LU Decomposition runs in " << setiosflags(ios::fixed) << setprecision(4) << end - begin << " seconds \n";
//			cout << "Matrix size  = " << n << endl;
//			cout << "Thread Number = " << size << endl;
//		}
//		else
//		{
//			cout << "The matrix is singular" << endl;
//		}
//
//		//All process delete matrix
//		delete[] a[0];
//		delete[] a;
//		delete[] L[0];
//		delete[] L;
//	}
//	MPI_Barrier(MPI_COMM_WORLD);
//	MPI_Finalize();
//	return 0;
//}