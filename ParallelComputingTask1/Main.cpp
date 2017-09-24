#include "stdlib.h"
#include <random>
#include <iostream>
#include <omp.h>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <complex>

using namespace std;



//int** matrix_multipl_par(int** matr1, int** matr2, int n, int threads_num = 1, int mode = 1) {
//	int** matr = zero_matrix(n);
//	if (mode == 4)
//		threads_num = int(threads_num / 2);
//	omp_set_num_threads(threads_num);
//	omp_set_nested(true);
//	double time1, time2;
//	time1 = omp_get_wtime();
//	int i = 0, j = 0, k = 0;
//#pragma omp parallel for if(mode == 2 || mode == 4) private(j, k)
//	for (i = 0; i < n; i++) {
//#pragma omp parallel for if(mode == 3 || mode == 4) private(k)
//		for (j = 0; j < n; j++) {
//			for (k = 0; k < n; k++) {
//				matr[i][j] += matr1[i][k] * matr2[k][j];
//			}
//		}
//	}
//	time2 = omp_get_wtime();
//	res_time = (time2 - time1) / (double)CLOCKS_PER_SEC;
//	switch (mode)
//	{
//	case 1:
//		cout << "Time standart " << threads_num << " : " << time2 - time1 << '\n';
//		break;
//	case 2:
//		cout << "Time paral i " << threads_num << " : " << time2 - time1 << '\n';
//		break;
//	case 3:
//		cout << "Time paral j " << threads_num << " : " << time2 - time1 << '\n';
//		break;
//	case 4:
//		cout << "Time paral ij " << threads_num << " : " << time2 - time1 << '\n';
//		break;
//	default:
//		break;
//	}
//	return matr;
//}
//
////int** matrix_multipl_par_two(int** matr1, int** matr2, int n, int threads_num){
////	int** matr = zero_matrix(n);
////	omp_set_num_threads(threads_num);
////	double time1, time2;
////	time1 = omp_get_wtime();
////	int i = 0;
////	int j = 0;
////#pragma omp parallel for  if(1 == 0) shared(matr, matr1, matr2) private(i)
////	for(i = 0; i < n; i++){
////		//#pragma omp parallel for shared(matr, matr1, matr2) private(j)
////		for(j = 0; j < n; j++){
////			for(int k = 0; k < n; k++){
////				matr[i][j] += matr1[i][k] * matr2[k][j];
////			}
////		}
////	}
////	time2 = omp_get_wtime();
////	res_time = ( time2 - time1 ) / (double) CLOCKS_PER_SEC;
////	cout << "Time paral(par_two) " << threads_num << ":" << time2 - time1 << '\n';
////	return matr;
////}
//
//int** matrix_multipl_block_par(int** matr1, int** matr2, int n, int block_size, int threads_num = 1, int mode = 1) {
//	int** matr = zero_matrix(n);
//	int block_num = n / block_size;
//	omp_set_num_threads(threads_num);
//	double time1, time2;
//	time1 = omp_get_wtime();
//	int i = 0, j = 0, k = 0, i1 = 0, j1 = 0, k1 = 0;
//#pragma omp parallel for if(mode == 2 || mode == 4)  private(j, k, i1, j1, k1)
//	for (i = 0; i < block_num; i++) {
//		for (j = 0; j < block_num; j++) {
//			for (k = 0; k < block_num; k++) {
//#pragma omp parallel  for if(mode == 3 || mode == 4) private(j1, k1)
//				for (i1 = i * block_size; i1 < (i + 1) * block_size; i1++) {
//					for (j1 = j * block_size; j1 < (j + 1) * block_size; j1++) {
//						for (k1 = k * block_size; k1 < (k + 1) * block_size; k1++) {
//							matr[i1][j1] += matr1[i1][k1] * matr2[k1][j1];
//						}
//					}
//				}
//			}
//		}
//	}
//	time2 = omp_get_wtime();
//	switch (mode)
//	{
//	case 1:
//		cout << "Time standart block " << threads_num << " block_num - " << block_num << " : " << time2 - time1 << '\n';
//		break;
//	case 2:
//		cout << "Time paral block i " << threads_num << " block_num - " << block_num << " : " << time2 - time1 << '\n';
//		break;
//	case 3:
//		cout << "Time paral block j " << threads_num << " block_num - " << block_num << " : " << time2 - time1 << '\n';
//		break;
//	case 4:
//		cout << "Time paral block ij " << threads_num << " block_num - " << block_num << " : " << time2 - time1 << '\n';
//		break;
//	default:
//		break;
//	}
//	return matr;
//}

void initMatr(int **matr1, int **matr2, int size)
{
	int left = -100;
	int right = 100;

	for (int i = 0; i < size; i++)
	{
		srand(time(0));

		for (int j = 0; j < size; j++)
		{
			matr1[i][j] = rand() % (right - left + 1) + left;
			matr2[i][j] = rand() % (right - left + 1) + left;
		}
	}
}

int  **multiply(int **matr1, int **matr2, int size, double &time)
{
	int **matrRes = new int *[size];
	for (int i = 0; i < size; i++)
	{
		matrRes[i] = new int[size];

		for (int j = 0; j < size; j++)
		{
			matrRes[i][j] = 0;
		}
	}

	double timeStart;
	double timeEnd;

	timeStart = clock();
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			for (int k = 0; k < size; k++)
			{
				matrRes[i][j] += matr1[i][k] * matr2[k][j];
			}
		}
	}
	timeEnd = clock();

	time = (timeEnd - timeStart) / 1000;

	return matrRes;
}

int **multiply(int **matr1, int **matr2, int size, int numberOfThreads, double &time)
{
	int **matrRes = new int *[size];
	for (int i = 0; i < size; i++)
	{
		matrRes[i] = new int[size];

		for (int j = 0; j < size; j++)
		{
			matrRes[i][j] = 0;
		}
	}

	numberOfThreads /= 2;

	omp_set_num_threads(numberOfThreads);
	omp_set_nested(true);

	double timeStart;
	double timeEnd;

	timeStart = omp_get_wtime();
	int i = 0;
	int j = 0;
	int k = 0;

#pragma omp parallel for private(j, k)
	for (i = 0; i < size; i++)
	{
#pragma omp parallel for private(k)
		for (j = 0; j < size; j++)
		{
			for (k = 0; k < size; k++)
			{
				matrRes[i][j] += matr1[i][k] * matr2[k][j];
			}
		}
	}
	timeEnd = omp_get_wtime();

	time = timeEnd - timeStart;

	return matrRes;
}

int **multiplyBlocks(int **matr1, int **matr2, int size, int blockSize, double &time)
{
	int **matrRes = new int *[size];
	for (int i = 0; i < size; i++)
	{
		matrRes[i] = new int[size];

		for (int j = 0; j < size; j++)
		{
			matrRes[i][j] = 0;
		}
	}

	int numberOfBlocks = size / blockSize;
	double timeStart;
	double timeEnd;

	timeStart = clock();
	for (int i = 0; i < numberOfBlocks; i++)
	{
		for (int j = 0; j < numberOfBlocks; j++)
		{
			for (int k = 0; k < numberOfBlocks; k++)
			{
				for (int i1 = i * blockSize; i1 < (i + 1) * blockSize; i1++)
				{
					for (int j1 = j * blockSize; j1 < (j + 1) * blockSize; j1++)
					{
						for (int k1 = k * blockSize; k1 < (k + 1) * blockSize; k1++)
						{
							matrRes[i1][j1] += matr1[i1][k1] * matr2[k1][j1];
						}
					}
				}
			}
		}
	}
	timeEnd = clock();

	time = (timeEnd - timeStart) / 1000;

	return matrRes;
}

int **multiplyBlocks(int **matr1, int **matr2, int size, int blockSize, int numberOfThreads, double &time)
{
	int **matrRes = new int *[size];
	for (int i = 0; i < size; i++)
	{
		matrRes[i] = new int[size];

		for (int j = 0; j < size; j++)
		{
			matrRes[i][j] = 0;
		}
	}

	int numberOfBlocks = size / blockSize;
	double timeStart;
	double timeEnd;

	omp_set_num_threads(numberOfThreads);
	
	timeStart = omp_get_wtime();
	int i = 0;
	int j = 0;
	int k = 0;
	int i1 = 0;
	int j1 = 0;
	int k1 = 0;
#pragma omp parallel for private(j, k, i1, j1, k1)
	for (i = 0; i < numberOfBlocks; i++) 
	{
		for (j = 0; j < numberOfBlocks; j++) 
		{
			for (k = 0; k < numberOfBlocks; k++)
			{
#pragma omp parallel  for private(j1, k1)
				for (i1 = i * blockSize; i1 < (i + 1) * blockSize; i1++) 
				{
					for (j1 = j * blockSize; j1 < (j + 1) * blockSize; j1++) 
					{
						for (k1 = k * blockSize; k1 < (k + 1) * blockSize; k1++)
						{
							matrRes[i1][j1] += matr1[i1][k1] * matr2[k1][j1];
						}
					}
				}
			}
		}
	}
	timeEnd = omp_get_wtime();

	time = timeEnd - timeStart;
	
	return matrRes;
}

void run(int **matr1, int **matr2, int size, string fileName)
{
	double time;
	int **matrRes;

	ofstream fout;
	fout.open(fileName);

	fout << "Size: " << size << endl;

	initMatr(matr1, matr2, size);
	matrRes = multiply(matr1, matr2, size, time);
	fout << "Point multiplication time: " << endl << time << endl;

	fout << "Point multiplication time parallel: " << endl;
	for (int numberOfThreads = 2; numberOfThreads <= 10; numberOfThreads += 2)
	{
		matrRes = multiply(matr1, matr2, size, numberOfThreads, time);
		fout << time << endl;
	}

	fout << "Block multiplication time: " << endl;
	for (int blockSize = 10; blockSize < size / 2; blockSize += 10)
	{
		matrRes = multiplyBlocks(matr1, matr2, size, blockSize, time);
		fout << time << endl;
	}

	fout << "Block multiplication time parallel: " << endl;
	for (int numberOfThreads = 2; numberOfThreads <= 10; numberOfThreads += 2)
	{
		fout << "number of threads: " << numberOfThreads << endl;
		for (int blockSize = 10; blockSize < size / 2; blockSize += 10)
		{
			matrRes = multiplyBlocks(matr1, matr2, size, blockSize, time);
			fout << time << endl;
		}
	}

	fout.close();

	for (int i = 0; i < size; i++)
	{
		delete matr1[i];
		delete matr2[i];
	}

	delete matr1;
	delete matr2;
}

int main()
{
	int size1 = 500;
	int size2 = 1500;

	int **matr1;
	int **matr2;

	string fileName1 = "Test1.csv";
	string fileName2 = "Test2.csv";

	matr1 = new int *[size1];
	matr2 = new int *[size1];
	for (int i = 0; i < size1; i++)
	{
		matr1[i] = new int[size1];
		matr2[i] = new int[size1];
	}
	run(matr1, matr2, size1, fileName1);

	matr1 = new int *[size2];
	matr2 = new int *[size2];
	for (int i = 0; i < size2; i++)
	{
		matr1[i] = new int[size2];
		matr2[i] = new int[size2];
	}
	run(matr1, matr2, size2, fileName2);

	cout << endl;
	system("pause");
	return 0;
}