#include "stdlib.h"
#include <random>
#include <iostream>
#include <omp.h>
#include <ctime>
#include <algorithm>
#include <fstream>
#include <iomanip>

using namespace std;

double startTime = -1;

int **randMatr(int size)
{
	int left = -100;
	int right = 100;
	int **matr = new int*[size];

	for (int i = 0; i < size; i++)
	{
		matr[i] = new int[size];

		srand(time(0));

		for (int j = 0; j < size; j++)
		{
			matr[i][j] = rand() % (right - left + 1) + left;
		}
	}

	return matr;
}

int **zeroMatr(int size)
{
	int **matr = new int*[size];

	for (int i = 0; i < size; i++)
	{
		matr[i] = new int[size];

		for (int j = 0; j < size; j++)
		{
			matr[i][j] = 0;
		}
	}

	return matr;
}

void removeMatr(int **matr, int size)
{
	for (int i = 0; i < size; i++)
	{
		delete matr[i];
	}

	delete matr;
}

void setStartTime()
{
	startTime = clock();
}

double getTime()
{
	return (clock() - startTime) / CLOCKS_PER_SEC;
}

double multiply(int size, int threadNum)
{
	double time = -1;
	int **matr1 = randMatr(size);
	int	**matr2 = randMatr(size);
	int	**resMatr = zeroMatr(size);

	omp_set_num_threads(threadNum);
	setStartTime();
#pragma omp parallel for shared(matr1, matr2, resMatr)
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			for (int k = 0; k < size; k++)
			{
				resMatr[i][j] += matr1[i][k] * matr2[j][k];
			}
		}
	}
	time = getTime();

	removeMatr(matr1, size);
	removeMatr(matr2, size);
	removeMatr(resMatr, size);

	return time;
}

double multiply(int size, int r, int threadNum)
{
	double time = -1;
	int **matr1 = randMatr(size);
	int **matr2 = randMatr(size);
	int **resMatr = zeroMatr(size);
	int q = ceil(1.0 * size / r);
	int ind1;
	int ind2;
	int ind3;
	int ind4;
	int ind5;

	omp_set_num_threads(threadNum);
	setStartTime();
#pragma omp parallel for shared(matr1, matr2, resMatr)
	for (int ij = 0; ij < q * q; ij++)
	{
		ind1 = ij / q;
		ind2 = ij % q;

		for (int k1 = 0; k1 < q; k1++)
		{
			ind3 = min((ind1 + 1)*r, size);
			ind4 = min((ind2 + 1)*r, size);
			ind5 = min((k1 + 1)*r, size);

			for (int i = ind1 * r; i < ind3; i++)
			{
				for (int j = ind2 * r; j < ind4; j++)
				{
					for (int k = k1 * r; k < ind5; k++)
					{
						resMatr[i][j] += matr1[i][k] * matr2[j][k];
					}
				}
			}
		}
	}
	time = getTime();

	removeMatr(matr1, size);
	removeMatr(matr2, size);
	removeMatr(resMatr, size);

	return time;
}

void run(vector<int> Ns, vector<int> Rs, vector<int> Threads, string fileName)
{
	double time;
	int **matrRes;

	ofstream fout;
	fout.open(fileName);

	for (int n : Ns) {
		fout << "Matrix size: " << n << endl;

		fout << "Threads";
		for (int threadNum : Threads)
		{
			fout << '\t' << threadNum;
		}
		fout << endl;

		for (int r : Rs)
		{
			if (r > n / 2 && r != n)
			{
				continue;
			}

			fout << r;
			for (int threadNum : Threads)
			{

				if (r == 0)
				{
					time = multiply(n, threadNum);
				}
				else
				{
					time = multiply(n, r, threadNum);
				}

				fout << '\t' << setprecision(3) << time;
			}
			fout << endl;
		}
		fout << endl;
	}

	fout.close();
}

int main()
{
	string fileName = "Test.txt";

	vector<int> Ns({ 500,1000,1500 });
	vector<int> Rs({ 0,1,2,5,10,50,100,200,300,400,500,600,700,800,900,1000 });
	vector<int> Threads({ 1,2,4,8,16 });

	run(Ns, Rs, Threads, fileName);

	cout << endl;
	system("pause");
	return 0;
}