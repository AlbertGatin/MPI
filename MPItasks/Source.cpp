#include<iostream>
#include<mpi.h>
#include "windows.h"

using namespace std;
//Hello world из всех процессов (2 балла)
//void task1(int argc, char *argv[])
//{
//	int size = 0;
//	int rank = 0;
//	MPI_Init(&argc, &argv);
//	MPI_Comm_size(MPI_COMM_WORLD, &size);
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//	printf("Hello world! size = %d, rank = %d \n", size, rank);
//	MPI_Finalize();
//}
////Max вектора  (3 балла)
//void task2(int argc, char *argv[], int vectSize)
//{
//	int *x = new int[vectSize];
//	int rank, size;
//	int result;
//	MPI_Init(&argc, &argv);
//	MPI_Comm_size(MPI_COMM_WORLD, &size);
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//	if (rank == 0)
//	{
//		for (int i = 0; i < vectSize; i++)
//		{
//			x[i] = rand() % 100;
//			printf("%d ", x[i]);
//		}
//	}
//	int partSize = vectSize / size;
//	int *partialVect = new int[partSize];
//	int privateMax;
//
//	MPI_Scatter(x, partSize, MPI_INT,
//		partialVect, partSize, MPI_INT, 0, MPI_COMM_WORLD);
//	privateMax = 0;
//	for (int i = 0; i < partSize; i++)
//		if (partialVect[i] > privateMax) privateMax = partialVect[i];
//	MPI_Reduce(&privateMax, &result, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
//	if (rank == 0)
//	{
//		printf("\nmax = %d\n", result);
//		delete x;
//		delete partialVect;
//	}
//	MPI_Finalize();
//	
//}
//
////Вычисление числа Пи методом Монте - Карло(3 баллов)
//void task3(int argc, char *argv[])
//{
//	const int pointCount = 50000;
//	int rank, size;
//	double point[2][pointCount];
//	double x, y;
//	int result = 0;
//
//	MPI_Init(&argc, &argv);
//	MPI_Comm_size(MPI_COMM_WORLD, &size);
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//	if (rank == 0)
//	{
//		for (int i = 0; i < pointCount; i++)
//		{
//			x = (double)(rand() % 1000) / 1000;
//			y = (double)(rand() % 1000) / 1000;
//			point[0][i] = x;
//			point[1][i] = y;
//		}
//	}
//	MPI_Bcast(point, 2*pointCount, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//
//	int cuttedELems = pointCount / size;
//	int i1 = rank * cuttedELems;
//	int i2 = (rank + 1) * cuttedELems;
//	if (rank == size - 1)
//		i2 = pointCount;
//
//	int count = 0;
//
//	for (int i = i1; i < i2; i++)
//	{
//		if ((point[0][i] * point[0][i]) + (point[1][i] * point[1][i]) < 1)
//			count++;
//	}
//
//	if (rank != 0)
//		MPI_Send(&count, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
//	else
//	{
//		int temp;
//		result = count;
//		for (int i = 1; i < size; i++)
//
//		{
//			MPI_Recv(&temp, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//			result += temp;
//		}
//		printf("\nPoint count = %d \n", pointCount);
//		printf("Point count in cirvle = %d \n", result);
//		printf("Pi = %f", 4.000 / pointCount * result);
//	}
//	MPI_Finalize();
//}
//Скалярное произведение (3 балла)
//void task5(int argc, char *argv[], int N)
//{
//	int *x = new int[N];
//	int *y = new int[N];
//	int rank, size;
//
//	MPI_Init(&argc, &argv);
//	MPI_Comm_size(MPI_COMM_WORLD, &size);
//	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
//	int partCount = N / size;
//	int *partialArr1 = new int[partCount];
//	int *partialArr2 = new int[partCount];
//
//	if (rank == 0) {
//		for (int i = 0; i < N; i++)
//		{
//			x[i] = rand() % 100;
//			printf("%d ", x[i]);
//		}
//		printf("\n");
//		for (int i = 0; i < N; i++)
//		{
//			y[i] = rand() % 100;
//			printf("%d ", y[i]);
//		}
//		printf("\n");
//	}
//	MPI_Scatter(x, partCount, MPI_INT,
//		partialArr1, partCount, MPI_INT, 0, MPI_COMM_WORLD);
//	MPI_Scatter(y, partCount, MPI_INT,
//		partialArr2, partCount, MPI_INT, 0, MPI_COMM_WORLD);
//
//	int sum = 0;
//	for (int i = 0; i < partCount; i++)
//		sum += partialArr1[i] * partialArr2[i];
//
//	if (rank != 0)
//	{
//		MPI_Send(&sum, 1, MPI_INT, 0, 0, MPI_COMM_WORLD);
//	}
//	else
//	{
//		int temp = 0;
//		for (int i = 1; i < size; i++)
//		{
//			MPI_Recv(&temp, 1, MPI_INT, MPI_ANY_SOURCE, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
//			sum += temp;
//		}
//		printf("result = %d\n", sum);
//	}
//	MPI_Finalize();
//}

//Среднее арифметическое среди положительных чисел массива(3 баллов)
void task4(int argc, char *argv[], int N)
{
	int *x = new int[N];
	int rank, size;
	int resultSum;
	int resultCount;


	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	int partCount = N / size;
	int *partialArr = new int[partCount];

	int sum;
	int count;
	
	if (rank == 0) {
		for (int i = 0; i < N; i++)
		{
			x[i] = rand() % 100;
			printf("%d ", x[i]);
		}
	}
	MPI_Scatter(x, partCount, MPI_INT,
		partialArr, partCount, MPI_INT, 0, MPI_COMM_WORLD);

	sum = 0;
	count = 0;
	for (int i = 0; i < partCount; i++)
		if (partialArr[i] > 0)
		{
			sum += partialArr[i];
			count++;
		}

	MPI_Reduce(&sum, &resultSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Reduce(&count, &resultCount, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	if (rank == 0)
	{
		printf("\navg = %f\n", (double)resultSum / resultCount);
		delete x;
	}
	MPI_Finalize();
}
//Maxmin матрицы (5 балла)  
void task6(int argc, char *argv[])
{

	int x[10][10];  //our matrix
	int ProcRank, ProcNum, N = 10;
	int resultMin;
	int resultMax;

	// trivial stuff
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	int elements_per_proc = N * N / ProcNum;
	int *subarr = new int[elements_per_proc];

	int localmax;
	int localmin;
	//	srand(time(NULL));
	if (ProcRank == 0) {
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++) {
				x[i][j] = rand() % 100;
				x[0][0] = 99;
				x[N - 1][N - 1] = 0;

				printf(" %d ", x[i][j]);
			}
			printf("\n");
		}
	}
	// Distribute the array
	MPI_Scatter(x, elements_per_proc, MPI_INT,
		subarr, elements_per_proc, MPI_INT, 0, MPI_COMM_WORLD);

	// Find the maximum element of the local subarray
	localmax = 0;
	for (int i = 0; i < elements_per_proc; i++)
		if (subarr[i] > localmax) localmax = subarr[i];
	localmin = 100;
	for (int i = 0; i < elements_per_proc; i++)
		if (subarr[i] < localmin) localmin = subarr[i];

	// Perform global max reduction
	MPI_Reduce(&localmin, &resultMin, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
	MPI_Reduce(&localmax, &resultMax, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
	if (ProcRank == 0)
		printf("\nmin = %d\nmax = %d", resultMin, resultMax);
	MPI_Finalize();

}

//через скаттерв и 1 редьюс
void task4alt(int argc, char *argv[], int N)
{
	int *x = new int[N];
	int rank, size;
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);
	int recv_buf[100];
	int *send_counts = (int*)malloc(size * sizeof(int));
	int *interval_begins = (int*)malloc(size * sizeof(int));
	int rem = N % size;
	int part_size = 0;
	int avg[2] = { 0 , 0 };
	int resAvg[2];
	if (rank == 0) {
		for (int i = 0; i < size; i++)
		{
			send_counts[i] = N / size;
			if (rem > 0)
			{
				send_counts[i]++;
				rem--;
			}

			interval_begins[i] = part_size;
			part_size += send_counts[i];
		}
		
		for (int i = 0; i < N; i++)
		{
			x[i] = rand() % 100;
			printf("%d ", x[i]);
		}
		printf("\n");
		
		for (int i = 0; i < size; i++)
		{
			printf("send_counts[%d] = %d\tinterval_begins[%d] = %d\n ", i, send_counts[i], i, interval_begins[i]);
		}
	}
	//MPI_Barrier(MPI_COMM_WORLD);
	MPI_Scatterv(x, send_counts, interval_begins, MPI_INT, recv_buf, 100, MPI_INT, 0, MPI_COMM_WORLD);

	avg[0] = 0;
	avg[1] = 0;
	for (int i = 0; i < 100; i++)
	{
			if (recv_buf[i] > 0)
			{
				avg[0] += recv_buf[i];
				avg[1] ++;
			}
	}
	MPI_Reduce(&avg, &resAvg, 2, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	if (rank == 0)
	{
		printf("\navg = %f\n", (double)resAvg[0]/resAvg[1]);
		delete x, recv_buf, send_counts, interval_begins;
	}
	MPI_Finalize();

	
}

//через скаттер в передача матрицы по строкам
void task6alt(int argc, char *argv[])
{
	int const heigth = 5;
	int const width = 5;
	int x[heigth][width];  
	int ProcRank, ProcNum;
	int maxmin;
	int minmax;

	int recv_buf[100];	

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	int *send_counts = (int*)malloc(ProcNum * sizeof(int));
	int *interval_begins = (int*)malloc(ProcNum * sizeof(int));
	int rem = heigth % ProcNum;
	int part_size = 0;
	int proc_maxmin = 0;
	int proc_minmax = 100;
	int min = 100;
	int max = 0;
	if (ProcRank == 0) {
		for (int i = 0; i < ProcNum; i++)
		{
			send_counts[i] = (heigth / ProcNum) * width;
			if (rem > 0)
			{
				send_counts[i]+=width;
				rem--;
			}

			interval_begins[i] = part_size;
			part_size += send_counts[i];
		}

		for (int i = 0; i < heigth; i++)
		{
			for (int j = 0; j < width; j++) {
				x[i][j] = rand() % 100;
				printf(" %d ", x[i][j]);
			}
			printf("\n");
		}

		printf("\n");

		for (int i = 0; i < ProcNum; i++)
		{
			printf("send_counts[%d] = %d\tinterval_begins[%d] = %d\n ", i, send_counts[i], i, interval_begins[i]);
		}
	}
	MPI_Scatterv(x, send_counts, interval_begins, MPI_INT, recv_buf, 100, MPI_INT, 0, MPI_COMM_WORLD);
    proc_maxmin = 0;
	proc_minmax = 100;
	min = 100;
	max = 0;
	for (int i = 0; i < 100; i+=width)
	{
		for (int j = i; j < width + i; j++)
		{
			if (recv_buf[j] < min)
				min = recv_buf[j];
			if (recv_buf[j] > max)
				max = recv_buf[j];
		}
		if (min > proc_maxmin)
			proc_maxmin = min;
		if (max < proc_minmax)
			proc_minmax = max;
		if (recv_buf[i + 2 * width] == NULL)
			break;
	}

	// Perform global max reduction
	MPI_Reduce(&proc_maxmin, &maxmin, 1, MPI_INT, MPI_MAX, 0, MPI_COMM_WORLD);
	MPI_Reduce(&proc_minmax, &minmax, 1, MPI_INT, MPI_MIN, 0, MPI_COMM_WORLD);
	if (ProcRank == 0)
		printf("\nmaxmin = %d\nminmax = %d", maxmin, minmax);
	MPI_Finalize();

}

//умножение матрицы на вектор
void task7(int argc, char *argv[])
{
	int x[10][10];
	int y[10];
	int ProcRank, ProcNum, N = 10;
	int resultMin;
	int resultMax;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	int elements_per_proc = N * N / ProcNum;
	int *subarr1 = new int[elements_per_proc];

	if (ProcRank == 0) {
		for (int i = 0; i < N; i++)
		{
			for (int j = 0; j < N; j++) {
				x[i][j] = rand() % 100;
				printf(" %d ", x[j][i]);
			}
			printf("\n");
		}
		
		for (int i = 0; i < N; i++)
		{
			y[i] = rand() % 100;
			printf(" %d ", y[i]);
		}
		printf("\n\n");
	}
	
	MPI_Scatter(x, elements_per_proc, MPI_INT,
		subarr1, elements_per_proc, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(y, N, MPI_INT, 0, MPI_COMM_WORLD);

	for (int i = 0; i < elements_per_proc; i++)
		subarr1[i] *= y[i % N];

	// Perform global max reduction
	MPI_Gather(subarr1, elements_per_proc, MPI_INT, x, elements_per_proc, MPI_INT, 0, MPI_COMM_WORLD);
	if (ProcRank == 0)
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < N; j++) {
				printf(" %d ", x[j][i]);
			}
			printf("\n");
		}
	MPI_Finalize();

}

//Scatter и Gather через Send и Recv
void task8(int argc, char *argv[])
{
	int x[10];
	int ProcRank, ProcNum, N = 10;

	// init and stuff
	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &ProcNum);
	MPI_Comm_rank(MPI_COMM_WORLD, &ProcRank);
	int elements_per_proc = N / ProcNum;
	int *subarr1 = new int[elements_per_proc];

	//srand(time(NULL));
	if (ProcRank == 0) {
		for (int i = 0; i < N; i++)
		{
			x[i] = rand() % 100;
			printf(" %d ", x[i]);
		}
		printf("\n");
	}

	if (ProcRank == 0) {
		MPI_Sendrecv(x, elements_per_proc, MPI_INT, 0, 0, subarr1, elements_per_proc, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		for (int i = elements_per_proc; i < N; i += elements_per_proc) {
			MPI_Send(x + i, elements_per_proc, MPI_INT, i / elements_per_proc, 0, MPI_COMM_WORLD);
		}
	}
	else {
		MPI_Recv(subarr1, elements_per_proc, MPI_INT, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}

	for (int i = 0; i < elements_per_proc; i++) {
		printf(" %d ", subarr1[i]);
	}
	printf(" : from process %d\n", ProcRank);




	MPI_Finalize();
}

int main(int argc, char *argv[])
{
	//Sleep(15000);
	//task1(argc, argv);
	//task2(argc, argv, 20);
	//task3(argc, argv);
	//task5(argc, argv, 10);
	//task4(argc, argv, 20);
	//task6(argc, argv);


	//task4alt(argc, argv, 10);
	task6alt(argc, argv);
	//task7(argc, argv);
	//task8(argc, argv);

	//system("pause");
	return 0;
}