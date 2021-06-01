#include "Gaus.h"
vector<double> Gaus(vector<vector<double>> matrix, vector<double> b)
{
	int SIZE = b.size();
	vector<double> answ(SIZE);
	vector<vector<double>> U(SIZE);
	vector<double> m(SIZE);
	for (int i = 0; i < SIZE; i++)
	{
		U[i].resize(SIZE);
		U[i][i] = 1.0;
	}
	//matrix = { {1,2,7},{1,6,4},{1,8,3} };
	for (int t = 0; t < SIZE; t++)
	{
		for (int i = t; i < SIZE; i++)
		{
			//if (fabs(matrix[t][i]) > fabs(matrix[t][t]))
			//{
			//	for(int j=0; j<SIZE;j++)
			//	swap(matrix[j][t], matrix[i][t]);
			//	swap(b[t], b[i]);
			//}
		}
		for (int i = t + 1; i < SIZE; i++)
		{
			m[i] = matrix[i][t] / matrix[t][t];
			U[i][t] = m[i];
			for (int j = t; j < SIZE; j++)
			{
				matrix[i][j] = matrix[i][j] - matrix[t][j] * m[i];
			}
		}
	}
	vector<double> y(SIZE);
	for (int i = 0; i < SIZE; i++)
	{
		double sum = 0;
		for (int j = 0; j < i; j++)
			sum += y[j] * matrix[j][i];
		y[i] = (b[i] - sum) / matrix[i][i];
	}
	for (int i = SIZE - 1; i >= 0; i--)
	{
		double sum = 0;
		for (int j = i + 1; j < SIZE; j++)
			sum += answ[j] * U[j][i];
		answ[i] = (y[i] - sum) / U[i][i];
	}
	return answ;
}