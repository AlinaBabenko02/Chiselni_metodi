#include <iostream>
#include <iomanip>
#include <vector>
#include <ctime>
#include <cmath>
#include <thread>
using namespace std;

int SIZE = 100;

#define eps 1.0/10000

ostream& operator<<(ostream& stream, vector<double> vec)
{
    stream << "(";
    for (int i = 0; i < SIZE - 1; i++)
    {
        stream << round(vec[i] / (eps)) * eps << ", ";
    }
    stream << round(vec[SIZE - 1] / (eps)) * eps << ")\n";
    return stream;
}


void PrintLine(ostream& stream)
{
    for (int i = 0; i < 120; i++)
    {
        stream << "-";
    }
    stream << endl;
}


vector<vector<double>> operator*(vector<vector<double>> vec1, vector<vector<double>> vec2)
{
    vector<vector<double>> res(SIZE);
    for (auto& c : res) c.resize(SIZE);
    for (int i = 0; i < SIZE; i++)
    {
        for (int j = 0; j < SIZE; j++)
        {
            res[j][i] = 0;
            for (int t = 0; t < SIZE; t++)
            {
                res[j][i] += vec1[t][i] * vec2[j][t];
            }
        }
    }
    return res;
}

ostream& operator<<(ostream& stream, const vector<vector<double>>& matrix)
{
    for (int i = 0; i < SIZE; i++)
    {
        for (int j = 0; j < SIZE; j++)
        {
            stream << setw(20) << setfill(' ') << matrix[j][i];
        }
        stream << endl;
    }
    PrintLine(stream);
    return stream;
}

double determ(vector<vector<double>> matrix, int size)
{
    int i, j;
    double det = 0;
    vector<vector<double>> matr(size - 1);
    if (size == 1)
    {
        det = matrix[0][0];
    }
    else if (size == 2)
    {
        det = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    }
    else
    {
        for (i = 0; i < size; ++i)
        {
            for (j = 0; j < size - 1; ++j)
            {
                if (j < i)
                    matr[j] = matrix[j];
                else
                    matr[j] = matrix[j + 1];
            }
            det += pow((double)-1, (i + j)) * determ(matr, size - 1) * matrix[i][size - 1];
        }
    }
    return det;
}

//метод Гауса
void Gaus(vector<vector<double>> matrix, vector<double> b)
{
    int SIZE = matrix.size();
    vector<double> answ(SIZE);
    vector<vector<double>> U(SIZE);
    vector<double> m(SIZE);
    unsigned int start = clock();
    for (int i = 0; i < SIZE; i++)
    {
        U[i].resize(SIZE);
        U[i][i] = 1.0;
    }
    //matrix = { {1,2,7},{1,6,4},{1,8,3} };
    for (int t = 0; t < SIZE; t++)
    {
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
    unsigned int end = clock();
    cout << "Gaus\n";
    cout << "Answer = " << answ
        << "Time = " << end - start << endl;
    PrintLine(cout);
}
//метод Якоби
void Jakobi(vector<vector<double>> matrix, vector<double> b)
{
    vector<double> answ(SIZE);
    unsigned int start = clock();
    for (int i = 0; i < SIZE; i++) answ[i] = matrix[i][i] / b[i];
    int counter = 0;
    double delta = 0;
    vector<double> new_answ(SIZE);
    do
    {
        counter++;
        for (int i = 0; i < SIZE; i++)
        {
            new_answ[i] = b[i];
            for (int j = 0; j < SIZE; j++)
            {
                if (i != j)
                    new_answ[i] -= matrix[j][i] * answ[j];
            }
            new_answ[i] /= matrix[i][i];
            delta = fabs(answ[0] - new_answ[0]);
            for (int i = 0; i < SIZE; i++)
            {
                if (fabs(answ[i] - new_answ[i]) > delta)
                    delta = fabs(answ[i] - new_answ[i]);
                answ[i] = new_answ[i];
            }
        }
    } while (delta >= eps);
    unsigned int end = clock();
    cout << "Jakobi" << endl
        << "Answer = " << answ
        << "Time = " << end - start << endl
        << "Counter = " << counter << endl;
    PrintLine(cout);
}
//метод Зейделя
void Zeidel(vector<vector<double>> matrix, vector<double> b)
{
    vector<double> answ(SIZE);
    unsigned int start = clock();
    for (int i = 0; i < SIZE; i++) answ[i] = matrix[i][i] / b[i];
    int counter = 0;
    double delta = 0;
    vector<double> new_answ = answ;
    do
    {
        counter++;
        for (int i = 0; i < SIZE; i++)
        {
            new_answ[i] = 0;
            for (int j = 0; j < SIZE; j++)
            {
                if (i != j)
                    new_answ[i] -= matrix[j][i] * new_answ[j];
            }
            new_answ[i] += b[i];
            new_answ[i] /= matrix[i][i];
            delta = fabs(answ[0] - new_answ[0]);
            for (int i = 0; i < SIZE; i++)
            {
                if (fabs(answ[i] - new_answ[i]) > delta)
                    delta = fabs(answ[i] - new_answ[i]);
                answ[i] = new_answ[i];
            }
        }
    } while (delta >= eps);
    unsigned int end = clock();
    cout << "Zeidel" << endl
        << "Answer = " << answ
        << "Time = " << end - start << endl
        << "Counter = " << counter << endl;
    PrintLine(cout);
}

pair<vector<vector<double>>, vector<double>> GenereteRandomSys(int n)
{
    srand(time(0));
    vector<vector<double>> res(n);
    vector<double> b(n);
    for (int i = 0; i < n; i++)
    {
        res[i].resize(n);
        for (int j = 0; j < n; j++)
            res[i][j] = rand() % 2001 - 1000;
        b[i] = rand() % 2001 - 1000;
    }
    return make_pair(res, b);
}

pair<vector<vector<double>>, vector<double>> GenereteRandomOk(int n)
{
    srand(time(0));
    vector<vector<double>> res(n);
    vector<double> b(n);
    for (int i = 0; i < n; i++) res[i].resize(n);
    for (int i = 0; i < n; i++)
    {
        double sum = 0;
        for (int j = 0; j < n; j++)
        {
            if (i != j)
            {
                res[j][i] = rand() % 2001 - 1000;
                sum += fabs(res[j][i]);
            }

        }
        res[i][i] = sum * pow(-1, rand() % 2) + rand() % 201 - 100;
        b[i] = rand() % 2001 - 1000;
    }
    return make_pair(res, b);
}

pair<vector<vector<double>>, vector<double>> GenereteGilbert(int n)
{
    vector<vector<double>> matr(n);
    vector<double> b(n);
    srand(time(0));
    for (int i = 0; i < n; i++)
    {
        b[i] = 0;
        matr[i].resize(n);
        for (int j = 0; j < n; j++)
        {
            matr[i][j] = 1.0 / (i + j + 1);
        }
        b[i] = (rand() % 2001 - 1000) / 10000.0;
    }

    return make_pair(matr, b);
}
int main()
{
    SIZE = 5;
    auto sys = GenereteRandomSys(SIZE);
    cout << sys.first << sys.second;
    Gaus(sys.first, sys.second);
    Jakobi(sys.first, sys.second);
    Zeidel(sys.first, sys.second);
    return 0;
}