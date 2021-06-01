
#include <iostream>
#include <cmath>
#include <vector>
#include <functional>
using namespace std;
#include "Gaus.h"
typedef function<double(double x0)> func;
typedef function<double(double)> Res;
double f(double x)
{
    return exp(x);
}
double fcalc(vector<double> koef, double x)
{
    double res = 0;
    for (int i = 0; i < koef.size(); i++) res += koef[i] * pow(x, i);
    return res;
}
double SLAR(vector<double> x, vector<double> y, double x0)
{
    const int size = x.size();
    vector < vector < double>> matrix(size);
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            matrix[j].resize(size);
            matrix[j][i] = pow(x[i], j);
        }
    }
    vector<double> ans = Gaus(matrix, y);
    return fcalc(ans, x0);
}
Res GenereteLagr(vector<double> x, vector<double> y)
{
    const int size = x.size();
    vector<func> f(size);
    for (int i = 0; i < size; i++)
    {
        f[i] = [i, x](double x0)
        {
            double sum = 1;
            for (int j = 0; j < x.size(); j++)
                if (j != i) sum *= (x0 - x[j]) / (x[i] - x[j]);
            return sum;
        };
    }
    Res L = [f, y](double x)
    {
        int size = f.size();
        double asw = 0;
        for (int i = 0; i < size; i++)
        {
            asw += f[i](x) * y[i];
        }
        return asw;
    };
    return L;
}
Res GenereteNtom(Res nton, vector<double> x, vector<double> y, int size)
{
    double A = 0;
    for (int i = 0; i < size; i++)
    {
        double  w = 1;
        for (int j = 0; j < size; j++)
            if (i != j) w *= (x[i] - x[j]);
        A += y[i] / w;
    }
    Res w = [x, size](double x0)
    {
        double res = 1;
        for (int i = 0; i < size - 1; i++) res *= (x0 - x[i]);
        return res;
    };
    Res res = [nton, A, w](double x)
    {
        double old = nton(x);
        double ww = w(x);
        double res = old + A * ww;
        return res;
    };
    return res;
}
int main()
{
    int a = 1, b = 10, n = 20;
    double x0 = 5;
    //cin >> a >> b >> n;
    vector<double> x(n);
    vector<double> y(n);
    for (int i = 0; i < n; i++)
    {
        x[i] = a + (double)(b - a) * i / (n - 1);
        y[i] = f(x[i]);
    }
    auto L = GenereteLagr(x, y);
    Res prev = [](double x) {return 0; };
    for (int i = 1; i < x.size(); i++)
    {
        prev = GenereteNtom(prev, x, y, i);
    }
    cout << SLAR(x, y, x0) << endl << L(x0) << endl << prev(x0) << endl;
}

