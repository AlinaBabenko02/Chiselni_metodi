#include <iostream>
#include <vector>
#include <cstring>
using namespace std;
double d = 0.85;
struct edge
{
    int FirstDot;
    int SecondDot;
};
struct dot
{
    int number;
    double rank;
};
ostream& operator<<(ostream& stream, vector<dot> dots)
{
    for (int i = 0; i < dots.size(); i++)
        stream << dots[i].number << " " << dots[i].rank << endl;
    return stream;
}
istream& operator>>(istream& stream, edge& e)
{
    stream >> e.FirstDot;
    stream >> e.SecondDot;
    return stream;
}

void Iteration(vector<dot>& dots, const vector<edge>& edges)
{
    vector<dot> res = dots;
    for (int i = 0; i < dots.size(); i++)
    {
        res[i].rank = 1 - d;
        for (int j = 0; j < edges.size(); j++)
        {
            int c = 0;
            if (edges[j].SecondDot == dots[i].number)
            {
                for (int g = 0; g < edges.size(); g++)
                {
                    if (edges[g].FirstDot == edges[j].FirstDot)
                        c++;
                }
                res[i].rank += d * dots[edges[j].FirstDot - 1].rank / c;
            }
        }
    }
    dots = res;
}


int main()
{
    setlocale(LC_ALL, "rus");
    int n;
    cout << "Введите количество вершин\n";
    cin >> n;
    vector<dot> dots(n);
    vector<edge> edges;
    for (int i = 0; i < n; i++)
    {
        dots[i].number = i + 1;
        dots[i].rank = 1;
    }
    cout << "Введите ребра в виде пар вершин(первая вершина - начало, вторая - конец)\n";
    int first;
    edge temp;
    while (cin >> temp)
    {
        edges.push_back(temp);
    }
    for (int i = 0; i < 50; i++)
        Iteration(dots, edges);
    cout << dots;
    return 0;
}