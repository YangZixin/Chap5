#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <fstream>
using namespace std;

const int LengthOfBox = 1;
const int GridInterval = 40;
const int GridCheck = LengthOfBox * GridInterval +1;
const double PI = 3.14159;
const double alpha = 2 / (1 + PI / LengthOfBox);
const int SetPlate = 0.3 * GridInterval;
double potential[GridCheck][GridCheck];
double temp_potential[GridCheck][GridCheck];
double delta_V;
double PermitErrorV = GridCheck * GridCheck * pow(10, -5);
double rhodensity[GridCheck][GridCheck];



void initializa_rho()
{
    rhodensity[GridCheck/2][GridCheck/2] = 1;
}

void update_V()
{
    delta_V = 0;
    for (int i = 1; i < GridCheck-1; i++)
    {
        for (int j = 1; j < GridCheck-1; j++)
        {
            temp_potential[i][j] = (potential[i-1][j]+potential[i+1][j]+potential[i][j-1]+potential[i][j+1]) / 4 + rhodensity[i][j] * GridInterval * GridInterval / 4;
            delta_V += sqrt((potential[i][j]-temp_potential[i][j]) * (potential[i][j]-temp_potential[i][j]));
        }

    }
    for (int i = 1; i < GridCheck-1; i++)
    {
        for (int j = 1; j < GridCheck-1; j++)
        {
            
                potential[i][j]=temp_potential[i][j];
        }
    }
}

void update_V_SOR()
{
    delta_V = 0;
    for (int i = 1; i < GridCheck-1; i++)
    {
        for (int j = 1; j < GridCheck-1; j++)
        {
            temp_potential[i][j] = (potential[i-1][j]+potential[i+1][j]+potential[i][j-1]+potential[i][j+1]) / 4 + rhodensity[i][j] * GridInterval * GridInterval / 4;
            potential[i][j] = alpha * (temp_potential[i][j] - potential[i][j]) + potential[i][j];
            delta_V += alpha * sqrt((potential[i][j]-temp_potential[i][j]) * (potential[i][j]-temp_potential[i][j]));
        }
    }
    for (int i = 1; i < GridCheck-1; i++)
    {
        for (int j = 1; j < GridCheck-1; j++)
        {
            
                potential[i][j]=temp_potential[i][j];
        }
    }
}

int laplace_calculate()
{
    int N_iter;
    do
    {
        update_V();
        N_iter++;
    }
    while(delta_V > PermitErrorV);
    return N_iter;
}

int laplace_calculate_SOR()
{
    int N_iter_SOR;
    do
    {
        update_V_SOR();
        N_iter_SOR++;
    }
    while(delta_V > PermitErrorV);
    return N_iter_SOR;
}
int main()
{
    initializa_rho();
    int N_iter = laplace_calculate_SOR();
    cout<<N_iter<<endl;
    for (int i = 0; i < GridCheck; ++i)
    {
        for (int j = 0; j < GridCheck; ++j)
        {
            cout<<potential[i][j]<<"    ";
        }
        cout<<endl;
    }

    return 0;
}