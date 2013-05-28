#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <fstream>
using namespace std;

const int LengthOfCapacitor = 2;
const int GridInterval = 10;
const int GridCheck = LengthOfCapacitor * GridInterval +1;
const double PI = 3.14159;
const double alpha = 2 / (1 + PI / LengthOfCapacitor);
const int SetPlate = 0.3 * GridInterval;
double potential[GridCheck][GridCheck];
double temp_potential[GridCheck][GridCheck];
double delta_V;
double PermitErrorV = GridCheck * GridCheck * pow(10, -5);

void initialize_V()
{
    for (int i = SetPlate; i < GridCheck-SetPlate; i++)
    {
        potential[i][SetPlate] = 1;
        temp_potential[i][SetPlate] = 1;
        potential[i][GridCheck-SetPlate-1] = -1;
        temp_potential[i][GridCheck-SetPlate-1] = -1;
    }
}

void update_V()
{
    delta_V = 0;
    for (int i = 1; i < GridCheck-1; i++)
    {
        for (int j = 1; j < GridCheck-1; j++)
        {
            if (j != SetPlate && j != GridCheck-SetPlate-1)
            {
                temp_potential[i][j] = (potential[i-1][j]+potential[i+1][j]+potential[i][j-1]+potential[i][j+1]) / 4;
                delta_V += sqrt((potential[i][j]-temp_potential[i][j]) * (potential[i][j]-temp_potential[i][j]));
            }
            else if (i < SetPlate || i >= GridCheck-SetPlate)
            {
                temp_potential[i][j] = (potential[i-1][j]+potential[i+1][j]+potential[i][j-1]+potential[i][j+1]) / 4;
                delta_V += sqrt((potential[i][j]-temp_potential[i][j]) * (potential[i][j]-temp_potential[i][j]));
            }
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
            if (j != SetPlate && j != GridCheck-SetPlate-1)
            {
                temp_potential[i][j] = (potential[i-1][j]+potential[i+1][j]+potential[i][j-1]+potential[i][j+1]) / 4;
                potential[i][j] = alpha * (temp_potential[i][j] - potential[i][j]) + potential[i][j];
                delta_V += alpha * sqrt((potential[i][j]-temp_potential[i][j]) * (potential[i][j]-temp_potential[i][j]));
            }
            else if (i < SetPlate || i >= GridCheck-SetPlate)
            {
                temp_potential[i][j] = (potential[i-1][j]+potential[i+1][j]+potential[i][j-1]+potential[i][j+1]) / 4;
                potential[i][j] = alpha * (temp_potential[i][j] - potential[i][j]) + potential[i][j];
                delta_V += alpha * sqrt((potential[i][j]-temp_potential[i][j]) * (potential[i][j]-temp_potential[i][j]));
            }
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
void mainPro()
{
	initialize_V();
    int N_iter = laplace_calculate();
    cout<<N_iter<<endl;
  /*  for (int i = 0; i < GridCheck; ++i)
    {
        for (int j = 0; j < GridCheck; ++j)
        {
            cout<<potential[i][j]<<"    ";
        }
        cout<<endl;
    }
*/}

void h1draw()
{
 
gStyle->SetPalette(51);
gStyle->SetOptStat(0);
gStyle->SetOptTitle(1);

    TCanvas *c1 = new TCanvas("c1","Potential Distribution",200,10,800,400);
    mainPro();
    TH2D *myMatrix =new TH2D("h1","Potential Distribution", GridCheck,-2,GridCheck,GridCheck,-2,GridCheck);
    for (int i = 0; i < GridCheck; ++i)
    {
        for (int j = 0; j < GridCheck; ++j)
        {
            myMatrix->SetBinContent(i,j,potential[i][j]);        
        }    
    }
    c1->Divide(2,1);
    c1->cd(1);
    myMatrix->SetTitle("Perspective plot of the potential");
    myMatrix->Draw("surf3");
    c1->cd(2);
    myMatrix->SetTitle("Equipotential lines");
    myMatrix->Draw("Cont1");

}