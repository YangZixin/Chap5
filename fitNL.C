#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <fstream>
using namespace std;

void fitNL()
{
	int Leng[] = {11,21,31,41,51,61,71,81,91,101};
	int Niter1[] = {38,140,330,574,860,1180,1526,1895,2280,2679};
	int Niter2[] = {30,93,104,227,351,472,593,713,831,949};
	TCanvas *c1 = new TCanvas("c","fitNL",0,0,500,400);
	TGraph *gr1 = new TGraph(10, Leng, Niter1);
	TGraph *gr2 = new TGraph(10, Leng, Niter2);
	gr1->GetXaxis()->SetTitle("Length of the Capacitor Plate");
	gr1->GetYaxis()->SetTitle("N_iter");
	gr1->SetTitle("N_iter ~ L");
	gr1->SetMarkerStyle(8);
	gr1->SetMarkerColor(3);
	gr1->SetLineColor(3);
	gr1->Draw("AP");

	gr2->SetMarkerStyle(8);
	gr2->SetMarkerColor(9);
	gr2->SetLineColor(9);
	gr2->Draw("P");
  	TLegend *legend = new TLegend(0.16,0.77,0.36,0.86);
   	legend->AddEntry(gr1,"Jacobi Method","L");
    legend->AddEntry(gr2,"SOR algorithm","L");
   	legend->Draw();

   	gr1->Fit("pol2");
   	gr2->Fit("pol1");




}