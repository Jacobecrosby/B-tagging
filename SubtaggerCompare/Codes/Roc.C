#include <iostream>
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TCanvas.h"
#include "TString.h"
#include "TStyle.h"
#include "THStack.h"
#include "TMinuit.h"
 





void Roc()
{ 

  bool save_plots = true;

  const int nfl = 3; 
  const char chfl[nfl] = {'l','c','b'};

  // # of bins
  const int nx = 5000;

  // Load MC templates

  TFile* ff_mc = new TFile("/home/jcrosby/WorkPlace/Histograms/Root_Files/mc_Pythia_r21.root");


  
  vector<double> v[nfl];
  double eff = 0;


  // Obtain efficiencies of l, b, c

  TH1* h[nfl];
  for(int fl = 0; fl<nfl; ++fl){
    TString sh = "jet_MV2c10_";
    sh += chfl[fl];
    TString final_name = "/subTagger/" + sh;
    TH1* mc = (TH1*)ff_mc->Get(final_name);
    h[fl] = new TH1D(sh, "", nx,0,nx); 
    for (int ix = 1; ix<=nx; ++ix) h[fl]->SetBinContent(ix,mc->GetBinContent(ix));
    int nyreal = mc->GetNbinsX();
    if (nyreal>nx) for (int ix = nx+1; ix<=nyreal; ++ix) h[fl]->AddBinContent(nx,mc->GetBinContent(ix));      
    for(int efb = 0; efb <= nx ; ++efb){
      eff = h[fl]->Integral(0,efb)/h[fl]->Integral(0,nx);
      v[fl].push_back(eff);
      // cout << eff << endl;
    }
  }
  

  // Create light jet rejection vector

  vector<double> l_rej; 
  double lr = 0;

  for(int i = 0; i<= nx; ++i){
    if(v[0].at(i)!=0){
      lr = 1/v[0].at(i);
      l_rej.push_back(lr);
      //cout << lr << endl;
    }
    else {
      cout << "Divided by 0" << endl;
      l_rej.push_back(0);
      }
  }


  TCanvas* c1 = new TCanvas("c1","", 10,10,800,800);;
  TGraph* gr = new TGraph(nx);
  for(uint i=0; i<l_rej.size(); ++i){
    if(l_rej[i]<1000){
      gr->SetPoint(i, v[2].at(i), l_rej[i]);
    }
  }
  gr->SetTitle("B-jet Efficiency vs Light Jet Rejection");
  gr->Draw("A*");
  c1->SetLogy();
  c1->Modified();

  if(save_plots){
    c1->SaveAs("Roc_beff_lrej.png");
  }
  
}  
