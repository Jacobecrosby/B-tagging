#include <vector>
#include <iostream>
#include <sstream>
#include <map>
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TLegend.h"
#include "TROOT.h"
#include "TTreeReader.h"
#include "TTreeReaderValue.h"
#include "TTreeReaderArray.h"
#include "TTree.h"
#include "TBranch.h"
#include "TMath.h"
#include "TImage.h"
#include "TCanvas.h"
#include "TString.h"
#include "TAttFill.h"
#include "TStyle.h"
#include "THStack.h"
#include "TLegend.h"
#include "TLatex.h"
#include <math.h>
#include <cmath>
#include "TMinuit.h"



void smtrebin(){
  const int nfl = 3;
  const int df = 4;
  const int bf = 2;
  const int effpar = 8;

  const char chfl[nfl] = {'l','c','b'};

  //Load in files

  TFile* ff = new TFile("/home/jcrosby/WorkPlace/Histograms/Root_Files/data.root");
  TFile* ff_mc = new TFile("/home/jcrosby/WorkPlace/Histograms/Root_Files/mc.root");
  
  // Number of bins NOTE: All mc and data histograms contain 200 bins

  TString sh = "jet_mu_pTrel_b";
  TH1* bin = (TH1*)ff_mc->Get(sh);
  const int nbins = bin->GetNbinsX();
  const int effbins = 7;
   
  //Create the new TFile 

  TString inputDir = "/home/jcrosby/WorkPlace/Histograms/Root_Files/";
  TFile f(inputDir + "mu_pTrel_rebin.root","recreate");

  //Create my b-jet histograms
  
  TString final_name = sh;
  TH1* b_raw = (TH1*)ff_mc->Get(final_name);
  cout << sh << endl;
  sh += "_eq";
  TH1* hd = new TH1D(sh,"",nbins,0,nbins);
  for (int ix = 1; ix<=nbins; ++ix) hd->SetBinContent(ix,b_raw->GetBinContent(ix));
  int nxreal = nbins;
  if (nxreal>nbins) for (int ix = nbins+1; ix<=nxreal; ++ix) hd->AddBinContent(nbins,b_raw->GetBinContent(ix));
  

 float checkk = b_raw->GetNbinsX();
 cout << "MC bbins: " << checkk << endl;
    
  vector<double> bcv;
  double bc = 0;
  double effn = 0;
  for(int i = 0; i<=nbins; ++i){
    effn = hd->Integral(i,nbins)/hd->Integral(0,nbins);

    if(effn < 0.0001 ){
      bc = hd->GetBinCenter(i);
      if(bc == 199.5){
	bcv.push_back(200.5);
	cout << effn << ": " << bc << endl;
	}
    }
    else if(effn>.28 && effn<.32){
      bc = hd->GetBinCenter(i);
      bcv.push_back(bc);
      cout << effn << ": " << bc << endl;
    }
    else if(effn > 0.45 && effn < 0.5){
      bc = hd->GetBinCenter(i);
      bcv.push_back(bc);
      cout << effn << ": " << bc << endl;
    }
    else if (effn>.54 && effn < .64){
      bc = hd->GetBinCenter(i);
      bcv.push_back(bc);
      cout << effn << ": " << bc << endl;
    }
    else if (effn>.65 && effn < .7){
      bc = hd->GetBinCenter(i);
      bcv.push_back(bc);
      cout << effn << ": " << bc << endl;
    }
    else if (effn>.74 && effn < .79){
      bc = hd->GetBinCenter(i);
      bcv.push_back(bc);
      cout << effn << ": " << bc << endl;
    }
    else if (effn>.81 && effn < .89){
      bc = hd->GetBinCenter(i);
      bcv.push_back(bc);
      cout << effn << ": " << bc << endl;
    }
    else if (effn == 1){
      bc = hd->GetBinCenter(i);
      if(bc == 0.5){
	bcv.push_back(bc);
	cout << effn << ": " << bc << endl;
      }
    }
  }	

  for(int re = 0; re<=7; ++re){
    cout << bcv.at(re) << endl;
  }
				     


  TH1* mcf[4];
  

  // Rebin the data and MC templates in order to apply cuts
         
  //Load in data

  TString ph = "jet_mu_pTrel_l";
  TString finall_name = ph;
  TH1* hd_raw = (TH1*)ff->Get(finall_name);
  ph += "_d";
  cout << ph << endl;
  TH1* hs = new TH1D(ph,"", effbins, 0, effbins);
      
      
  //fill efficiency bins with data
  
  for(int ix = 1; ix <=effbins; ++ix){
    float check = hd_raw->GetNbinsX();
    //cout << "Data bins: " << check << endl;
    
    Double_t up_integral = hd_raw->Integral(bcv.at(ix-1)-0.5, nbins);
    float low_integral = hd_raw->Integral(bcv.at(ix)-0.5, nbins);
    float wanted_val = up_integral - low_integral;
    cout << "Data value: " << wanted_val << endl;
    
    hs->SetBinContent(ix,wanted_val);
  }
  int xreal = hd_raw->GetNbinsX();
  if (xreal>nbins) for (int ix = effbins+1; ix<=xreal; ++ix) hs->AddBinContent(effbins,hd_raw->GetBinContent(ix));
  
  
  // Loop through Monte Carlo templates to rebin
  
  for(int flav = 0; flav < nfl; ++flav){
    TString wh = "jet_mu_pTrel_";
    wh += chfl[flav];
    TString final_name = wh;
    cout << final_name << endl;
    TH1* mc_raw = (TH1*)ff_mc->Get(final_name);
    mcf[flav] = new TH1D(wh, "", effbins, 0, effbins);
    
    //Fill efficiency bins with MC templates
    
    for(int ix = 1; ix <=effbins; ++ix){ 
      float checkk = mc_raw->GetNbinsX();
      // cout << "MC bins: " << checkk << endl;
      
      float upp_integral = mc_raw->Integral(bcv.at(ix-1)-0.5,nbins);
      float loww_integral = mc_raw->Integral(bcv.at(ix)-0.5, nbins);
      float wanted_vall = upp_integral - loww_integral;
      cout << "MC Value: " << wanted_vall << endl;
      
      mcf[flav]->SetBinContent(ix,wanted_vall);
    }
    
    int nxreal = mc_raw->GetNbinsX();
    if (nxreal>nbins) for (int ix = effbins+1; ix<=nxreal; ++ix) mcf[flav]->AddBinContent(effbins,mc_raw->GetBinContent(ix));
    
  }
  
  
  //Write histograms to root file
  
  f.cd();
  
  hs->SetStats(0);
  hs->Write();
  
  

  for(int mcw = 0; mcw<nfl; ++mcw){
    mcf[mcw]->SetStats(0);
    mcf[mcw]->Write();
  }
  
  //f.Write();
  f.Close();
  
  
}


