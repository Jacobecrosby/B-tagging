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



void Effrebin(){
  const int nfl = 3;
  const int comp = 2;
  const int ip = 2;
  const int df = 4;
  const int bf = 2;

  const TString chip[ip] = {"IP2D","IP3D"};
  const TString bu[comp] = {"_logpbpu","_logpcpu"};
  const char chfl[nfl] = {'l','c','b'};

  //Load in files

  TFile* ff = new TFile("/home/jcrosby/WorkPlace/Histograms/Root_Files/data_r21.root");
  TFile* ff_mc = new TFile("/home/jcrosby/WorkPlace/Histograms/Root_Files/mc_Pythia_r21.root");
  
  // Number of bins NOTE: All mc and data histograms contain 70 bins

  TH1* bin = (TH1*)ff_mc->Get("/subTagger/jet_IP2D_logpbpu_b");
  const int nbins = bin->GetNbinsX();
  const int effbins = 7;
   
  //Create the new TFile 

  TString inputDir = "/home/jcrosby/WorkPlace/Histograms/Root_Files/";
  TFile f(inputDir + "DT_rebinned_mcp_d.root","create");
  
  //Create my b-jet histograms

  TH1* bjr[bf];
  TH1* b_raw = 0;
  for(int ips =0;ips <ip; ++ips){
    TString ph = chip[ips];
    ph += "_logpbpu";
    TString sh = "jet_" + ph + "_b";
    TString final_name = "/subTagger/" + sh;
    TH1* b_raw = (TH1*)ff_mc->Get(final_name);
    sh += "_eq";
    cout << sh << endl;
    bjr[ips] = new TH1D(sh,"",nbins,0,nbins);
    for (int ix = 1; ix<=nbins; ++ix) bjr[ips]->SetBinContent(ix,b_raw->GetBinContent(ix));
    int nxreal = nbins;
    if (nxreal>nbins) for (int ix = nbins+1; ix<=nxreal; ++ix) bjr[ips]->AddBinContent(nbins,b_raw->GetBinContent(ix));
  }

  vector<double> bcv[bf];
  double bc = 0;
  double effn = 0;
  for(int bfd = 0; bfd < bf; ++bfd){
    for(int i = 0; i<=nbins; ++i){
      effn = bjr[bfd]->Integral(i,nbins)/bjr[bfd]->Integral(0,nbins);

      if(effn < 0.0001 ){
	bc = bjr[bfd]->GetBinCenter(i);
	if(bc == 69.5){
	  bcv[bfd].push_back(70);
	  cout << effn << ": " << bc << endl;
	}
      }
     
      else if(effn>.28 && effn<.32){
	bc = bjr[bfd]->GetBinCenter(i);
	bcv[bfd].push_back(bc);
	cout << effn << ": " << bc << endl;
      }
      else if(effn > 0.47 && effn < 0.53){
	bc = bjr[bfd]->GetBinCenter(i);
	bcv[bfd].push_back(bc);
	cout << effn << ": " << bc << endl;
      }
      else if (effn>.57 && effn < .63){
	bc = bjr[bfd]->GetBinCenter(i);
	bcv[bfd].push_back(bc);
	cout << effn << ": " << bc << endl;
      }
      else if (effn>.67 && effn < .73){
	bc = bjr[bfd]->GetBinCenter(i);
	bcv[bfd].push_back(bc);
	cout << effn << ": " << bc << endl;
      }
      else if (effn>.74 && effn < .79){
	bc = bjr[bfd]->GetBinCenter(i);
	bcv[bfd].push_back(bc);
	cout << effn << ": " << bc << endl;
      }
      else if (effn>.83 && effn < .87){
	bc = bjr[bfd]->GetBinCenter(i);
	bcv[bfd].push_back(bc);
	cout << effn << ": " << bc << endl;
      }
      else if (effn == 1){
	bc = bjr[bfd]->GetBinCenter(i);
	if(bc == 0.5){
	  bcv[bfd].push_back(bc);
	  cout << effn << ": " << bc << endl;
	}
      }
    }
  }
  
  TH1* mcf[12];
  TH1* datf[df];

  // Rebin the data and MC templates in order to apply cuts
  
  for(int ips = 0; ips < ip; ++ips){
    TString ph = chip[ips];
    for(int com = 0; com < comp; ++ com){
      ph += bu[com];
      
      //Load in data
      
      TString sh = "jet_" + ph + "_l";
      TString final_name = "/subTagger/" + sh;
      TH1* hd_raw = (TH1*)ff->Get(final_name);
      sh += "_d";
      cout << sh << endl;
      datf[ips*2+com] = new TH1D(sh,"", effbins, 0, effbins);
      
      
      //fill efficiency bins with data

      for(int ix = 1; ix <=effbins; ++ix){
	float check = hd_raw->GetNbinsX();
	//cout << "Data bins: " << check << endl;
	  
	Double_t up_integral = hd_raw->Integral(bcv[ips].at(ix-1)-0.5, nbins);
	float low_integral = hd_raw->Integral(bcv[ips].at(ix)-0.5, nbins);
	float wanted_val = up_integral - low_integral;
	//cout << "Data value: " << wanted_val << endl;
	  
	datf[ips*2+com]->SetBinContent(ix,wanted_val);
      }
      int nxreal = hd_raw->GetNbinsX();
      if (nxreal>nbins) for (int ix = effbins+1; ix<=nxreal; ++ix) datf[ips*2+com]->AddBinContent(effbins,hd_raw->GetBinContent(ix));
	
	
      // Loop through Monte Carlo templates to rebin
	
      for(int flav = 0; flav < nfl; ++flav){
	TString wh = "jet_" + ph + "_" + chfl[flav];
	TString final_name = "/subTagger/" + wh;
	cout << final_name << endl;
	TH1* mc_raw = (TH1*)ff_mc->Get(final_name);
	mcf[flav+3*com+6*ips] = new TH1D(wh, "", effbins, 0, effbins);

	//Fill efficiency bins with MC templates
	  
	for(int ix = 1; ix <=effbins; ++ix){ 
	  float checkk = mc_raw->GetNbinsX();
	  // cout << "MC bins: " << checkk << endl;
	    
	  float upp_integral = mc_raw->Integral(bcv[ips].at(ix-1)-0.5,nbins);
	  float loww_integral = mc_raw->Integral(bcv[ips].at(ix)-0.5, nbins);
	  float wanted_vall = upp_integral - loww_integral;
	  //cout << "MC Value: " << wanted_vall << endl;
	    
	  mcf[flav+3*com+6*ips]->SetBinContent(ix,wanted_vall);
	}
	  
	int nxreal = mc_raw->GetNbinsX();
	if (nxreal>nbins) for (int ix = effbins+1; ix<=nxreal; ++ix) mcf[flav+3*com+6*ips]->AddBinContent(effbins,mc_raw->GetBinContent(ix));
	    
      }
	
      ph = chip[ips];

    }
  }

  //Write histograms to root file

  f.cd();

  for(int dfw=0; dfw<df; ++dfw){
    datf[dfw]->SetStats(0);
    datf[dfw]->Write();
    //cout << dfw << endl;
  }

  for(int mcw = 0; mcw<=11; ++mcw){
    mcf[mcw]->SetStats(0);
    mcf[mcw]->Write();
    // cout << mcw << endl;
  }
  
  //f.Write();
  f.Close();
  
}

