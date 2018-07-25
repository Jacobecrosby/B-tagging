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



void SubtaggerCompare(){

  bool save_plots = true;

  //Load Files

  TFile* mc = new TFile("/home/jcrosby/WorkPlace/Histograms/Root_Files/mc_Pythia_r21.root");

  //Create a legend

  auto ld1 = new TLegend(0.6,0.6,0.9,0.9);
  ld1->SetHeader("Legend","C");

  //Obtain Histograms from file

  TH1* h_IP2D_bu_l = (TH1*)mc->Get("/subTagger/jet_IP2D_logpbpu_l");
  TH1* h_IP2D_bu_b = (TH1*)mc->Get("/subTagger/jet_IP2D_logpbpu_b");

  TH1* h_IP3D_bu_l = (TH1*)mc->Get("/subTagger/jet_IP3D_logpbpu_l");
  TH1* h_IP3D_bu_b = (TH1*)mc->Get("/subTagger/jet_IP3D_logpbpu_b");

  //Get number of bins

  int nbinsx2dl = h_IP2D_bu_l->GetNbinsX();
  int nbinsx2db = h_IP2D_bu_b->GetNbinsX();
  int nbinsy2d = h_IP2D_bu_l->GetNbinsY();
  int nbinsx3dl = h_IP3D_bu_l->GetNbinsX();
  int nbinsx3db = h_IP3D_bu_b->GetNbinsX();
  int nbinsy3d = h_IP3D_bu_l->GetNbinsY();

  //Setup arrays to store histogram data

  vector<double> v_ip2d_l;
  vector<double> v_ip2d_b;

  vector<double> v_ip3d_l;
  vector<double> v_ip3d_b;

  vector<double> v_inv_ip2d_l;
  vector<double> v_inv_ip3d_l;
  double binval = 0;


  for (int lx = 0; lx <=nbinsx2dl; ++lx){
    binval = h_IP2D_bu_l->GetBinContent(lx);
    v_ip2d_l.push_back(binval);
    //cout << "IP2D_u: " << v_ip2d_l[lx] << endl;
    //ip2d_l[lx] = binval;
  }


 for (int lx = 0; lx <=nbinsx2dl; ++lx){
    binval = h_IP2D_bu_b->GetBinContent(lx);
    v_ip2d_b.push_back(binval);
    //cout << "IP2D_b: " << v_ip2d_b[lx] << endl;
    //ip2d_l[lx] = binval;
  }

 for (int lx = 0; lx <=nbinsx3dl; ++lx){
    binval = h_IP3D_bu_l->GetBinContent(lx);
    v_ip3d_l.push_back(binval);
    // cout << "IP3D_u: " << v_ip3d_l[lx] << endl;
    //ip2d_l[lx] = binval;
  }

 for (int lx = 0; lx <=nbinsx3db; ++lx){
    binval = h_IP3D_bu_b->GetBinContent(lx);
    v_ip3d_b.push_back(binval);
    //  cout << "IP3D_b: " << v_ip3d_b[lx] << endl;
    //ip2d_l[lx] = binval;
  }



 //Make inverse vectors

 float inv = 0;
 for(int iv =0; iv<= nbinsx2dl; ++iv){
   float x = v_ip2d_l.at(iv);
   inv = 1/x;
   v_inv_ip2d_l.push_back(inv);
 }

 for(int iv =0; iv<= nbinsx3dl; ++iv){
   float y = v_ip3d_l.at(iv);
   y = v_ip3d_l.at(iv);
   inv = 1/y;
   v_inv_ip3d_l.push_back(inv);
 }




 //Create 2D Histograms



 TH2* h1 = new TH2D("IP2D mc b-jet vs light jet", "IP2D bvu", nbinsx2dl, -50,100, nbinsy2d , -50,50);
 TH2* h2 = new TH2D("IP3D mc b-jet vs light jet", "IP2D bvu", nbinsx3dl, -50,100, nbinsy3d, -50, 50);

 //Fill Histograms

 float px, py;
 for(int ix = 0; ix<=nbinsx2dl; ++ix){
   px = v_ip2d_b.at(ix);
   py = v_inv_ip2d_l.at(ix);
   h1->Fill(px,py);

 }

 for(int ix = 0; ix<=nbinsx3dl; ++ix){
   px = v_ip3d_b.at(ix);
   py = v_inv_ip3d_l.at(ix);
   h2->Fill(px,py);

 }

 //Color palette

 int palette[5]={622, 623,625,633,635};
 gStyle->SetPalette(5,palette);

 //Make Canvas and draw histograms

 TCanvas* c1 = new TCanvas("c1","", 10,10,800,800);
 //c1->SetLogy();
 h1->GetXaxis()->SetTitle("b-jet");
 h1->GetYaxis()->SetTitle("light jet");
 h1->Draw("COLZ");
 c1->Modified();
 c1->Print("c1.png");


 TCanvas* c2 = new TCanvas("c2","",10,10,800,800);
 //c2->SetLogy();
 h2->GetXaxis()->SetTitle("b-jet");
 h2->GetYaxis()->SetTitle("light jet");
 h2->Draw("COLZ");
 c2->Modified();
 c2->Print("c2.png");

 if (save_plots){
  
   c1->SaveAs("/home/jcrosby/Histograms/");
   c2->SaveAs("/home/jcrosby/Histograms/");
 }


   // TString inputDir = "/home/jcrosby/Histograms/


}
