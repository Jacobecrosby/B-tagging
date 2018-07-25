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


//To Do: Makes 3 Histograms from jet_pt

//One: x < 200 GeV -> 200 < x < 500 GeV -> x > 500 GeV






void MyHistoMaker()
{
  
  //Global Variables
  
  
  Double_t nbins1 = 1000;
  Double_t nbins2 = 500;
  Double_t nbins3 = 500;
  //  Double_t cut = 1;
  
  
  //Define the 1D jet_pt histograms. Varying nbins and range for demonstration
  
  
  auto jet_pt1a = new TH1D("jet_pt1a","Jet pt < 50 GeV;X axis;Y axis",nbins1,0,50e3);
  auto jet_pt1b = new TH1D("jet_pt1b","Jet pt < 100 Gev",nbins1,0,100e3);
  auto jet_pt1c = new TH1D("jet_pt1c","Jet pt < 200 GeV",nbins1,0,200e3);
  
  auto jet_pt2a = new TH1D("jet_pt2a","Jet pt < 300 GeV",nbins2,0,300e3);
  auto jet_pt2b = new TH1D("jet_pt2b","Jet pt < 400 GeV",nbins2,0,400e3);
  auto jet_pt2c = new TH1D("jet_pt2c","Jet pt < 500 GeV",nbins2,0,500e3);
  
  auto jet_pt3a = new TH1D("jet_pt3a","Jet pt < 600 GeV",nbins3,0,600e3);
  auto jet_pt3b = new TH1D("jet_pt3b","Jet pt < 700 GeV",nbins3,0,700e3);
  auto jet_pt3c = new TH1D("jet_pt3c","Jet pt < 800 GeV",nbins3,0,800e3);



  //Few Coding tips for styles:

  // Color and line  codes are here : https://root.cern.ch/doc/master/classTAttLine.html
  
  //TH1 Class Reference: https://root.cern.ch/doc/master/classTH1.html

  // jet_pt1a->GetXaxis()->SetTitle("x axis title");
  // jet_pt1a->SetLineColor(2);
  // jet_pt1a->SetLineWidth(3);
  // jet_pt1a->SetLineStyle(7);
  // jet_pt2c->SetFillColor(23);


  //retrieve file from
  
  TString inputDir = "/home/jcrosby/Histograms/Root_Files/";
  
  
  
  
  //Open and access to the file
  
  
  auto myFile = TFile::Open(inputDir+"sample_ttbar.root");
  
  
  TTreeReader myReader("bTag_AntiKt4EMTopoJets",myFile);
  
  
  
  //Fetch the variables from the jet_pt file
  
  TTreeReaderValue<std::vector<float>> jet_pt(myReader,"jet_pt");
  
  
  //loop over jet_pt
  
  while(myReader.Next()){
    
    //y axis (= 1 due to normalization). Used within "Fill"
    
    double y = 1;
    
    
    
    for(uint i = 0; i< (*jet_pt).size(); i++)
      {
	
	jet_pt1a->Fill((*jet_pt).at(i),y);
	jet_pt1b->Fill((*jet_pt).at(i),y);
	jet_pt1c->Fill((*jet_pt).at(i),y);

	jet_pt2a->Fill((*jet_pt).at(i),y);
	jet_pt2b->Fill((*jet_pt).at(i),y);
	jet_pt2c->Fill((*jet_pt).at(i),y);

	jet_pt3a->Fill((*jet_pt).at(i),y);
	jet_pt3b->Fill((*jet_pt).at(i),y);
	jet_pt3c->Fill((*jet_pt).at(i),y);

      }//end of "for" loop


  }//end of "while" loop


  jet_pt1a->Sumw2();
  jet_pt1b->Sumw2();

  //Write histograms to root file and save into a directory

  TFile f(inputDir+"sample_jet_pt"+"_histo.root","recreate");

  jet_pt1a->Write();
  jet_pt1b->Write();
  jet_pt1c->Write();

  jet_pt2a->Write();
  jet_pt2b->Write();
  jet_pt2c->Write();

  jet_pt3a->Write();
  jet_pt3b->Write();
  jet_pt3c->Write();

  f.Write();
  f.Close();


}//end of code
