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
#include <TLine.h>

// working arrays
double* sfit = 0; // templates
double* dfit = 0; // data
// fit parameters
int nsfit = 0; // number of templates
int nxfit = 0; // number of bins/template
int ixfit = 0; // number of light template scale factors to fit

void fill(int ansfit, int anxfit, int aixfit, TH1** hmc, TH1* hdata)
{
  nsfit = ansfit;
  nxfit = anxfit;
  ixfit = aixfit;
  if (sfit) delete sfit;
  if (dfit) delete dfit;
  sfit = new double[nsfit*nxfit];
  dfit = new double[nxfit];
  for (int ix = 0; ix<nxfit; ++ix) { // template bins
    for (int is = 0; is<nsfit; ++is) { // templates (l,c,b)
      int i = is*nxfit+ix;
      sfit[i] = hmc[is]->GetBinContent(ix+1);
    }
    dfit[ix] = hdata->GetBinContent(ix+1);
  }
}

void xlog(int& np, double* deriv, double& f, double par[], int flag)
{
  double sum = 0;
  for (int ix = 0; ix<nxfit; ++ix) { // template bins
    double pred = 0;
    for (int is = 0; is<nsfit; ++is) { // templates (l,c,b)
      int i = is*nxfit+ix;
      double dpred = sfit[i]*par[is];
      if (is==0 && ix<ixfit) dpred *= par[nsfit+ix];
      pred += dpred;
    }
    if (pred>0) sum += pred - dfit[ix]*log(pred);
  }
  f = 2.*sum;
}

void IPunbin()
{
  bool save_plots = false;

  const int nfl = 3;
  const int comp = 2;
  const int ip = 2;
  const int df = 4;
  const int nsf = 2;

  const TString chip[ip] = {"IP2D","IP3D"};
  const TString bu[comp] = {"_logpbpu","_logpcpu"};
  const char chfl[nfl] = {'l','c','b'};
  const char ipflv = 6;

  int icol[nfl] = {2,3,4};

  int i = 0;

  // Load in Files

  TFile* ff = new TFile("/home/jcrosby/WorkPlace/Histograms/Root_Files/data_r21.root");
  TFile* ff_mc = new TFile("/home/jcrosby/WorkPlace/Histograms/Root_Files/mc_Pythia_r21.root");

  // Set # of bins 
  int nbins = 0;
  TH1* bin = (TH1*)ff->Get("/subTagger/jet_IP2D_logpbpu_l");
  nbins = bin->GetNbinsX();
  //int nbins = 7;

  //Set Upper and Lower Bounds

  double ub = 10;
  double lb = -10;

  //Start the Big Loop!

  for(int ips = 0; ips<ip; ++ips){
    TString ph = chip[ips];
    for(int com=0; com<comp; ++com){
      ph += bu[com];
	
      // Load in Data
	
      TString sh = "jet_" + ph + "_l";
      TString final_name = "/subTagger/" + sh;
      //cout << final_name << endl;
      TH1* hd_raw = (TH1*)ff->Get(final_name);
      sh += "_eq";
      TH1* hd = new TH1D(sh,"",nbins,lb,ub);
      for (int ix = 1; ix<=nbins; ++ix) hd->SetBinContent(ix,hd_raw->GetBinContent(ix));
      int nxreal = nbins;
      if (nxreal>nbins) for (int ix = nbins+1; ix<=nxreal; ++ix) hd->AddBinContent(nbins,hd_raw->GetBinContent(ix));

      //Legend Creation

      TLegend* l1 = new TLegend(0.6,0.6,0.72,0.9);
      l1->SetBorderSize(0);
      l1->SetFillStyle(0);
      l1->SetTextFont(42);
      l1->SetMargin(0.4);
      
      //total histogram used for ratio

      TH1* tot =  new TH1D("","",nbins,lb,ub);

      //Stacked Histogram Creation
      
      TH1* h[nfl];
      THStack* ths = new THStack("",ph+ "_PreFit");
      for(int flav = 0; flav<nfl; ++flav){
	TString sh = "jet_" + ph +"_"+ chfl[flav];
	TString final_name = "/subTagger/" + sh;
	//cout << final_name <<endl;
	TH1* mc_raw = (TH1*)ff_mc->Get(final_name);
	h[flav] = new TH1D(sh,"",nbins,lb,ub);
	for (int ix = 1; ix<=nbins; ++ix) h[flav]->SetBinContent(ix,mc_raw->GetBinContent(ix));
	int nyreal = mc_raw->GetNbinsX();
	if (nyreal>nbins) for (int ix = nbins+1; ix<=nyreal; ++ix) h[flav]->AddBinContent(nbins,mc_raw->GetBinContent(ix));
	h[flav]->SetFillColor(icol[flav]);
	ths->Add(h[flav]);
	tot->Add(h[flav]);
	TString sl = chfl[flav]; sl+= " jets";
	l1->AddEntry(h[flav],sl,"f");
      }
      l1->AddEntry(hd,"Data","pl");


      //prefit plots
    
      TString sc = "c"; sc += ips; sc += com;
      TCanvas* c1 = new TCanvas(sc,sc);
      //c1->SetCanvasSize(1000,900);
      //c1->SetWindowSize(1000,900);

      //upper plot
      TPad* pad1 = new TPad("pad1","pad1",0,0.3,1,1.0);
      pad1->SetBottomMargin(.05);
      //pad1->SetGridx();         // Vertical grid
      pad1->Draw();             // Draw the upper pad: pad1
      pad1->cd();               // pad1 becomes the current pad
      //pad1->SetLogy();
      
      //put stack w/ data on pad 1

      ths->Draw();
      //ths->GetXaxis()->SetTitle("Arbitrary Units");
      ths->GetYaxis()->SetTitle("Arbitrary Units");
      //ths->GetYaxis()->SetTitleSize(25);
      ths->SetMinimum(1e2);
      hd->SetMarkerStyle(8);
      hd->SetMarkerSize(1.2);
      hd->Draw("esame");
      l1->Draw();
      
      //create pad 2
      c1->cd();
      TPad* pad2 = new TPad("pad2","pad2",0,0,1,0.3);
      pad2->SetTopMargin(0.05);
      pad2->SetBottomMargin(0.3);
      //pad2->SetGridx(); // vertical grid
      pad2->Draw();
      pad2->cd(); 

      //put ratio plot on pad 2

      TH1D *scf = (TH1D*)hd->Clone("sf");
      scf->SetLineColor(kBlack);
      scf->Sumw2();
      scf->Divide(tot);
      scf->SetMarkerStyle(21);
      scf->SetMarkerStyle(8);
      scf->SetMarkerSize(0.8);
      scf->SetMinimum(0.6);
      scf->SetMaximum(1.5);
      scf->GetXaxis()->SetTitle("Arbitrary Units");
      scf->SetStats(0);
       
      // Y axis ratio plot settings
      scf->GetYaxis()->SetNdivisions(505);
      // scf->GetYaxis()->SetTitleSize(25);
      // scf->GetYaxis()->SetTitleFont(43);
      // scf->GetYaxis()->SetTitleOffset(1.35); 
      scf->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
      scf->GetYaxis()->SetLabelSize(15);
       
      // X axis ratio plot settings
      scf->GetXaxis()->SetTitleSize(15);
      scf->GetXaxis()->SetTitleFont(43);
      scf->GetXaxis()->SetTitleOffset(3);
      scf->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel 
      scf->GetXaxis()->SetLabelSize(15);
       
      scf->Draw("esame");
         
	TLine *line1 = new TLine(lb,1,ub,1);
	line1->SetLineColor(kBlack);
	line1->Draw();
  
      if (save_plots) {
	TString sp = chip[ips]; sp += bu[com]; sp += "_prefit.png";
	c1->Print(sp);
      }

      //do the fit

      fill(nfl, nbins, nsf, h, hd);

      const int npar = nfl+nsf;
      TMinuit* minuit = new TMinuit(npar);
      //minuit->SetPrintLevel(-1); // suppress minuit output
      minuit->SetFCN(xlog);
      for (int ipar = 0; ipar<npar; ++ipar) {
	TString sp;
	if (ipar<nfl) { sp = "SF_frac_"; sp += chfl[ipar]; }
	else { sp = "SFl_bin_"; sp += ipar-nfl; }
	minuit->DefineParameter(ipar,sp,1.,0.1,0,0);
      }
       
      minuit->Migrad();
       
      double par[npar]; // fit parameters
      double err[npar]; // fit parameter errors
       
      for (int ipar = 0; ipar<npar; ++ipar) {
	minuit->GetParameter(ipar,par[ipar],err[ipar]);
      }
       
      // prepare scaled histograms
       
      TLegend* l2 = new TLegend(0.6,0.6,0.72,0.9);
      l2->SetBorderSize(0);
      l2->SetFillStyle(0);
      l2->SetTextFont(42);
      l2->SetMargin(0.4);
       
      // total scaled mc templates for ratio
       
      TH1* totf = new TH1D("","",nbins,lb,ub);
       
      double* sff= par; // flavor fraction scale factors
      double* sfl = par+nfl; // light scale factors in sensible bins
       
      cout << "0: " << sff[0] << endl;
      cout << "1: " << sff[1] << endl;
      cout << "2: " << sff[2] << endl;

      TH1* hs[nfl];
      THStack* thss = new THStack("", ph + "_PostFit");
      for (int jfl = 0; jfl<nfl; ++jfl) {
	sh = h[jfl]->GetName();// sh += "_scaled";
	hs[jfl] = (TH1*)h[jfl]->Clone(sh);
	for (int ix = 1; ix<=nbins; ++ix) {
	  double val = h[jfl]->GetBinContent(ix);
	  val *= sff[jfl];
	  if (ix<=nsf) val *= sfl[ix-1];
	  hs[jfl]->SetBinContent(ix,val);
	}
	hs[jfl]->SetFillColor(icol[jfl]);
	thss->Add(hs[jfl]);
	totf->Add(hs[jfl]);
	TString sl = chfl[jfl]; sl += " jets";
	l2->AddEntry(hs[jfl],sl,"f");
      }
      l2->AddEntry(hd,"Data","pl");



      // postfit plots
       
      TString sf = "c2"; sf += ips; sf += com;
      TCanvas* c2 = new TCanvas(sf,sf);
      //c2->SetCanvasSize(1000,900);
      //c2->SetWindowSize(1000,900);

      //upper plot
      TPad* pad3 = new TPad("pad3","pad3",0,0.3,1,1.0);
      pad3->SetBottomMargin(.05);
      //pad3->SetGridx();         // Vertical grid
      pad3->Draw();             // Draw the upper pad: pad1
      pad3->cd();               // pad1 becomes the current pad
      // pad3->SetLogy();
      
      //put stack w/ data on pad 3

      thss->Draw();
      //thss->GetXaxis()->SetTitle("Arbitrary Units");
      thss->GetYaxis()->SetTitle("Arbitrary Units");
      //thss->GetYaxis()->SetTitleSize(25);
      thss->SetMinimum(1e2);
      hd->SetMarkerStyle(8);
      hd->SetMarkerSize(1.2);
      hd->Draw("esame");
      l2->Draw();
       
      //create pad 4
      c2->cd();
      TPad* pad4 = new TPad("pad4","pad4",0,0,1,0.3);
      pad4->SetTopMargin(0.05);
      pad4->SetBottomMargin(0.3);
      //pad4->SetGridx(); // vertical grid
      pad4->Draw();
      pad4->cd(); 
       
      //put ratio plot on pad 4

      TH1D *scff = (TH1D*)hd->Clone("sff");
      scff->SetLineColor(kBlack);
      scff->Sumw2();
      scff->Divide(totf);
      scff->SetMarkerStyle(21);
      scff->SetMarkerStyle(8);
      scff->SetMarkerSize(0.8);
      scff->SetMinimum(0.6);
      scff->SetMaximum(1.5);
      scff->GetXaxis()->SetTitle("Arbitrary Units");
      scff->SetStats(0);
       
      // Y axis ratio plot settings
      scff->GetYaxis()->SetNdivisions(505);
      // scff->GetYaxis()->SetTitleSize(25);
      // scff->GetYaxis()->SetTitleFont(43);
      // scff->GetYaxis()->SetTitleOffset(1.35); 
      scff->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
      scff->GetYaxis()->SetLabelSize(15);
       
      // X axis ratio plot settings
      scff->GetXaxis()->SetTitleSize(15);
      scff->GetXaxis()->SetTitleFont(43);
      scff->GetXaxis()->SetTitleOffset(3);
      scff->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel 
      scff->GetXaxis()->SetLabelSize(15);
       
      scff->Draw("esame");
            
	TLine *line2 = new TLine(lb,1,ub,1);
	line2->SetLineColor(kBlack);
	line2->Draw();
	c2->Modified();
      
      if (save_plots) {
	TString sp = chip[ips]; sp += bu[com]; sp += "_postfit.png";
	c2->Print(sp);
      }
       
       
       
      ph = chip[ips];
       
    }
  }
}
