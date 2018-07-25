
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

void fit_temp()
{
  bool save_plots = true;

  const int nfl = 3;
  const char chfl[nfl] = {'l','c','b'};
  int icol[nfl] = {2,3,4};

  // number of light template bins to fit
  const int nsf = 2;

  TFile* ff = new TFile("/home/jcrosby/Histograms/Root_Files/data_r21.root");
  TFile* ff_mc = new TFile("/home/jcrosby/Histograms/Root_Files/mc_Pythia_r21.root");
  //TFile* ff_mc = new TFile("data/mc_HERWIG_r21.root");
  //TFile* ff_mc = new TFile("data/mc_Shrpa_r21.root");

  //TFile* ff = new TFile("data/data_r20.7.root");
  //TFile* ff_mc = new TFile("data/mc_Pythia_r20.7.root");



  //Retrieving information from files
  const int npt = 8;
  const int neta = 2;
  double xpt[npt+1] = {20, 60, 100, 200, 300, 500, 750, 1000, 3000};
  const int nhsf = nfl-1+nsf; // c, b + light template bins
  TH1* hsf[neta*nhsf];
  TH1* hfr0[neta*nhsf];
  TH1* hfr1[neta*nhsf];
  for (int ieta = 1; ieta<=neta; ++ieta) {
    for (int ihsf = 0; ihsf<nhsf; ++ihsf) {
      int ih = (ieta-1)*nhsf+ihsf;
      //sh = looping to retrieve names of files
      TString sh = "sf"; sh += ieta; sh += ihsf;
      hsf[ih] = new TH1F(sh,"",npt,xpt);
      sh = "fr0"; sh += ieta; sh += ihsf;
      hfr0[ih] = new TH1F(sh,"",npt,xpt);
      sh = "fr1"; sh += ieta; sh += ihsf;
      hfr1[ih] = new TH1F(sh,"",npt,xpt);
    }
  }

  // number of efficiency bins
  const int nx = 7;

  for (int ipt = 1; ipt<=npt; ++ipt) {
    for (int ieta = 1; ieta<=neta; ++ieta) {
  //{ int ipt = 3; { int ieta = 1;
      TString sh = "w_";
      sh += "_pt"; sh += ipt;
      sh += "eta"; sh += ieta;
      sh += "_MV2c10";
      TH1* hd_raw = (TH1*)ff->Get(sh);
      sh += "_eq";
      TH1* hd = new TH1D(sh,"",nx,0,nx);
      for (int ix = 1; ix<=nx; ++ix) hd->SetBinContent(ix,hd_raw->GetBinContent(ix));
      int nxreal = hd_raw->GetNbinsX();
      if (nxreal>nx) for (int ix = nx+1; ix<=nxreal; ++ix) hd->AddBinContent(nx,hd_raw->GetBinContent(ix));


      //Legend creation


      TLegend* l1 = new TLegend(0.8,0.6,0.92,0.9);
      l1->SetBorderSize(0);
      l1->SetFillStyle(0);
      l1->SetTextFont(42);
      l1->SetMargin(0.4);


      //Histogram Creation


      TH1* h[nfl];
      THStack* ths = new THStack();
      for (int jfl = 0; jfl<nfl; ++jfl) {
	sh = "w_"; sh += chfl[jfl];
	sh += "_pt"; sh += ipt;
	sh += "eta"; sh += ieta;
	sh += "_MV2c10";
	TH1* h_raw = (TH1*)ff_mc->Get(sh);
	sh += "_eq";
	h[jfl] = new TH1D(sh,"",nx,0,nx);
	for (int ix = 1; ix<=nx; ++ix) h[jfl]->SetBinContent(ix,h_raw->GetBinContent(ix)); 
	int nxreal = h_raw->GetNbinsX();
	if (nxreal>nx) for (int ix = nx+1; ix<=nxreal; ++ix) h[jfl]->AddBinContent(nx,h_raw->GetBinContent(ix));
	h[jfl]->SetFillColor(icol[jfl]);
	ths->Add(h[jfl]);
	TString sl = chfl[jfl]; sl += " jets";
	l1->AddEntry(h[jfl],sl,"f");
      }
      l1->AddEntry(hd,"Data","pl");

      // prefit plots

      TString sc = "c"; sc += ipt; sc += ieta;
      TCanvas* c1 = new TCanvas(sc,sc);
      c1->SetLogy();
      ths->Draw();
      ths->GetXaxis()->SetTitle("Template bins");
      ths->GetYaxis()->SetTitle("Arbitrary units");
      ths->SetMinimum(1e2);
      hd->SetMarkerStyle(8);
      hd->SetMarkerSize(1.2);
      hd->Draw("esame");
      l1->Draw();
      c1->Modified();

      if (save_plots) {
	TString sp = "pt"; sp += ipt; sp += "eta"; sp += ieta; sp += "_prefit.pdf";
	c1->Print(sp);
      }

      // do the fit

      fill(nfl, nx, nsf, h, hd);

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

      TLegend* l2 = new TLegend(0.8,0.6,0.92,0.9);
      l2->SetBorderSize(0);
      l2->SetFillStyle(0);
      l2->SetTextFont(42);
      l2->SetMargin(0.4);

      double* sff= par; // flavor fraction scale factors
      double* sfl = par+nfl; // light scale factors in sensible bins
      TH1* hs[nfl];
      THStack* thss = new THStack();
      for (int jfl = 0; jfl<nfl; ++jfl) {
	sh = h[jfl]->GetName(); sh += "_scaled";
	hs[jfl] = (TH1*)h[jfl]->Clone(sh);
	for (int ix = 1; ix<=nx; ++ix) {
	  double val = h[jfl]->GetBinContent(ix);
	  val *= sff[jfl];
	  if (ix<=nsf) val *= sfl[ix-1];
	  hs[jfl]->SetBinContent(ix,val);
	}
	hs[jfl]->SetFillColor(icol[jfl]);
	thss->Add(hs[jfl]);
	TString sl = chfl[jfl]; sl += " jets";
	l2->AddEntry(hs[jfl],sl,"f");
      }
      l2->AddEntry(hd,"Data","pl");

      // postfit plots

      sc = "c"; sc += ipt; sc += ieta; sc += "s";
      TCanvas* c2 = new TCanvas(sc,sc);
      c2->SetLogy();
      thss->Draw();
      thss->SetMinimum(1e2);
      thss->GetXaxis()->SetTitle("Template bins");
      thss->GetYaxis()->SetTitle("Arbitrary units");
      hd->SetMarkerStyle(8);
      hd->SetMarkerSize(1.2);
      hd->Draw("esame");
      l2->Draw();
      c2->Modified();

      if (save_plots) {
	TString sp = "pt"; sp += ipt; sp += "eta"; sp += ieta; sp += "_postfit.pdf";
	c2->Print(sp);
      }

      // keep the results

      for (int ipar = 1; ipar<npar; ++ipar) {
	int ih = (ieta-1)*nhsf+ipar-1;
	hsf[ih]->SetBinContent(ipt,par[ipar]);
	hsf[ih]->SetBinError(ipt,err[ipar]);
	double val = 0;
	if (ipar<=2) { // c/b-frac
	  val = h[ipar]->Integral()/(h[0]->Integral()+h[1]->Integral()+h[2]->Integral());
	} else { // light tag rate
	  val = h[0]->GetBinContent(ipar-2)/h[0]->Integral();
	}
	hfr0[ih]->SetBinContent(ipt,val);
	hfr1[ih]->SetBinContent(ipt,val*par[ipar]);
	hfr1[ih]->SetBinError(ipt,val*err[ipar]);
      }
    }
  }

  //Histogram Canvas titles/creation

  for (int ieta = 1; ieta<=neta; ++ieta) {
    for (int ihsf = 0; ihsf<nhsf; ++ihsf) {
      int ih = (ieta-1)*nhsf+ihsf;

      TString sc = "c"; sc += ieta; sc += ihsf; sc += "sf";
      TCanvas* c3 = new TCanvas(sc,sc);
      c3->SetGridx();
      c3->SetGridy();
      c3->SetLogx();
      hsf[ih]->SetMarkerStyle(8);
      hsf[ih]->SetMarkerSize(1.2);
      hsf[ih]->SetMinimum(0.2);
      hsf[ih]->SetMaximum(2.2);
      hsf[ih]->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
      TString sy;
      if (ihsf==0) sy = "c-scale factor";
      else if (ihsf==1) sy = "b-scale factor";
      else { sy = "l"; sy += ihsf-2; sy += "-scale factor"; }
      hsf[ih]->GetYaxis()->SetTitle(sy);
      hsf[ih]->Draw("e0");
      c3->Modified();
      c3->SaveAs("/home/jcrosby/Histograms/");

      sc = "c"; sc += ieta; sc += ihsf; sc += "fr";
      TCanvas* c4 = new TCanvas(sc,sc);
      c4->SetGridx();
      c4->SetGridy();
      c4->SetLogx();
      hfr1[ih]->SetMarkerStyle(8);
      hfr1[ih]->SetMarkerSize(1.2);
      hfr1[ih]->GetXaxis()->SetTitle("Jet p_{T} [GeV]");
      if (ihsf==0) sy = "c fraction";
      else if (ihsf==1) sy = "b fraction";
      else { sy = "l"; sy += ihsf-2; sy += " tag rate"; }
      hfr1[ih]->GetYaxis()->SetTitle(sy);
      hfr1[ih]->Draw("e0");

      hfr0[ih]->SetLineColor(2);
      hfr0[ih]->SetLineWidth(2);
      hfr0[ih]->Draw("same");

      if (hfr1[ih]->GetMinimum()>hfr0[ih]->GetMinimum()) hfr1[ih]->SetMinimum(0.9*hfr0[ih]->GetMinimum()); 
      if (hfr1[ih]->GetMaximum()<hfr0[ih]->GetMaximum()) hfr1[ih]->SetMaximum(1.1*hfr0[ih]->GetMaximum()); 
      c4->Modified();

      if (save_plots) {
	TString sp = "SF_";
	if (ihsf==0) sp += "c";
	else if (ihsf==1) sp += "b";
	else { sp += "l"; sp += ihsf-2; }
	sp += "_eta"; sp += ieta;
	sp += ".pdf";
	c3->Print(sp);

	sp = "rate_";
	if (ihsf==0) sp += "c";
	else if (ihsf==1) sp += "b";
	else { sp += "l"; sp += ihsf-2; }
	sp += "_eta"; sp += ieta;
	sp += ".pdf";
	c4->Print(sp);
      }
    }
  }
}
