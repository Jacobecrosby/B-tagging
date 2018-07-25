//#include <errno.h>
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

void SubtaggerFit()
{
  bool save_plots = false;
  
  const int nfl = 3;
  const char chfl[nfl] = {'l','c','b'};
  int icol[nfl] = {2,3,4};
  
  // number of light template bins to fit
  
  const int nsf = 2;
  
  TFile* ff = new TFile("/home/jcrosby/WorkPlace/Histograms/Root_Files/data_r21.root");
  TFile* ff_mc = new TFile("/home/jcrosby/WorkPlace/Histograms/Root_Files/mc_Pythia_r21.root");
  
  
  
  //change bins to original
  // const int nx = 7;
  
  //Set # of bins to tagger file
  TH1* mc_bins = (TH1*)ff_mc->Get("subTagger/jet_MV2c10_l");
  int nx = mc_bins->GetNbinsX();
  int i = 0;
  
  
  //retrive the data
  
  TH1* h_raw = (TH1*)ff->Get("/subTagger/jet_MV2c10_l");
  TH1* hs = new TH1D("data","",nx,0,nx);
  for (int ix = 1; ix<=nx; ++ix) hs->SetBinContent(ix,h_raw->GetBinContent(ix));
  int nxreal = h_raw->GetNbinsX();
  if (nxreal>nx) for (int ix = nx+1; ix<=nxreal; ++ix) hs->AddBinContent(nx,h_raw->GetBinContent(ix));
      
      

  //Legend creation

  
  TLegend* l1 = new TLegend(0.65,0.7,0.95,0.9);
  l1->SetBorderSize(0);
  l1->SetFillStyle(0);
  l1->SetTextFont(42);
  l1->SetMargin(0.4);

  //total histogram used for ratio

  TH1* tot= new TH1D("","",nx,0,nx);


  //Stacked Histogram Creation 

  TH1* h[nfl];
  THStack* ths = new THStack();
  for ( int flav = 0; flav < nfl; ++flav){	
    TString sh = "jet_MV2c10_";
    sh += chfl[flav];	
    TString final_name = "/subTagger/" + sh;
    TH1* hd_raw = (TH1*)ff_mc->Get(final_name);
    h[flav] = new TH1D(sh,"",nx,0,nx);
    for (int ix = 1; ix<=nx; ++ix) h[flav]->SetBinContent(ix,hd_raw->GetBinContent(ix));
    int nyreal = hd_raw->GetNbinsX();
    if (nyreal>nx) for (int ix = nx+1; ix<=nyreal; ++ix) h[flav]->AddBinContent(nx,hd_raw->GetBinContent(ix));
    tot->Add(h[flav]);
    h[flav]->SetFillColor(icol[flav]);
    ths->Add(h[flav]);
    TString sl = chfl[flav]; sl+= " jets";
    l1->AddEntry(h[flav],sl,"f");
  }

 
  l1->AddEntry(hs,"Data","pl");
  
  // prefit plot w/ ratio
  
  TString sc = "c1" + to_string(i);
  TCanvas* c1 = new TCanvas(sc,sc);
  c1->SetCanvasSize(1000,900);
  c1->SetWindowSize(1000,900);

  //upper plot
  TPad *pad1 = new TPad("pad1","pad1", 0,0.3,1,1.0);
  pad1->SetBottomMargin(0); // Upper and lower plot are joined
  pad1->SetGridx();         // Vertical grid
  pad1->Draw();             // Draw the upper pad: pad1
  pad1->cd();               // pad1 becomes the current pad
  pad1->SetLogy();

  ths->Draw();
  ths->GetXaxis()->SetTitle("Template bins");
  ths->GetYaxis()->SetTitle("Arbitrary units");
  // ths->GetYaxis()->SetTitleSize(13); //ERROR: GETS RID OF TITLE
  ths->GetYaxis()->SetTitleOffset(1.25);
  ths->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  ths->GetYaxis()->SetLabelSize(15);
  ths->SetTitle("MV2c10 MC vs l-jet Data {Pre-fit}");
  ths->SetMinimum(1e2);
  hs->SetMarkerStyle(8);
  hs->SetMarkerSize(1.2);
  hs->Draw("esame");
  l1->Draw();
  

  c1->cd();
  TPad *pad2 = new TPad("pad2","pad2", 0, .01, 1, 0.4);
  pad2->SetTopMargin(0.04);
  pad2->SetBottomMargin(0.4);
  pad2->SetGridx(); // vertical grid
  pad2->Draw();
  pad2->cd(); 
  //pad2->SetLogy(); 

  //pre-fit ratio plot on pad 2
  
  TH1D *scf = (TH1D*)hs->Clone("sf");
  scf->SetLineColor(kBlack);
  scf->Sumw2();
  scf->Divide(tot);
  scf->SetMarkerStyle(21);
  scf->SetMarkerStyle(8);
  scf->SetMarkerSize(0.8);
  scf->SetMinimum(0.8);
  scf->SetMaximum(1.3);
  scf->GetXaxis()->SetTitle("Template Bins");
  scf->GetYaxis()->SetTitle("u-data/MC");
  scf->SetStats(0);

  // Y axis ratio plot settings
  scf->GetYaxis()->SetNdivisions(505);
  scf->GetYaxis()->SetTitleSize(25);
  scf->GetYaxis()->SetTitleFont(43);
  scf->GetYaxis()->SetTitleOffset(1.35);
  scf->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  scf->GetYaxis()->SetLabelSize(15);
    
  // X axis ratio plot settings
  scf->GetXaxis()->SetTitleSize(25);
  scf->GetXaxis()->SetTitleFont(43);
  scf->GetXaxis()->SetTitleOffset(4.5);
  scf->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel 
  scf->GetXaxis()->SetLabelSize(15);
   
  scf->Draw("esame");

  //Draw line at 1
  TLine *line1 = new TLine(0,1,nx,1);
  line1->SetLineColor(kBlack);
  line1->Draw();
  
 
   // save plot  
   if (save_plots) {
     TString save_name = "Prefit.png";
     c1->SaveAs(save_name);
  
   }
  
   // do the fit

   fill(nfl, nx, nsf, h, hs);

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
   //legend for scalded histograms
   TLegend* l2 = new TLegend(0.65,0.7,0.95,0.9);
   l2->SetBorderSize(0);
   l2->SetFillStyle(0);
   l2->SetTextFont(42);
   l2->SetMargin(0.4);
   
   double* sff= par; // flavor fraction scale factors
   double* sfl = par+nfl; // light scale factors in sensible bins
   
   cout << "0: " << sff[0] << endl;
   cout << "1: " << sff[1] << endl;
   cout << "2: " << sff[2] << endl;
   
   //total fit for ratio
   TH1* totf = new TH1D("","",nx,0,nx);

   TH1* hd[nfl];
   THStack* thss = new THStack();
   for (int jfl = 0; jfl<nfl; ++jfl) {
      TString sh = "jet_MV2c10_";
      sh += chfl[jfl];	
      // sh = h[jfl]->GetName();
      sh += "_scaled";
      hd[jfl] = (TH1*)h[jfl]->Clone(sh);
      for (int ix = 1; ix<=nx; ++ix) {
       double val = h[jfl]->GetBinContent(ix);
       val *= sff[jfl];
       if (ix<=nsf) val *= sfl[ix-1];
       hd[jfl]->SetBinContent(ix,val);
      
     }
      totf->Add(hd[jfl]);
      hd[jfl]->SetFillColor(icol[jfl]);
      thss->Add(hd[jfl]);
      TString sl = chfl[jfl]; sl += " jets";
      l2->AddEntry(hd[jfl],sl,"f");
   }
   l2->AddEntry(hs,"Data","pl");
   
  
 // post-fit plot w/ ratio
  
  TString su = "cf"+to_string(i);
  TCanvas* c2 = new TCanvas(su,su);
  c2->SetCanvasSize(1000,900);
  c2->SetWindowSize(1000,900);

  //upper plot
  TPad *pad3 = new TPad("pad3","pad3", 0,0.3,1,1.0);
  pad3->SetBottomMargin(0); // Upper and lower plot are joined
  pad3->SetGridx();         // Vertical grid
  pad3->Draw();             // Draw the upper pad: pad1
  pad3->cd();               // pad1 becomes the current pad
  pad3->SetLogy();

  thss->Draw();
  thss->GetXaxis()->SetTitle("Template bins");
  thss->GetYaxis()->SetTitle("Arbitrary units");
  // ths->GetYaxis()->SetTitleSize(20); //ERROR: GETS RID OF TITLE
  thss->GetYaxis()->SetTitleOffset(1.25);
  thss->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  thss->GetYaxis()->SetLabelSize(15);
  thss->SetTitle("MV2c10 MC vs l-jet Data {Post-fit}");
  thss->SetMinimum(1e2);
  hs->SetMarkerStyle(8);
  hs->SetMarkerSize(1.2);
  hs->Draw("esame");
  l2->Draw();
  
  
  c2->cd();
  TPad *pad4 = new TPad("pad4","pad4", 0, .01, 1, 0.4);
  pad4->SetTopMargin(0.04);
  pad4->SetBottomMargin(0.4);
  pad4->SetGridx(); // vertical grid
  pad4->Draw();
  pad4->cd(); 
 
  //pre-fit ratio plot on pad 2
  
  TH1D *scff = (TH1D*)hs->Clone("sff");
  scff->SetLineColor(kBlack);
  scff->Sumw2();
  scff->Divide(totf);
  scff->SetMarkerStyle(21);
  scff->SetMarkerStyle(8);
  scff->SetMarkerSize(0.8);
  scff->SetMinimum(0.8);
  scff->SetMaximum(1.3);
  scff->GetXaxis()->SetTitle("Template Bins");
  scff->GetYaxis()->SetTitle("u-data/MC");
  scff->SetStats(0);

  // Y axis ratio plot settings
  scff->GetYaxis()->SetNdivisions(505);
  scff->GetYaxis()->SetTitleSize(25);
  scff->GetYaxis()->SetTitleFont(43);
  scff->GetYaxis()->SetTitleOffset(1.35);
  scff->GetYaxis()->SetLabelFont(43); // Absolute font size in pixel (precision 3)
  scff->GetYaxis()->SetLabelSize(15);
    
  // X axis ratio plot settings
  scff->GetXaxis()->SetTitleSize(25);
  scff->GetXaxis()->SetTitleFont(43);
  scff->GetXaxis()->SetTitleOffset(4.5);
  scff->GetXaxis()->SetLabelFont(43); // Absolute font size in pixel 
  scff->GetXaxis()->SetLabelSize(15);
   
  scff->Draw("esame");

  //Draw line at 1
  TLine *line2 = new TLine(0,1,nx,1);
  line2->SetLineColor(kBlack);
  line2->Draw();
  
 // save plot  
   if (save_plots) {
     TString save_name = "Postfit.png";
     c2->SaveAs(save_name);
   }
}
