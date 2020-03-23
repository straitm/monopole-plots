/*
  This macro (originally written by Enhao Song), takes the coverage n-tuples and
  plots them on a 2D histogram.
 */

#include "MFRoot.hh"

#include <TColor.h>
#include <TCanvas.h>
#include <TFile.h>
#include <TH2.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMath.h>
#include <TNtuple.h>
#include <TStyle.h>
#include <TTree.h>

#include <cmath>
#include <cstdlib>
#include <iostream>

const double textsize = 0.057;

TCanvas *c1 = new TCanvas("c1", "Coverage");

void fill_histogram(TNtuple *ntuple, TH2F *hist, bool verbose)
{
  // extract branches from the n-tuple:
  float beta(0.0), mass(0.0), angle(0.0);
  ntuple->SetBranchAddress("beta", &beta);
  ntuple->SetBranchAddress("mass", &mass);
  ntuple->SetBranchAddress("angle", &angle);

  // fill histogram:
  float last_beta(0), last_mass(0);
  for (int i = 0; i < ntuple->GetEntries(); i++)
  {
    ntuple->GetEntry(i);
    
    int bin = hist->FindBin(TMath::Log10(mass), TMath::Log10(beta));
    // the only angles that we consider here are:
    //   * angle = 0 (traversing entire planet)
    //   * angle = 90 (entering detector horizontally)
    // this yields only two possible values for content: 0.5 or 1.0
    double content =
      (4 * TMath::Pi() - 2 * TMath::Pi() * (1 - std::cos(angle/180.0 * TMath::Pi()))) / (4 * TMath::Pi());
    hist->SetBinContent(bin, content);

    if (verbose)
      std::cout << "data\t"
		<< beta << "\t"
		<< mass << "\t"
		<< angle << std::endl;

    // fill out remaining column:
    if (last_mass == mass)
    {
      for (double current_beta = last_beta;
    	   current_beta / 1.1220184543 > beta;
    	   current_beta /= 1.1220184543)
      {
    	int bin_2 =
    	  hist->FindBin(TMath::Log10(mass), TMath::Log10(current_beta));
        hist->SetBinContent(bin_2, content);
      }
    } else {
      for (double current_beta = 1;
    	   current_beta / 1.1220184543 > beta;
    	   current_beta /= 1.1220184543)
      {
    	int bin_2 =
    	  hist->FindBin(TMath::Log10(mass), TMath::Log10(current_beta));
        hist->SetBinContent(bin_2, content);
      }
    }
    
    last_beta = beta;
    last_mass = mass;
  }
}

void draw(TH2F* h)
{
  gStyle->SetOptStat(0);

  TAxis *x = h->GetXaxis();
  TAxis *y = h->GetYaxis();
  TAxis *z = h->GetZaxis();

  x->SetTitle("log_{10}(Monopole mass (GeV))");
  x->SetTitleOffset(1.1);
  x->SetTitleSize(textsize);
  x->SetLabelSize(textsize);
  x->CenterTitle();
  x->SetNdivisions(15, 0, 0, true);
  x->SetRangeUser(4.5, 18);
  
  y->SetTitle("log_{10}#beta");
  y->SetTitleOffset(0.9);
  y->SetTitleSize(textsize);
  y->SetLabelSize(textsize);
  y->CenterTitle();
  y->SetNdivisions(5, 0, 0, true);
  y->SetRangeUser(-4, 0);

  z->SetRangeUser(0, 1.6);

  c1->SetCanvasSize(600, 400);
  
  c1->SetTickx();
  c1->SetTicky();

  c1->SetRightMargin(0.09);
  c1->SetTopMargin(0.03);
  c1->SetLeftMargin(0.11);
  c1->SetBottomMargin(0.14);

  
  {
    const int NRGBs = 2, NCont = 512;
    gStyle->SetNumberContours(NCont);
    Double_t stops[NRGBs] = { 0.00, 1.00 };
    Double_t red[NRGBs]   = { 1.00, 0.20 };
    Double_t green[NRGBs] = { 1.00, 0.20 };
    Double_t blue[NRGBs]  = { 1.00, 0.20 };
    TColor::CreateGradientColorTable(NRGBs, stops, red, green, blue, NCont);
  }

  h->Draw("COL");

  // add omega labels
  TLatex region_labels;
  region_labels.SetTextFont(42);
  region_labels.SetTextSize(textsize);
  region_labels.SetTextAlign(22);
  region_labels.DrawLatex(8.8 , -0.5, "#Omega = 2#pi");
  region_labels.DrawLatex(15.0, -0.5, "#Omega = 4#pi");

  // add slow monopole beta range lines
  TLine line;
  line.SetLineColor(kBlack);
  line.SetLineStyle(kDashed);
  line.SetLineWidth(3);
  line.DrawLine(4.4, -2.1, 18, -2.1);
  line.DrawLine(4.4, -3.5, 18, -3.5);
  line.DrawLine(log10(5e8) , -3.5, log10(5e8) , -2.1);
  line.DrawLine(log10(2e15), -3.5, log10(2e15), -2.1);

  const double arrowy = (-2.1 -3.5)/2, arrowdx = 1.0;
  TArrow * ahalf = new TArrow(log10(5e8), arrowy, log10(5e8)+arrowdx, arrowy, 0.02, "|>");
  ahalf->SetLineWidth(2);
  ahalf->Draw();

  TArrow * afull = new TArrow(log10(2e15), arrowy, log10(2e15)+arrowdx, arrowy, 0.02, "|>");
  afull->SetLineWidth(2);
  afull->Draw();

  TLatex line_labels;
  line_labels.SetTextSize(textsize);
  line_labels.SetTextFont(42);
  line_labels.SetTextAlign(22);
  line_labels.DrawLatex(18.8, -2.1 + 0.02, "#minus2.1");
  line_labels.DrawLatex(18.8, -3.5 + 0.02, "#minus3.5");
}


void coverage()
{
  std::string ntuple_name ="data/"
    "beta_mass_angle_initial_beta_1e-4_10_meter_overburden_standard_dedx.root";
  
  TFile * file = new TFile(ntuple_name.c_str(), "read");
  
  TNtuple *ntuple = dynamic_cast<TNtuple*>(file->Get("ntuple"));

  TH2F *hist= new TH2F("accessible_region", "",
		       139, 4, 18,
		       40, -4,  0);

  fill_histogram(ntuple, hist, false);

  draw(hist);

  c1->SaveAs("coverage.pdf");
}
