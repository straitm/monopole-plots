/*
  This program calculates the 90% C.L. upper limit on the monopole flux
  for a given velocity point.
 */

#include "Constants.hh"
#include "Event_Info.hh"
#include "MFRoot.hh"

#include <TFile.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TTree.h>
#include <TStyle.h>
#include <TArrow.h>

#include <iostream>
#include <fstream>
#include <map>
#include <string>
#include <vector>

#include <cmath>
#include <iomanip>

#define DRAWHEAVY

using namespace MFRoot;

static double textsize = 0.057;

/*
  This class keeps track of the cumulative information needed to calculate
  the final 90% C.L. upper limits.
 */

class Limit_Info
{
public:
  Limit_Info();
  Limit_Info(double beta);

  void   add_area(double area);
  double area() const { return area_; }
  double beta() const;
  double efficiency() const;
  void   increment_counts(bool is_signal);
  double limit(std::string coverage) const;
  double log_beta() const;
  double n_monopoles_total() const;
  double n_monopoles_signal() const;
  void   print() const;
  double projected_area_total() const;
  
private:
  double area_;
  double beta_;
  double n_;
  double n_signal_;
};

Limit_Info::Limit_Info()
{
}

Limit_Info::Limit_Info(double beta) :
  area_(0.0),
  beta_(beta),
  n_(0.0),
  n_signal_(0.0)
{
}


void Limit_Info::add_area(double area)
{
  area_ += area;
}


double Limit_Info::beta() const
{
  return beta_;
}


double Limit_Info::efficiency() const
{
  return n_signal_ / n_;
}


void Limit_Info::increment_counts(bool is_signal)
{
  ++n_;

  if (is_signal)
    ++n_signal_;
}


double Limit_Info::limit(std::string coverage) const
{
  double result = 0.0;

  double omega = M_PI;
  if (coverage == "half")
    omega *= 2;
  else if (coverage == "full")
    omega *= 4;
  else
    std::cerr << "Coverage value (" << coverage << ") is not valid."
	      << std::endl;

  double T = LIVE_TIME_DDT_CORRECTED;

  if (projected_area_total() > 0)
  {
    // convert area to cm^2
    double A_cm2 = projected_area_total() * 1e4;

    double eA    = A_cm2 / n_monopoles_total();
    result       = 2.3 / (omega * T * eA);
  }

  return result;
}


double Limit_Info::log_beta() const
{
  return std::log10(beta_);
}


double Limit_Info::n_monopoles_total() const
{
  return n_;
}


double Limit_Info::n_monopoles_signal() const
{
  return n_signal_;
}


void Limit_Info::print() const
{
  std::cout << std::setw(12) << log_beta() << ", "
	    << std::setw(12) << efficiency() << ", "
	    << std::setw(12) << limit("half") << ", "
	    << std::setw(12) << limit("full")
	    << std::endl;
}


double Limit_Info::projected_area_total() const
{
  return area_;
}

typedef std::map<double, Limit_Info> lim_t;

void draw_limits(const lim_t & lims)
{
  std::map<std::string, TGraph*> g;
  g["half"] = new TGraph();
  g["full"] = new TGraph();

  std::ifstream infile("limitdata.txt");
  if(infile.is_open()){
    double logbeta, half, full;
    while(infile >> logbeta >> half >> full){
      g["half"]->SetPoint(g["half"]->GetN(), pow(10, logbeta), half);
      g["full"]->SetPoint(g["full"]->GetN(), pow(10, logbeta), full);
    }
  }
  else{
    int n_point = 0;
    std::ofstream outfile("limitdata.txt");
    if(!outfile.is_open()){
      fprintf(stderr, "Could not open limitdata.txt for writing\n");
      exit(1);
    }
    for(const auto & lim : lims)
    {
      const Limit_Info & info = lim.second;

      if (info.log_beta() > -3.75 && info.log_beta() < -1.05)
      {
        if(info.limit("full")*1e15 > 10)
          printf("  $%3.1f$ & %02.0lf & %.0lf\\phantom{.0} \\\\\n",
            info.log_beta(), info.limit("half")*1e15, info.limit("full")*1e15);
        else
          printf("  $%3.1f$ & %02.0lf & \\phantom{0}%.1lf \\\\\n",
            info.log_beta(), info.limit("half")*1e15, info.limit("full")*1e15);

        g.at("half")->SetPoint(n_point, info.beta(), info.limit("half"));
        g.at("full")->SetPoint(n_point, info.beta(), info.limit("full"));
        ++n_point;

        outfile << info.log_beta() << " " << info.limit("half") << " "
                << info.limit("full") << std::endl;
      }
    }
  }

  TCanvas *can = new TCanvas;

  can->SetCanvasSize(600, 400);
  can->SetRightMargin(0.015);
  can->SetTopMargin(0.03);
  can->SetLeftMargin(0.13);
  can->SetBottomMargin(0.13);

  can->SetFrameLineWidth(2);

  can->SetLogy();

  can->SetLogx();

  const double xmin = pow(10, -4.4);

  TH2D dum("dum", "", 1, xmin, 1, 1, 1e-18, 1e-12);
  
  dum.Draw();
  
  g.at("half")->SetLineWidth(2);
  g.at("half")->SetFillStyle(1001);
  g.at("half")->SetFillColor(kYellow);
  g.at("half")->SetLineColor(kBlack);
  g.at("full")->SetLineWidth(2);

  TAxis *x = dum.GetXaxis();
  x->SetTitle("Monopole speed (#beta)");
  x->SetTitleSize(textsize);
  x->SetLabelSize(textsize);
  x->CenterTitle();
  x->SetNdivisions(6, 5, 0, true);

  TAxis *y = dum.GetYaxis();
  y->SetTitle("90% C.L. flux limit (cm^{#minus2}#kern[-0.5]{ }"
              "s^{#minus1}#kern[-0.5]{ }sr^{#minus1}#kern[-0.9]{ })");
  y->CenterTitle();
  y->SetTitleOffset(1.15);
  y->SetLabelOffset(0.001);
  x->SetTitleOffset(1.12);
  y->SetTitleSize(textsize);
  y->SetLabelSize(textsize);
  y->SetRangeUser(1e-16, 1e-12);

  x->SetTickSize(0.025); // Smaller than default (0.03)
  y->SetTickSize(0.025); 

  const double thisworky = 0.76;
  const double thisworkx = 0.335;
  TLegend *l = new TLegend(thisworkx,      thisworky,
                           thisworkx+0.15, thisworky +0.18
  );
  l->SetTextSize(textsize);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->SetTextFont(42);
  l->SetTextAlign(22);
  l->AddEntry((TH1D*)NULL, "NOvA", "");
  l->AddEntry((TH1D*)NULL, "(This work)", "");
  l->AddEntry
    (g.at("half"), "> 5#kern[-0.5]{ }#times#kern[-0.9]{ }10^{8}#kern[-0.3]{ }GeV",
  "");

  TGraph slimlight;
  slimlight.SetPoint(slimlight.GetN(), 1.00, 1e-11);
  slimlight.SetPoint(slimlight.GetN(), 0.05, 1e-11);
  slimlight.SetPoint(slimlight.GetN(), 0.05, 1.3e-15);
  slimlight.SetPoint(slimlight.GetN(), 1.00, 1.3e-15);

  slimlight.SetLineWidth(2);

  slimlight.SetFillStyle(1001);
  slimlight.SetFillColor(kGreen);
 
  TGraph slimheavy; // irrelevant
  slimheavy.SetPoint(slimheavy.GetN(), 0.05, 0.65e-15);
  slimheavy.SetPoint(slimheavy.GetN(), 0.00, 0.65e-15);

  slimheavy.SetLineWidth(2);

  // Extracted from the plot in the MACRO paper
  TGraph macroheavy;
  macroheavy.SetPoint(macroheavy.GetN(), (4.e-05), 3.111271086997062e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (4.63003e-05), 2.7495747867546176e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (5.05009e-05), 2.559717811217955e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (6.11617e-05), 2.3293333163141835e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (8.0847e-05), 2.1544346900318954e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (9.05598e-05), 2.1128017944843972e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (0.0001), 2.0991037201085632e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (0.0001), 1.5843778311510706e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (0.00011), 1.5869563145973281e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (0.00011), 1.272095668992179e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (0.00012), 1.2700287704932164e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (0.00014), 1.29084978338489e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (0.00018), 1.3597985579403563e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (0.00018), 1.3356916613034106e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (0.00020), 1.3686721650920103e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (0.000256597), 1.4701880261868488e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (0.000301783), 1.5792334259956746e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (0.000541903), 1.5588223004671087e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (0.000607082), 1.553760872996211e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (0.000735782), 1.5588223004671087e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (0.0009235), 1.566445373014912e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (0.0012), 1.5715481206644581e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (0.0012), 1.4324301262444822e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (0.00165697), 1.3911103503173643e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (0.00298785), 1.3400427210994613e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (0.00361065), 1.6234988203559148e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (0.00404453), 1.5973122800602588e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (0.005), 1.5689946724208885e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (0.005), 1.6608827826277233e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (0.006), 1.63143817898912e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (0.00689217), 1.59212587731796e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (0.0100762), 1.5237355654616932e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (0.0203527), 1.4441241135747835e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (0.0301526), 1.4116190614769622e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (0.044674), 1.3911103503173643e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (0.0622653), 1.3820912681235753e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (0.0998089), 1.3776036785441396e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (0.100726), 1.4535480021259694e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (0.22115), 1.4488283956730566e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (0.659272), 1.44647434219323e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (1), 1.44647434219323e-16);
  macroheavy.SetPoint(macroheavy.GetN(), (1), 1e-12);
  macroheavy.SetPoint(macroheavy.GetN(), (xmin), 1e-12);

  TGraph macrolight;
  for(int i = 0; i < macroheavy.GetN(); i++)
    macrolight.SetPoint(i, macroheavy.GetX()[i], macroheavy.GetY()[i]*2);

  macrolight.SetPoint(macrolight.GetN(),
    macroheavy.GetX()[macroheavy.GetN()-1], 1e-11);

  macrolight.SetPoint(macrolight.GetN(),
    macroheavy.GetX()[0], 1e-11);

  macroheavy.SetLineWidth(2);
  macroheavy.SetFillStyle(1001);
  macroheavy.SetLineColor(kBlack);
  macroheavy.SetFillColor(kCyan-10);


  macrolight.SetLineWidth(2);
  macrolight.SetFillStyle(1001);
  macrolight.SetLineColor(kBlack);
  macrolight.SetFillColor(kCyan);

  TGraph icecube;
  icecube.SetPoint(icecube.GetN(), (0.995), 1e-11); // just for drawing
  icecube.SetPoint(icecube.GetN(), (0.8), 1e-11); // just for drawing
  icecube.SetPoint(icecube.GetN(), (0.8), 5.57e-18);
  icecube.SetPoint(icecube.GetN(), (0.9), 3.56e-18);
  icecube.SetPoint(icecube.GetN(), (0.995), 3.38e-18);

  icecube.SetFillStyle(1001);
  icecube.SetFillColor(kRed);
  icecube.SetLineColor(kBlack);
  icecube.SetLineWidth(2);


#ifdef DRAWHEAVY
  macroheavy.Draw("lf"); // 1e16
#endif
  macrolight   .Draw("lf"); // 1e10
  g.at("half")->Draw("lf"); // 5e8
  icecube      .Draw("lf"); // 1e6
  slimlight    .Draw("lf"); // 1e5


  l->Draw();

  {
    const double icex = 0.725,
                 icey = 0.173;
    TLegend *l = new TLegend(icex,      icey,
                             icex+0.15, icey+0.12);
    l->SetTextSize(textsize);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextFont(42);
    l->AddEntry((TH1D*)NULL, "IceCube", "");
    l->AddEntry((TH1D*)NULL, "> 10^{6}#kern[-0.3]{ }GeV", "");
    l->Draw();

    TArrow * a = new TArrow(0.36, 1.0e-17,
                            0.75, 1.0e-17, 0.011, "|>");
    a->SetLineWidth(2);
    a->Draw();
  }
  {
    const double macrox = 0.36,
                 macroy = 0.54;
    TLegend *l = new TLegend(macrox,      macroy,
                             macrox+0.15, macroy+0.12);
    l->SetTextSize(textsize);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextFont(42);
    l->SetTextAlign(22);
    l->AddEntry((TH1D*)NULL, "MACRO", "");
    l->AddEntry(&macrolight, "> 10^{10}#kern[-0.3]{ }GeV", "");
    l->Draw();
  }

  #ifdef DRAWHEAVY
  {
    const double macrox = 0.36,
                 macroy = 0.25;
    TLegend *l = new TLegend(macrox,      macroy,
                             macrox+0.15, macroy+0.12
                             );
    l->SetTextSize(textsize);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextFont(42);
    l->SetTextAlign(22);
    l->AddEntry((TH1D*)NULL, "MACRO", "");
    l->AddEntry((TH1D*)NULL, "> 10^{16}#kern[-0.3]{ }GeV", "");
    l->Draw();

    TArrow * a = new TArrow(1.85e-3, 6.0e-17,
                            1.85e-3, 2.0e-16, 0.011, "|>");
    a->SetLineWidth(2);
    a->Draw();
  }
  #endif

  {
    const double slimy = 0.69,
                 slimx = 0.76;
    TLegend *l = new TLegend(slimx,      slimy,
                             slimx+0.15, slimy+0.12);
    l->SetTextSize(textsize);
    l->SetBorderSize(0);
    l->SetTextAlign(22);
    l->SetFillStyle(0);
    l->SetTextFont(42);
    l->AddEntry((TH1D*)NULL, "SLIM", "");
    l->AddEntry(&slimlight, "> 10^{5}#kern[-0.3]{ }GeV",
#if 0
    "l");
    l->AddEntry(&slimheavy,
      "> 5 #times 10^{13} GeV", "l");
#else
    "");
#endif
    l->Draw();
  }

  can->RedrawAxis();
  can->SaveAs("limit_plot.pdf");
}


lim_t extract_limits(const std::string & file_name)
{
  std::cout << "extracting limits from " << file_name << std::endl;

  TFile *file = new TFile(file_name.c_str(), "read");
  TTree *tree = dynamic_cast<TTree*>(file->Get("mono/Event"));
  tree_t t;
  std::string *input_file_name = new std::string("invalid");
  tree->SetBranchAddress("input_file_name", &input_file_name);

  for(unsigned int i = 0; i < nbranch; i++)
    tree->SetBranchAddress(BRANCH_NAMES[i], &t[BRANCH_NAMES[i]]);

  lim_t lims;

  for (int entry = 0; entry != tree->GetEntries(); ++entry)
  {
    if (entry % 100000 == 0)
      std::cout << "processing entry " << entry << " of " << tree->GetEntries()
                << " ..." << std::endl;

    tree->GetEntry(entry);

    #if 0
    if (*input_file_name == "invalid"){
      std::cerr << "\n\n\t\tThe input_file_name branch is not set!\n" << std::endl;
      _exit(1);
    }
    #endif

    Event_Info e(t, *input_file_name);
    const double beta_f = e.beta_from_file_name();

    if (lims.find(beta_f) == lims.end())
      lims[beta_f] = Limit_Info(beta_f);

    lims.at(beta_f).increment_counts(e.is_signal());

    if (e.is_signal())
      lims.at(beta_f).add_area(e.area_projected());
  }

  return lims;
}

int main()
{
  gStyle->SetOptStat(0);
 
  std::map<std::string, lim_t> l;
  {
    std::ifstream infile("limitdata.txt");
    if(!infile.is_open())
      l["0.9dEdx"] = extract_limits(MC_RECO_FILE);
  }

  draw_limits(l["0.9dEdx"]);
}
