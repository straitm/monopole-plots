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

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <cmath>
#include <iomanip>

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

  int n_point = 0;
  for(const auto & lim : lims)
  {
    const Limit_Info & info = lim.second;

    printf("log beta = %f\n", info.log_beta());

    if (info.log_beta() > -3.6 && info.log_beta() < -2.2)
    {
      g.at("half")->SetPoint(n_point, info.log_beta(), info.limit("half"));
      g.at("full")->SetPoint(n_point, info.log_beta(), info.limit("full"));
      ++n_point;
    }
  }

  TCanvas *can = new TCanvas;

  can->SetCanvasSize(600, 400);
  can->SetRightMargin(0.025);
  can->SetTopMargin(0.03);
  can->SetLeftMargin(0.14);
  can->SetBottomMargin(0.14);

  can->SetLogy();
  can->SetTickx();
  can->SetTicky();
  
  g.at("half")->Draw("APC");
  g.at("full")->Draw("SAME PC");
  
  g.at("half")->SetMarkerStyle(kOpenSquare);
  g.at("full")->SetMarkerStyle(kFullSquare);

  g.at("half")->SetLineStyle(kDashed);
  g.at("full")->SetLineStyle(kSolid);
  
  TAxis *x = g.at("half")->GetXaxis();
  x->SetTitle("Monopole Velocity (#beta)");
  x->SetTitleOffset(1.1);
  x->SetTitleSize(textsize);
  x->SetLabelSize(textsize);
  x->CenterTitle();
  x->SetNdivisions(6, 5, 0, true);

  TAxis *y = g.at("half")->GetYaxis();
  y->SetTitle("90% C.L. Flux Limit (cm^{#minus2}s^{#minus1}sr^{#minus1})");
  y->CenterTitle();
  y->SetTitleOffset(1.2);
  y->SetTitleSize(textsize);
  y->SetLabelSize(textsize);
  y->SetRangeUser(1e-16, 1e-12);

  TLegend *l = new TLegend(0.3, 0.77, 0.71, 0.91);
  l->SetTextSize(textsize);
  l->SetBorderSize(0);
  l->SetFillStyle(0);
  l->SetTextFont(42);
  l->AddEntry
    (g.at("half"), "m > 5 #times 10^{8} GeV/c^{2} (#Omega = 2#pi)", "pl");
  l->AddEntry
    (g.at("full"), "m > 2 #times 10^{15} GeV/c^{2} (#Omega = 4#pi)", "pl");
  l->Draw();

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
  std::map<std::string, lim_t> l;
  l["0.9dEdx"] = extract_limits(MC_RECO_FILE);

  draw_limits(l["0.9dEdx"]);
}
