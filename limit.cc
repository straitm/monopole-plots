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

  double T = Constants::LIVE_TIME_DDT_CORRECTED;

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

void draw_limits(lim_t lims)
{
  std::map<std::string, TGraph*> g;
  g["half"] = new TGraph();
  g["full"] = new TGraph();

  int n_point = 0;
  for(unsigned int i = 0; i < lims.size(); i++)
  {
    Limit_Info info = lims[i];

    if (info.log_beta() > -3.6 && info.log_beta() < -2.2)
    {
      g.at("half")->SetPoint(n_point, info.log_beta(), info.limit("half"));
      g.at("full")->SetPoint(n_point, info.log_beta(), info.limit("full"));
      ++n_point;
    }
  }

  TCanvas *can = new TCanvas;

  can->SetLeftMargin(0.12);
  can->SetRightMargin(0.05);
  can->SetBottomMargin(0.12);
  can->SetTopMargin(0.09);
  
  g.at("half")->Draw("APC");
  g.at("full")->Draw("SAME PC");
  
  g.at("half")->SetMarkerStyle(kOpenSquare);
  g.at("full")->SetMarkerStyle(kFullSquare);

  g.at("half")->SetLineStyle(kDashed);
  g.at("full")->SetLineStyle(kSolid);
  
  TAxis *x = g.at("half")->GetXaxis();
  x->SetTitle("Monopole Velocity (#beta)");
  x->SetTitleOffset(1.2);
  x->SetTitleSize(0.04);
  x->CenterTitle();
  x->SetNdivisions(6, 5, 0, true);
  x->SetLabelOffset(0.012);

  TAxis *y = g.at("half")->GetYaxis();
  y->SetTitle("90% C.L. Upper Flux Limit (cm^{-2}s^{-1}sr^{-1})");
  y->CenterTitle();
  y->SetTitleOffset(1.3);
  y->SetTitleSize(0.04);
  y->SetRangeUser(1e-16, 1e-12);

  TLegend *l = new TLegend(0.5, 0.73, 0.91, 0.87);
  l->AddEntry
    (g.at("half"), "#Omega = 2#pi, m > 5 #times 10^{8} GeV/c^{2}", "pl");
  l->AddEntry
    (g.at("full"), "#Omega = 4#pi, m > 2 #times 10^{15} GeV/c^{2}", "pl");
  l->Draw();

  can->SaveAs("limit_plot.pdf");
}


lim_t extract_limits(std::string file_name)
{
  std::cout << "extracting limits from " << file_name << std::endl;

  TFile *file = new TFile(file_name.c_str(), "read");
  TTree *tree = dynamic_cast<TTree*>(file->Get("mono/Event"));
  tree_t t;
  std::string *input_file_name = new std::string("invalid");
  tree->SetBranchAddress("input_file_name", &input_file_name);

  for(unsigned int i = 0; i < Constants::nbranch; i++)
    tree->SetBranchAddress(Constants::BRANCH_NAMES[i], &t[Constants::BRANCH_NAMES[i]]);

  lim_t lims;
  
  for (int entry = 0; entry != tree->GetEntries(); ++entry)
  // for (int entry = 0; entry != 10000; ++entry)
  {
    if (entry % 100000 == 0)
      std::cout << "processing entry " << entry << " ..." << std::endl;

    tree->GetEntry(entry);

    if (*input_file_name == "invalid")
    {
      std::cerr << "\n\n\t\tThe input_file_name branch is not set!\n"
		<< std::endl;
      
      assert(false);
    }

    Event_Info e(t, *input_file_name);

    double beta_f = e.beta_from_file_name();
    if (lims.find(beta_f) == lims.end())
      lims[beta_f] = Limit_Info(beta_f);

    lims.at(beta_f).increment_counts(e.is_signal());

    if (e.is_signal())
      lims.at(beta_f).add_area(e.area_projected());
  }

  return lims;
}

void limit()
{
  std::map<std::string, lim_t> l;
  // l["low"] = extract_limits(Constants::MC_RECO_FILE_LOW_ENERGY);
  // l["high"] = extract_limits(Constants::MC_RECO_FILE_HIGH_ENERGY);
  l["0.9dEdx"] = extract_limits(Constants::MC_RECO_FILE);

  // draw_limits(l["low"]);
  // draw_limits(l["high"]);
  draw_limits(l["0.9dEdx"]);
  // for (auto const& lim : l["0.9dEdx"])
  //   lim.second.print();

  // determine the maximum limit
  // for (auto const& lim : l["low"])
  // {
  //   double beta = lim.first;
  //   if (l["low"].at(beta).limit("half") > l["high"].at(beta).limit("half"))
  //     l["max"][beta] = l["low"].at(beta);
  //   else
  //     l["max"][beta] = l["high"].at(beta);
  // }
  // draw_limits(l["max"], "Max");

  // std::cout << "\n\nFinal Limits" << std::endl;
  // for (auto const& lim : l["max"])
  //   lim.second.print();
}
