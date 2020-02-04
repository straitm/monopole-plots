#include "MFRoot.hh"

#include <TLatex.h>
#include <TPaveStats.h>

#include <iostream>
#include <sstream>



double MFRoot::calculate_integral(TH1* hist, double x_min, double x_max)
{
  TAxis *x = hist->GetXaxis();
  int min_bin = x->FindBin(x_min);
  int max_bin = x->FindBin(x_max);

  return hist->Integral(min_bin, max_bin);
}



void MFRoot::copy(hist_t & h, std::string from, std::string to)
{
  h[to] = dynamic_cast<TH1D*>(get(h, from)->Clone(to.c_str()));
}



void MFRoot::fill(tree_t const& t, hist_t const& h, std::string const name)
{
  h.at(name)->Fill(t.at(name));
}



TH1* MFRoot::get(TFile* file, std::string hist_name)
{
  TH1* result = dynamic_cast<TH1*>(file->Get(hist_name.c_str()));

  if (!result)
  {
    std::cerr << "\n" << hist_name << " does not exist in " 
	      << file->GetName() << " \n" << std::endl;
    assert(false);
  }

  result->SetLineWidth(2);

  return result;
}



TH2* MFRoot::get2D(TFile* file, std::string hist_name)
{
  TH2* result = dynamic_cast<TH2*>(file->Get(hist_name.c_str()));

  if (!result)
  {
    std::cerr << "\n" << hist_name << " does not exist in " 
	      << file->GetName() << " \n" << std::endl;
    assert(false);
  }

  return result;
}



TH1* MFRoot::integral_plot(TH1 *hist)
{
  std::string name(hist->GetName());
  name += "_integral_plot";

  TH1* result = dynamic_cast<TH1D*>(hist->Clone(name.c_str()));
  result->SetTitle("Accumulation Plot (High #rightarrow Low)");

  double integral = 0.0;
  int overflow_bin = hist->GetNbinsX() + 1;
  for (int bin = overflow_bin; bin >= 1; --bin)
  {
    integral += hist->GetBinContent(bin);
    result->SetBinContent(bin, integral);
  }

  return result;
}



TFile* MFRoot::open(std::string file_name)
{
  TFile* result = new TFile(file_name.c_str(), "READ");

  if (!result->IsOpen())
  {
    std::cerr << "\n\t" << file_name << " does not exist.\n" 
	      << std::endl;
    assert(false);
  }

  return result;
}



void MFRoot::resize_stats_box(TH1* hist,
			      double const x1, double const x2,
			      double const y1, double const y2)
{
  /*
    The Draw() function must be called before and after this method!
    The first call creates the statistics box and the second call applies
    the new settings.
  */
  
  gPad->Update();
  TPaveStats *stats_box = dynamic_cast<TPaveStats*>(hist->FindObject("stats"));
  if (!stats_box)
  {
    std::cerr << "\n\tStatistics box not found!"
	      << "\n\tCall hist->Draw() before AND after this method.\n"
	      << std::endl;
    assert(false);
  }

  stats_box->SetX1NDC(x1);
  stats_box->SetX2NDC(x2);

  stats_box->SetY1NDC(y1);
  stats_box->SetY2NDC(y2);
}



void MFRoot::wilson_interval
(double const& k, double const& n, double & center, double & error)
{
  // See https://hep.baylor.edu/elog/mfrank/125 or search 
  // "Wilson Score Interval" for a derivation of the below.
  // We pick z = 1 for a 68.3% confidence interval

  double p;

  if (n > 0)
  {
    p = k/n;
    center = (p + (1/(2*n))) / (1 + (1/n));
    error = sqrt( (1/n) * (p*(1-p) + (1/(4*n))) ) / (1+(1/n));
  } else {
    center = 0;
    error = 0;
  }
}



TH1* MFRoot::wilson_interval
(TH1* numerator, TH1* denominator, std::string const& result_hist_name)
{
  TH1* result = dynamic_cast<TH1*>(numerator->Clone(result_hist_name.c_str()));

  for (int bin = 1; bin <= numerator->GetNbinsX(); ++bin)
  {
    double k = numerator->GetBinContent(bin);
    double n = denominator->GetBinContent(bin);
    double center(0.0), error(0.0);

    wilson_interval(k, n, center, error);

    result->SetBinContent(bin, center);
    result->SetBinError(bin, error);
  }

  return result;
}



TH1* MFRoot::divide_with_poisson_errors
(TH1* numerator, TH1* denominator, std::string const& result_hist_name)
{
  TH1* result = dynamic_cast<TH1*>(numerator->Clone(result_hist_name.c_str()));

  for (int bin = 1; bin <= numerator->GetNbinsX(); ++bin)
  {
    double k = numerator->GetBinContent(bin);
    double n = denominator->GetBinContent(bin);

    double center(0.0), error(0.0);
    if (n > 0)
    {
      center = k            / n;
      error  = std::sqrt(k) / n;
    }
    
    result->SetBinContent(bin, center);
    result->SetBinError(bin, error);
  }

  return result;
}



//
// Canvas Manager
//

MFRoot::Canvas_Manager::Canvas_Manager() : n_canvas_(0)
{
}



void MFRoot::Canvas_Manager::new_canvas
(std::string canvas_name, bool logx, bool logy, bool gridx, bool gridy,
 int format)
{
  std::string name;
  if (!canvas_name.empty())
  {
    name = canvas_name;
    if (canvases_.find(name) != canvases_.end())
      name += "_" + std::to_string(n_canvas_);
  } else {
    name = "canvas_" + std::to_string(n_canvas_);
  }

  canvases_[name] = new TCanvas(name.c_str(), name.c_str(), format);
  if (logx)
    canvases_[name]->SetLogx();
  if (logy)
    canvases_[name]->SetLogy();

  if (gridx)
  {
    canvases_[name]->SetGridx();
    canvases_[name]->SetTickx();
  }
  if (gridy)
  {
    canvases_[name]->SetGridy();
    canvases_[name]->SetTicky();
  }

  
  canvases_[name]->cd();

  ++n_canvas_;
}



void MFRoot::Canvas_Manager::add_canvas
(TCanvas *canvas, std::string canvas_name)
{
  std::string name;
  if (!canvas_name.empty())
  {
    name = canvas_name;
    if (canvases_.find(name) != canvases_.end())
      name += "_" + std::to_string(n_canvas_);
  } else {
    name = "canvas_" + std::to_string(n_canvas_);
  }

  canvases_[name] = canvas;
  // canvases_[name]->cd();
  ++n_canvas_;
}



TCanvas* MFRoot::Canvas_Manager::get(std::string canvas_name)
{
  auto result_it = canvases_.find(canvas_name);

  if (result_it == canvases_.end())
  {
    std::cerr << "Cannot get " << canvas_name 
	      << " because it does not exist."
	      << std::endl;
    
    assert(false);
  }

  return result_it->second;
}



void MFRoot::Canvas_Manager::save
(std::string canvas_name, std::string path, std::string extension)
{
  if (canvases_.find(canvas_name) == canvases_.end())
  {
    std::cerr << "Cannot save " << canvas_name 
	      << " because it does not exist."
	      << std::endl;
  } else {
    std::ostringstream stream;
    stream << path << canvas_name << "." << extension;
    canvases_[canvas_name]->SaveAs(stream.str().c_str());
  }
}



void MFRoot::Canvas_Manager::save_all(std::string path)
{
  std::vector<std::string> extensions = { "pdf", "gif" };

  for (auto const& ext : extensions)
    save_all(path, ext);

  // save_all(path, "pdf");
  // save_all(path, "gif");
  // save_all(path, "eps");
}



void MFRoot::Canvas_Manager::save_all(std::string path, std::string extension)
{
  for (auto const& can : canvases_)
    save(can.first, path, extension);
}



MFRoot::Tree_Info::Tree_Info()
{
}



MFRoot::Tree_Info::Tree_Info(tree_t const& t) : t_(t)
{
}



double MFRoot::Tree_Info::get(std::string branch_name) const
{
  return MFRoot::get(t_, branch_name);
}



void nova::Preliminary()
{
  // Put a "NOvA Preliminary" tag in the corner
  TLatex* prelim = new TLatex(.9, .95, "NOvA Preliminary");
  prelim->SetTextColor(kBlue);
  prelim->SetNDC();
  prelim->SetTextSize(2/30.);
  prelim->SetTextAlign(32);
  prelim->Draw();
}



void nova::Simulation()
{
  // Put a "NOvA Simulation" tag in the corner
  TLatex* prelim = new TLatex(.9, .95, "NOvA Simulation");
  prelim->SetTextColor(kGray+1);
  prelim->SetNDC();
  prelim->SetTextSize(2/30.);
  prelim->SetTextAlign(32);
  prelim->Draw();
}
