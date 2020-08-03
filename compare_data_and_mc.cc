/*
  Run over Data and MC monopole trees and generate reconstruction histograms 
  comparing Data with MC.
*/

#include "Constants.hh"
#include "Event_Info.hh"
#include "Event_List.hh"

#include <TROOT.h>
#include <TArrow.h>
#include <TFile.h>
#include <TH1.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TStyle.h>
#include <TTree.h>

#include <fstream>
#include <iostream>
#include <string>

static double textsize = 0.057;

typedef std::map<std::string, TH1*> hist_t;
typedef std::map<std::string, hist_t> hists_t;
typedef std::map<std::string, double> tree_t;
typedef std::map<std::string, tree_t> trees_t;

template<typename T>
T get(std::map<std::string, T> const& m, const std::string element)
{
  if (m.find(element) == m.end())
  {
    std::cerr << "\n\n\t\tThe map does not contain the following element: "
              << element << "\n\n"
              << std::endl;
    _exit(1);
  }

  return m.at(element);
}

struct Hist_Info
{
  Hist_Info
  (std::string name_, int n_bins_, double low_bin_, double high_bin_ ) :
    name(name_), n_bins(n_bins_), low_bin(low_bin_), high_bin(high_bin_)
  {}

  std::string name;
  int n_bins;
  double low_bin, high_bin;
};

void declare_hist         (hists_t & h, const std::string type,
			   const std::string name,
			   const int n_bins, const double low_bin,
			   const double high_bin);
void declare_hists        (hists_t & h, const std::string type);
void draw_beta            (hists_t const& h, std::string hist_name);
void draw_dplane          (hists_t const& h, std::string hist_name);
void draw_n_hits          (hists_t const& h, std::string hist_name);
void draw_n_tracks        (hists_t const& h, std::string hist_name);
void draw_projected_length(hists_t const& h, std::string hist_name);
void draw_r2              (hists_t const& h, std::string hist_name);
void draw_theta           (hists_t const& h, std::string hist_name);
void draw_time_gap        (hists_t const& h, std::string hist_name);
void draw_track_length    (hists_t const& h);
void fill_hists           (hist_t & h, tree_t const& t, TTree* tree,
			   std::string type);


void compare_data_and_mc()
{
  gROOT->LoadMacro("Event_List.cc+");
  gROOT->LoadMacro("Event_Info.cc+");

  gStyle->SetOptStat("");

  std::map<std::string, TFile*> files;
  files[DATA_SAMPLE_NAME] = new TFile(DATA_RECO_FILE.c_str(), "read");
  files["MC"] = new TFile(MC_RECO_FILE.c_str(), "read");

  std::map<std::string, TTree*> trees;
  for (auto const& file : files)
    trees[file.first] = dynamic_cast<TTree*>(file.second->Get("mono/Event"));

  trees_t t;
  for (auto const& tree : trees)
    for (auto const& name : BRANCH_NAMES)
      tree.second->SetBranchAddress(name, &(t[tree.first])[name]);


  hists_t h;
  for (auto const& tree : trees) declare_hists(h, tree.first);

  for (auto const& tree : trees) {
    std::string type = tree.first;
    std::cout << "Looping over " << type << " tree ("
	      << tree.second->GetEntries() << " entries)." << std::endl;

    fill_hists(h.at(type), t.at(type), tree.second, type);
  }

  // Let's scale the MC histograms to the same size as the Data histograms
  for (auto & hist : h.at(DATA_SAMPLE_NAME)) {
    std::string hist_name = hist.first;
    TH1* data_hist = hist.second;
    double data_integral  = data_hist->Integral();

    TH1* mc_hist = h.at("MC").at(hist_name);
    double mc_integral = mc_hist->Integral();

    double scale_factor = data_integral / mc_integral;
    if(data_integral > 0) mc_hist->Scale(scale_factor);
  }
  
  for (auto const& hists : h) {
    for (auto const& hist : hists.second) {
      hist.second->SetLineWidth(2);
      hist.second->SetLineColor(kBlack);
    }
  }

  draw_r2(h, "r2 min low gap and slow");
  draw_time_gap(h, "time gap max");
}


void declare_hist
(hists_t & h,
 const std::string type, const std::string name,
 const int n_bins, const double low_bin, const double high_bin)
{
  std::string TH1_name = type + " " + name;

  (h[type])[name] = new TH1D(TH1_name.c_str(), "", n_bins, low_bin, high_bin);
}


void declare_hists(hists_t & h, const std::string type)
{
  std::vector<Hist_Info> infos =
    {
      {"dplane xz", 30, 0, 300},
      {"dplane yz", 30, 0, 300},
      {"dplane min",150, 0, 300},

      {"nhits"   , 50, 0, 500},
      {"nhits xz", 50, 0, 500},
      {"nhits yz", 50, 0, 500},

      {"Number of Reco Tracks",         5, -0.5, 4.5},
      {"First Preselected Reco Track",  5, -0.5, 4.5},

      {"r2 xt" ,                  40, 0, 1},
      {"r2 yt" ,                  40, 0, 1},
      {"r2 min",                  40, 0, 1},
      {"r2 min low gap",          40, 0, 1},
      {"r2 min low gap and slow", 40, 0, 1},
      
      {"Track #Deltax", 32, 0, 16},
      {"Track #Deltay", 32, 0, 16},
      {"Track #Deltaz", 35, 0, 70},
      {"Track Length" , 35, 0, 70},

      {"#theta_{xz}", 37, -92.5, 92.5},
      {"#theta_{yz}", 37, -92.5, 92.5},

      {"time gap xz" , 50, 0, 1},
      {"time gap yz" , 50, 0, 1},
      {"time gap max", 50, 0, 1}
    };

  for (auto const& i : infos)
    declare_hist(h, type, i.name, i.n_bins, i.low_bin, i.high_bin);

  // histograms with special binning:
  double xbins[] = { 1e-6, 3.16e-6, 
		     1e-5, 3.16e-5, 
		     1e-4, 3.16e-4, 
		     1e-3, 3.16e-3, 
		     1e-2, 3.16e-2, 
		     1e-1, 3.16e-1, 1 };
  
  (h[type])["Beta"] =
    new TH1D((type + " Beta").c_str(), "", 12, xbins);
  (h[type])["Beta Linear Events"] =
    new TH1D((type + " Beta Linear Events").c_str(), "", 12, xbins);
}


void fill_hists(hist_t & h, tree_t const& t, TTree* tree, std::string type)
{
  std::ifstream inr2("r2histdat.txt");
  std::ifstream ingap("gaphistdat.txt");
  if(inr2.is_open() && ingap.is_open()){
    printf("Reading in data from cached histogram files r2histdat.txt\n"
           "and gaphistdat.txt and ignoring ROOT files.  Remove the cache\n"
           "files to regenerate\n");

    std::string what;
    double val;
    int r2bin = 1;
    while(inr2 >> what >> val){
      if(type == what)
        get(h, "r2 min low gap and slow")->SetBinContent(r2bin++, val);
    }
    int gapbin = 1;
    while(ingap >> what >> val){
      if(type == what)
        get(h, "time gap max")->SetBinContent(gapbin++, val);
    }

    return;
  }

  std::string *input_file_name = new std::string("invalid");
  if (type != DATA_SAMPLE_NAME)
    tree->SetBranchAddress("input_file_name", &input_file_name);

  Event_List elist(type);
  for (int entry = 0; entry < tree->GetEntries(); ++entry){
  // for (int entry = 0; entry < 1000000; ++entry) {
    tree->GetEntry(entry);
    if(entry%0x10000 == 0) printf("Entry %d\n", entry);
    Event_Info e(t, *input_file_name);

    if (!e.beta_value_matches(MC_BETA_NAME)) continue;

    if (elist.is_duplicate(e.run(), e.subrun(), e.event())) continue;

    get(h, "Number of Reco Tracks")->Fill(e.n_tracks());

    if (!e.is_preselected_reco()) continue;

    const int n_track = e.first_preselected_reco_track();

    get(h, "First Preselected Reco Track")->Fill(n_track);
    get(h, "Beta")       ->Fill(e.get("beta",      n_track));
    get(h, "dplane xz")  ->Fill(e.get("dplane_x",  n_track));
    get(h, "dplane yz")  ->Fill(e.get("dplane_y",  n_track));
    get(h, "nhits")      ->Fill(e.get("n_hits",    n_track));
    get(h, "nhits xz")   ->Fill(e.get("n_hits_x",  n_track));
    get(h, "nhits yz")   ->Fill(e.get("n_hits_y",  n_track));
    get(h, "#theta_{xz}")->Fill(e.reco_theta("xz", n_track));
    get(h, "#theta_{yz}")->Fill(e.reco_theta("yz", n_track));

    const double r2_xt = e.get("r2_xt", n_track);
    const double r2_yt = e.get("r2_yt", n_track);
    const double r2min = std::min(r2_xt, r2_yt);
    get(h, "r2 xt")->Fill(r2_xt);
    get(h, "r2 yt")->Fill(r2_yt);
    get(h, "r2 min")->Fill(r2min);

    const double dt = e.get("dt", n_track);
    const double gap_xz = dt == 0?1:e.get("max_time_gap_xz", n_track) / dt;
    const double gap_yz = dt == 0?1:e.get("max_time_gap_yz", n_track) / dt;
    const double gap_max = std::max(gap_xz, gap_yz);

    if(e.high_r2_min() && e.is_slow()){
      get(h, "time gap xz")->Fill(gap_xz);
      get(h, "time gap yz")->Fill(gap_yz);
      get(h, "time gap max")->Fill(gap_max);
    }

    if(dt == 0)
      std::cerr << "This track has zero time duration!!!" << std::endl;

    get(h, "Track #Deltax")->Fill(e.get("dx"    , n_track) / 100);
    get(h, "Track #Deltay")->Fill(e.get("dy"    , n_track) / 100);
    get(h, "Track #Deltaz")->Fill(e.get("dz"    , n_track) / 100);
    get(h, "Track Length") ->Fill(e.get("length", n_track) / 100);

    if (e.low_gap_max()) {
      get(h, "r2 min low gap")->Fill(std::min(r2_xt, r2_yt));

      if (e.is_slow())
        get(h, "r2 min low gap and slow")->Fill(std::min(r2_xt, r2_yt));
    }

    if (e.is_linear()) {
      int n_track = e.first_preselected_reco_track();
      double beta = e.get("beta", n_track);

      get(h, "Beta Linear Events")->Fill(beta);

      if (type == DATA_SAMPLE_NAME && e.is_slow())
        printf("Selected data: run %5.0f ev %7.0f beta %6.4f "
               "fmax %5.3f r2min %5.3f\n",
          get(t, "run_number"), get(t, "event_number"),
          beta, gap_max, r2min);
    }
  }
}


void draw_beta(hists_t const& h, std::string hist_name)
{
  TH1* data = get(h.at(DATA_SAMPLE_NAME), hist_name);
  TH1* mc   = get(h.at("MC"), hist_name);
  TAxis *x  = data->GetXaxis();
  TAxis *y  = data->GetYaxis();

  x->SetTitle("#beta_{reco}");
  y->SetRangeUser(0, 1.1 * std::max(data->GetMaximum(), mc->GetMaximum()));
  
  mc->SetLineColor(kRed);
  
  data->Draw();
  mc->Draw("SAMES");
}


void draw_dplane(hists_t const& h, std::string hist_name)
{
  TH1* data = get(h.at(DATA_SAMPLE_NAME), hist_name);
  TH1* mc   = get(h.at("MC"), hist_name);
  TAxis *x  = data->GetXaxis();
  TAxis *y  = data->GetYaxis();

  std::string axis_name = "#DeltaPlane (xz-View)";
  std::string canvas_name = "DPlane XZ";
  if (hist_name.find("yz") != std::string::npos) {
    axis_name = "#DeltaPlane (yz-View)";
    canvas_name = "DPlane YZ";
  }

  x->SetTitle(axis_name.c_str());
  x->SetRangeUser(0, 300);
  x->CenterTitle();
  
  y->SetTitle("Number of Tracks / 10 Planes");
  y->CenterTitle();
  y->SetRangeUser(0, 1.1 * std::max(data->GetMaximum(), mc->GetMaximum()));

  mc->SetLineColor(kRed);
  
  data->Draw("E HIST");



  TLegend *l = new TLegend(0.5, 0.7, 0.7, 0.88);
  l->SetTextSizePixels(20);
  l->AddEntry(data, "Data", "l");
  std::string mc_legend_title =
    "MC, #beta = " + MC_BETA_NICE_NAME;
  l->AddEntry(mc, mc_legend_title.c_str(), "l");
  l->Draw();  
}


void draw_n_hits(hists_t const& h, std::string hist_name)
{
  TH1* data = get(h.at(DATA_SAMPLE_NAME), hist_name);
  TH1* mc   = get(h.at("MC"), hist_name);
  TAxis *x  = data->GetXaxis();
  TAxis *y  = data->GetYaxis();

  std::string axis_name = "N^{hits}_{Total}";
  std::string canvas_name = "Number of Hits";
  if (hist_name.find("xz") != std::string::npos) {
    axis_name = "Number of Hits (xz-View)";
    canvas_name = "Number of Hits XZ";
  } else if (hist_name.find("yz") != std::string::npos) {
    axis_name = "Number of Hits (yz-View)";
    canvas_name = "Number of Hits YZ";
  }

  x->SetTitle(axis_name.c_str());
  x->CenterTitle();
  
  y->SetTitle("Number of Tracks / 10 Hits");
  y->CenterTitle();
  y->SetRangeUser(0, 1.1 * std::max(data->GetMaximum(), mc->GetMaximum()));

  mc->SetLineColor(kRed);
  

  data->Draw("E HIST");
  mc->Draw("SAMES");

  TLegend *l = new TLegend(0.6, 0.7, 0.7, 0.88);
  l->SetTextSizePixels(20);
  l->AddEntry(data, "Data", "l");
  l->AddEntry(mc, "MC", "l");
  l->Draw();  
}


void draw_n_tracks(hists_t const& h, std::string hist_name)
{
  TH1* data = get(h.at(DATA_SAMPLE_NAME), hist_name);
  TH1* mc   = get(h.at("MC"), hist_name);
  TAxis *x  = data->GetXaxis();
  TAxis *y  = data->GetYaxis();

  x->SetTitle(hist_name.c_str());
  x->CenterTitle();
  x->SetNdivisions(5, 0, 0);
  
  y->SetRangeUser(0, 1.1 * std::max(data->GetMaximum(), mc->GetMaximum()));

  mc->SetLineColor(kRed);

  data->Draw();
  mc->Draw("SAMES");
}


void draw_projected_length(hists_t const& h, std::string hist_name)
{
  TH1* data = get(h.at(DATA_SAMPLE_NAME), hist_name);
  TH1* mc   = get(h.at("MC"), hist_name);
  TAxis *x  = data->GetXaxis();
  TAxis *y  = data->GetYaxis();

  x->SetTitle((hist_name + " [m]").c_str());
  y->SetRangeUser(0, 1.1 * std::max(data->GetMaximum(), mc->GetMaximum()));

  mc->SetLineColor(kRed);
  
  data->Draw();

  mc->Draw("SAMES");
}


void draw_r2(hists_t const& h, std::string hist_name)
{
  TH1* data = get(h.at(DATA_SAMPLE_NAME), hist_name);
  TH1* mc   = get(h.at("MC"), hist_name);
  TAxis *x  = data->GetXaxis();
  TAxis *y  = data->GetYaxis();

  if(0 != access("r2histdat.txt", F_OK)){
    std::ofstream outfile("r2histdat.txt");
    if(!outfile.is_open()){
      fprintf(stderr, "Could not open r2histdat.txt for writing\n");
    }

    for(int i = 1; i <= data->GetNbinsX(); i++)
      outfile << "Data " << data->GetBinContent(i) << endl;
    for(int i = 1; i <= mc->GetNbinsX(); i++)
      outfile << "MC " << mc->GetBinContent(i) << endl;
    outfile.close();
  }

  std::string axis_name = "Correlation coefficient (xt)";
  std::string canvas_name = "r2_XT";
  if (hist_name.find("yt") != std::string::npos) {
    axis_name = "Correlation coefficient (yt)";
    canvas_name = "r2_YT";
  } else if (hist_name.find("min") != std::string::npos) {
    axis_name = "Correlation coefficient r^{2}_{min}";
    canvas_name = "r2_Min";

    if (hist_name.find("low gap") != std::string::npos)
      canvas_name = "r2_Min_Low_Gap";

    if (hist_name.find("low gap and slow") != std::string::npos)
      canvas_name = "r2_Min_Low_Gap_and_Slow";
  }
  
  x->SetTitle(axis_name.c_str());
  x->CenterTitle();
  
  y->SetTitle("Number of events");
  y->CenterTitle();
  y->SetRangeUser(5e-1, 1e5);

  x->SetTitleOffset(1.06);
  y->SetTitleOffset(1.1);
  x->SetLabelOffset(0.011);
  y->SetLabelOffset(0.015);
  
  mc->SetLineColor(kRed);
  mc->SetMarkerColor(kRed);
  mc->SetLineStyle(7);
  
  TCanvas * c1 = new TCanvas;

  c1->SetCanvasSize(600, 400);
  c1->SetRightMargin(0.025);
  c1->SetTopMargin(0.03);
  c1->SetLeftMargin(0.14);
  c1->SetBottomMargin(0.14);

  y->SetTickSize(0.015);
  x->SetDecimals();

  data->Draw("hist");
  mc->Draw("SAME hist");

  c1->SetLogy();
  c1->SetLogz();
  c1->SetTickx();
  c1->SetTicky();

  x->SetTitleSize(textsize);
  y->SetTitleSize(textsize);
  x->SetLabelSize(textsize);
  y->SetLabelSize(textsize);

  TLegend *leg = new TLegend(0.25, 0.75, 0.55, 0.93);
  leg->SetTextSize(textsize);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);

  leg->AddEntry(data, "Data", "l");
  std::string mc_legend_title =
    "MC, #beta = " + MC_BETA_NICE_NAME + ", normalized to data";
  leg->AddEntry(mc, mc_legend_title.c_str(), "l");
  leg->Draw();

  if (hist_name.find("min") != std::string::npos) {
    TLine *line = new TLine(0.95, 0.5, 0.95, 1e5);
    line->SetLineColor(kGreen + 2);
    line->SetLineWidth(2);
    line->SetLineStyle(2);
    line->Draw();

    TArrow *arrow = new TArrow(0.95, 3.5e4, 0.985, 3.5e4, 0.013, "|>");
    arrow->SetLineWidth(2);
    arrow->SetLineColor(kGreen + 2);
    arrow->SetFillColor(kGreen + 2);
    arrow->Draw();
  }

  c1->SaveAs("r2min-n-1.pdf");
}


void draw_theta(hists_t const& h, std::string hist_name)
{
  TH1* data = get(h.at(DATA_SAMPLE_NAME), hist_name);
  TH1* mc   = get(h.at("MC"), hist_name);
  TAxis *x  = data->GetXaxis();
  TAxis *y  = data->GetYaxis();

  x->SetTitle((hist_name).c_str());
  y->SetRangeUser(0, 1.1 * std::max(data->GetMaximum(), mc->GetMaximum()));

  mc->SetLineColor(kRed);
  
  std::string canvas_name = "Theta XZ";
  if (hist_name.find("yz") != std::string::npos) canvas_name = "Theta YZ";

  data->Draw();

  mc->Draw("SAMES");
}


void draw_time_gap(hists_t const& h, std::string hist_name)
{
  TH1* data = get(h.at(DATA_SAMPLE_NAME), hist_name);
  TH1* mc   = get(h.at("MC"), hist_name);
  TAxis *x  = data->GetXaxis();
  TAxis *y  = data->GetYaxis();

  if(0 != access("gaphistdat.txt", F_OK)){
    std::ofstream outfile("gaphistdat.txt");
    if(!outfile.is_open()){
      fprintf(stderr, "Could not open gaphistdat.txt for writing\n");
    }
    for(int i = 1; i <= data->GetNbinsX(); i++)
      outfile << "Data " << data->GetBinContent(i) << endl;
    for(int i = 1; i <= mc->GetNbinsX(); i++)
      outfile << "MC " << mc->GetBinContent(i) << endl;
    outfile.close();
  }

  std::string axis_name = "Time Gap Fraction (xt)";
  std::string canvas_name = "Time Gap XT";
  if (hist_name.find("yz") != std::string::npos) {
    axis_name = "Time Gap Fraction (yt)";
    canvas_name = "Time Gap YT";
  } else if (hist_name.find("max") != std::string::npos) {
    axis_name = "Time gap fraction f_{max}";
    canvas_name = "Time Gap Max";
  }

  x->SetTitle(axis_name.c_str());
  x->CenterTitle();
  
  y->SetTitle("Number of events");
  y->CenterTitle();
  auto y_max = 1.5;
  y->SetRangeUser(0, y_max);

  mc->SetLineColor(kRed);
  mc->SetMarkerColor(kRed);
  mc->SetLineStyle(7);
  
  TCanvas * c1 = new TCanvas;

  c1->SetCanvasSize(600, 400);
  c1->SetRightMargin(0.025);
  c1->SetTopMargin(0.03);
  c1->SetLeftMargin(0.14);
  c1->SetBottomMargin(0.14);

  c1->SetTickx();
  c1->SetTicky();

  x->SetTitleSize(textsize);
  y->SetTitleSize(textsize);
  x->SetLabelSize(textsize);
  y->SetLabelSize(textsize);
  y->SetTickSize(0.015);
  x->SetDecimals();
  x->SetTitleOffset(1.06);
  y->SetTitleOffset(1.1);
  x->SetLabelOffset(0.011);
  y->SetLabelOffset(0.015);

  data->Draw("hist");
  mc->Draw("SAME hist ][");

  TLegend *leg = new TLegend(0.35, 0.75, 0.65, 0.93);
  leg->SetTextSize(textsize);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSizePixels(20);
  leg->AddEntry(data, "Data", "l");
  std::string mc_legend_title =
    "MC, #beta = " + MC_BETA_NICE_NAME + ", normalized to data";
  leg->AddEntry(mc, mc_legend_title.c_str(), "l");
  leg->Draw();

  if (hist_name.find("max") != std::string::npos) {
    TLine *line = new TLine(0.2, 0, 0.2, y_max);
    line->SetLineColor(kGreen + 2);
    line->SetLineWidth(2);
    line->SetLineStyle(2);
    line->Draw();

    TArrow *arrow =
      new TArrow(0.2, 0.8 * y_max, 0.15, 0.8 * y_max, 0.012, "|>");
    arrow->SetLineWidth(2);
    arrow->SetLineColor(kGreen + 2);
    arrow->SetFillColor(kGreen + 2);
    arrow->Draw();
  }

  c1->SaveAs("fmax-n-1.pdf");
}


void draw_track_length(hists_t const& h)
{
  std::string hist_name = "Track Length";
  TH1* data = get(h.at(DATA_SAMPLE_NAME), hist_name);
  TH1* mc   = get(h.at("MC"), hist_name);
  TAxis *x  = data->GetXaxis();
  TAxis *y  = data->GetYaxis();

  x->SetTitle("Track Length [m]");
  x->CenterTitle();
  
  y->SetTitle("Number of Tracks / 2 m");
  y->CenterTitle();
  y->SetRangeUser(0, 1.1 * std::max(data->GetMaximum(), mc->GetMaximum()));
    
  mc->SetLineColor(kRed);

  data->Draw("E HIST");
  mc->Draw("SAMES");

  TLegend *l = new TLegend(0.6, 0.7, 0.7, 0.88);
  l->SetTextSizePixels(20);
  l->AddEntry(data, "Data", "l");
  l->AddEntry(mc, "MC", "l");
  l->Draw();
}
