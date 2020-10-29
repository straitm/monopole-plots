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
  y->SetRangeUser(5e-1, 3e4);

  x->SetTitleOffset(1.06);
  y->SetTitleOffset(1.1);
  x->SetLabelOffset(0.011);
  y->SetLabelOffset(0.015);
  
  TCanvas * c1 = new TCanvas;

  c1->SetCanvasSize(600, 400);
  c1->SetRightMargin(0.025);
  c1->SetTopMargin(0.03);
  c1->SetLeftMargin(0.14);
  c1->SetBottomMargin(0.14);

  y->SetTickSize(0.015);
  x->SetDecimals();

  c1->SetLogy();
  c1->SetLogz();
  c1->SetTickx();
  c1->SetTicky();

  x->SetTitleSize(textsize);
  y->SetTitleSize(textsize);
  x->SetLabelSize(textsize);
  y->SetLabelSize(textsize);

  TH1D * mc5e3 = new TH1D("mc5e3", "", 40, 0, 1);

  mc5e3->SetLineColor(kBlue);
  mc5e3->SetLineColor(kBlue);

  mc5e3->SetBinContent(1,38);
  mc5e3->SetBinContent(2,8);
  mc5e3->SetBinContent(3,4);
  mc5e3->SetBinContent(4,4);
  mc5e3->SetBinContent(5,3);
  mc5e3->SetBinContent(6,7);
  mc5e3->SetBinContent(7,2);
  mc5e3->SetBinContent(9,2);
  mc5e3->SetBinContent(10,1);
  mc5e3->SetBinContent(12,1);
  mc5e3->SetBinContent(13,3);
  mc5e3->SetBinContent(15,2);
  mc5e3->SetBinContent(16,2);
  mc5e3->SetBinContent(17,1);
  mc5e3->SetBinContent(18,1);
  mc5e3->SetBinContent(19,1);
  mc5e3->SetBinContent(20,1);
  mc5e3->SetBinContent(21,3);
  mc5e3->SetBinContent(22,2);
  mc5e3->SetBinContent(23,1);
  mc5e3->SetBinContent(24,3);
  mc5e3->SetBinContent(25,2);
  mc5e3->SetBinContent(26,6);
  mc5e3->SetBinContent(27,4);
  mc5e3->SetBinContent(28,4);
  mc5e3->SetBinContent(29,6);
  mc5e3->SetBinContent(30,4);
  mc5e3->SetBinContent(31,5);
  mc5e3->SetBinContent(32,3);
  mc5e3->SetBinContent(33,5);
  mc5e3->SetBinContent(34,10);
  mc5e3->SetBinContent(35,5);
  mc5e3->SetBinContent(36,18);
  mc5e3->SetBinContent(37,18);
  mc5e3->SetBinContent(38,37);
  mc5e3->SetBinContent(39,61);
  mc5e3->SetBinContent(40,5176);


  mc5e3->Scale(data->Integral()/mc5e3->Integral());

  mc->SetLineColor(kRed);
  mc->SetMarkerColor(kRed);
  mc->SetLineStyle(9);

  mc5e3->SetLineColor(kBlue);
  mc5e3->SetMarkerColor(kBlue);
  mc5e3->SetLineStyle(7);

  data ->SetLineWidth(2);
  mc   ->SetLineWidth(2);
  mc5e3->SetLineWidth(2);

  
  data->Draw("hist");
  mc   ->Draw("SAME hist");
  mc5e3->Draw("SAME hist");

  TLegend *leg = new TLegend(0.25, 0.93 - 3*(1.1*textsize), 0.53, 0.93);
  leg->SetTextSize(textsize);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);

  leg->AddEntry(data, "Data", "l");
  leg->AddEntry(mc,    "MC, #beta = 1#kern[-0.1]{ }#times10^{#minus3}, normalized to data", "l");
  leg->AddEntry(mc5e3, "MC, #beta = 5#kern[-0.5]{ }#times10^{#minus3}", "l");
  leg->Draw();

  if (hist_name.find("min") != std::string::npos) {
    TLine *line = new TLine(0.95, 0.5, 0.95, 3e4);
    line->SetLineColor(kGreen + 2);
    line->SetLineWidth(2);
    line->SetLineStyle(2);
    line->Draw();

    TArrow *arrow = new TArrow(0.95, 1.5e4, 0.985, 1.5e4, 0.013, "|>");
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

  TH1D * mc5e3 = new TH1D("mc5e3", "", 100, 0, 1);
  TH1D * mc3e4 = new TH1D("mc3e4", "", 100, 0, 1);

  mc->SetLineColor(kRed);
  mc->SetMarkerColor(kRed);
  mc->SetLineStyle(9);
  
  mc5e3->SetMarkerColor(kBlue);
  mc5e3->SetLineColor(kBlue);
  mc5e3->SetLineStyle(7);

  mc3e4->SetMarkerColor(kBlack);
  mc3e4->SetLineColor(kBlack);
  mc3e4->SetLineStyle(kDashed);
  
  mc5e3->SetBinContent(2,41);
  mc5e3->SetBinContent(3,218);
  mc5e3->SetBinContent(4,358);
  mc5e3->SetBinContent(5,405);
  mc5e3->SetBinContent(6,493);
  mc5e3->SetBinContent(7,543);
  mc5e3->SetBinContent(8,612);
  mc5e3->SetBinContent(9,486);
  mc5e3->SetBinContent(10,440);
  mc5e3->SetBinContent(11,313);
  mc5e3->SetBinContent(12,289);
  mc5e3->SetBinContent(13,252);
  mc5e3->SetBinContent(14,208);
  mc5e3->SetBinContent(15,137);
  mc5e3->SetBinContent(16,137);
  mc5e3->SetBinContent(17,92);
  mc5e3->SetBinContent(18,91);
  mc5e3->SetBinContent(19,66);
  mc5e3->SetBinContent(20,56);
  mc5e3->SetBinContent(21,34);
  mc5e3->SetBinContent(22,34);
  mc5e3->SetBinContent(23,31);
  mc5e3->SetBinContent(24,20);
  mc5e3->SetBinContent(25,13);
  mc5e3->SetBinContent(26,14);
  mc5e3->SetBinContent(27,17);
  mc5e3->SetBinContent(28,9);
  mc5e3->SetBinContent(29,10);
  mc5e3->SetBinContent(30,9);
  mc5e3->SetBinContent(31,8);
  mc5e3->SetBinContent(32,8);
  mc5e3->SetBinContent(33,4);
  mc5e3->SetBinContent(34,5);
  mc5e3->SetBinContent(36,2);
  mc5e3->SetBinContent(37,2);
  mc5e3->SetBinContent(38,3);
  mc5e3->SetBinContent(39,3);
  mc5e3->SetBinContent(41,1);
  mc5e3->SetBinContent(45,1);
  mc5e3->SetBinContent(46,1);

  mc3e4->SetBinContent(2,64);
  mc3e4->SetBinContent(3,427);
  mc3e4->SetBinContent(4,844);
  mc3e4->SetBinContent(5,1017);
  mc3e4->SetBinContent(6,869);
  mc3e4->SetBinContent(7,545);
  mc3e4->SetBinContent(8,348);
  mc3e4->SetBinContent(9,214);
  mc3e4->SetBinContent(10,109);
  mc3e4->SetBinContent(11,81);
  mc3e4->SetBinContent(12,30);
  mc3e4->SetBinContent(13,20);
  mc3e4->SetBinContent(14,7);
  mc3e4->SetBinContent(15,2);
  mc3e4->SetBinContent(16,5);
  mc3e4->SetBinContent(17,3);
  mc3e4->SetBinContent(18,4);
  mc3e4->SetBinContent(19,1);
  mc3e4->SetBinContent(21,1);


  mc5e3->Rebin(2);
  mc5e3->Scale(data->Integral()/mc5e3->Integral());
  mc5e3->Draw("hist same ][");

  mc3e4->Rebin(2);
  mc3e4->Scale(data->Integral()/mc3e4->Integral());
  //mc3e4->Draw("hist same ][");

  mc->SetLineWidth(2);
  mc5e3->SetLineWidth(2);
  mc3e4->SetLineWidth(2);
  data->SetLineWidth(2);

  TLegend *leg = new TLegend(0.33, 0.93-3*(1.1*textsize), 0.61, 0.93);
  leg->SetTextSize(textsize);
  leg->SetBorderSize(0);
  leg->SetFillStyle(0);
  leg->SetTextFont(42);
  leg->SetTextSizePixels(20);
  leg->AddEntry(data, "Data", "l");
  leg->AddEntry(mc,    "MC, #beta = 1#kern[-0.1]{ }#times10^{#minus3}, normalized to data", "l");
  leg->AddEntry(mc5e3, "MC, #beta = 5#kern[-0.5]{ }#times10^{#minus3}", "l");
  //leg->AddEntry(mc3e4, "MC, #beta = 3#kern[-0.5]{ }#times10^{#minus4}", "l");
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


  // This is to provide to the PRD referee, who asked for the 1D
  // fmax plot without the r2 cut.

  delete mc5e3;

  data->SetBinContent(1,2210);
  data->SetBinContent(2,210);
  data->SetBinContent(3,191);
  data->SetBinContent(4,383);
  data->SetBinContent(5,770);
  data->SetBinContent(6,1321);
  data->SetBinContent(7,2077);
  data->SetBinContent(8,2761);
  data->SetBinContent(9,3442);
  data->SetBinContent(10,4064);
  data->SetBinContent(11,4537);
  data->SetBinContent(12,5315);
  data->SetBinContent(13,5776);
  data->SetBinContent(14,6191);
  data->SetBinContent(15,6491);
  data->SetBinContent(16,6552);
  data->SetBinContent(17,6897);
  data->SetBinContent(18,6782);
  data->SetBinContent(19,6897);
  data->SetBinContent(20,6601);
  data->SetBinContent(21,6574);
  data->SetBinContent(22,6542);
  data->SetBinContent(23,6291);
  data->SetBinContent(24,6221);
  data->SetBinContent(25,6170);
  data->SetBinContent(26,6170);
  data->SetBinContent(27,5586);
  data->SetBinContent(28,5284);
  data->SetBinContent(29,5059);
  data->SetBinContent(30,4583);
  data->SetBinContent(31,4449);
  data->SetBinContent(32,4167);
  data->SetBinContent(33,4017);
  data->SetBinContent(34,3670);
  data->SetBinContent(35,3395);
  data->SetBinContent(36,3396);
  data->SetBinContent(37,3165);
  data->SetBinContent(38,2975);
  data->SetBinContent(39,2872);
  data->SetBinContent(40,2781);
  data->SetBinContent(41,2638);
  data->SetBinContent(42,2643);
  data->SetBinContent(43,2744);
  data->SetBinContent(44,2832);
  data->SetBinContent(45,3015);
  data->SetBinContent(46,3494);
  data->SetBinContent(47,4610);
  data->SetBinContent(48,5459);
  data->SetBinContent(49,6425);
  data->SetBinContent(50,6361);

  mc->SetBinContent(1,38566.89);
  mc->SetBinContent(2,105115);
  mc->SetBinContent(3,38173.11);
  mc->SetBinContent(4,12925.12);
  mc->SetBinContent(5,5929.804);
  mc->SetBinContent(6,2617.452);
  mc->SetBinContent(7,1227.655);
  mc->SetBinContent(8,903.3686);
  mc->SetBinContent(9,694.8989);
  mc->SetBinContent(10,463.2659);
  mc->SetBinContent(11,393.776);
  mc->SetBinContent(12,324.2861);
  mc->SetBinContent(13,138.9798);
  mc->SetBinContent(14,69.48989);
  mc->SetBinContent(15,115.8165);
  mc->SetBinContent(16,162.1431);
  mc->SetBinContent(17,69.48989);
  mc->SetBinContent(18,648.5723);
  mc->SetBinContent(19,162.1431);
  mc->SetBinContent(20,138.9798);
  mc->SetBinContent(21,162.1431);
  mc->SetBinContent(22,324.2861);
  mc->SetBinContent(23,231.633);
  mc->SetBinContent(24,671.7356);
  mc->SetBinContent(26,92.65319);
  mc->SetBinContent(27,162.1431);
  mc->SetBinContent(28,69.48989);
  mc->SetBinContent(29,393.776);
  mc->SetBinContent(30,92.65319);
  mc->SetBinContent(31,46.32659);
  mc->SetBinContent(33,208.4697);
  mc->SetBinContent(34,138.9798);
  mc->SetBinContent(35,69.48989);
  mc->SetBinContent(36,23.1633);
  mc->SetBinContent(37,69.48989);
  mc->SetBinContent(38,185.3064);
  mc->SetBinContent(39,46.32659);
  mc->SetBinContent(40,162.1431);
  mc->SetBinContent(43,208.4697);
  mc->SetBinContent(44,69.48989);
  mc->SetBinContent(46,115.8165);
  mc->SetBinContent(47,208.4697);
  mc->SetBinContent(48,254.7963);
  mc->SetBinContent(49,69.48989);
  mc->SetBinContent(50,138.9798);

  mc->Scale(data->Integral()/mc->Integral());

  mc->Scale(1e-3);
  data->Scale(1e-3);

  data->Draw("hist ][");
  mc->Draw("same hist ][");

  y_max = mc->GetMaximum()*1.1;

  {
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

    TLegend *leg = new TLegend(0.33, 0.93-3*(1.1*textsize), 0.61, 0.93);
    leg->SetTextSize(textsize);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSizePixels(20);
    leg->AddEntry((TH1D*)NULL, "No r^{2}_{min} cut", "");
    leg->AddEntry(data, "Data", "l");
    leg->AddEntry(mc,    "MC, #beta = 1#kern[-0.1]{ }#times10^{#minus3}, normalized to data", "l");
    //leg->AddEntry(mc3e4, "MC, #beta = 3#kern[-0.5]{ }#times10^{#minus4}", "l");
    leg->Draw();

  }

  y->SetRangeUser(0, mc->GetMaximum()*1.1);
  y->SetTitle("Number of events #times 10^{3}");

  c1->SaveAs("fmax-n-2.pdf");
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
