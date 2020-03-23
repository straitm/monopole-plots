/*
  Run over Data and MC monopole trees and generate 2D scatter plots/histograms.
*/

#include "Constants.hh"
#include "Event_Info.hh"
#include "Event_List.hh"
#include "MFRoot.hh"

#include <TArrow.h>
#include <TFile.h>
#include <TGraph.h>
#include <TH2.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TStyle.h>
#include <TTree.h>

#include <limits.h>

#include <fstream>
#include <iostream>
#include <map>

static double textsize = 0.057;

typedef std::map<std::string, TGraph*> graph_t;
typedef std::map<std::string, graph_t> graphs_t;

using namespace MFRoot;

void fill_graphs(graph_t & g, tree_t const& t, TTree* tree, std::string type)
{
  std::cout << "Filling graphs for type: " << type << std::endl;

  std::string *input_file_name = new std::string("invalid");
  if (type != DATA_SAMPLE_NAME)
    tree->SetBranchAddress("input_file_name", &input_file_name);

  const int maxentries = type == "MC"? 10000: INT_MAX;

  if(type == "Data"){
    std::ifstream infile("nose_data.txt");
    if(infile.is_open()){
      printf("Reading nose data from text file and not filling anything else\n");

      double x, y;
      while(infile >> x >> y){
        TGraph * G = g.at("Beta_vs_r2_Low_Gap");
        G->SetPoint(G->GetN(), x, y);
      }
      return;
    }
  }

  Event_List elist(type);
  for (int entry = 0; entry < maxentries && entry < tree->GetEntries(); ++entry)
  {
    if (entry % 100000 == 0)
      std::cout << "processing entry " << entry << " ..." << std::endl;

    tree->GetEntry(entry);
    Event_Info e(t, *input_file_name);

    if (!e.beta_value_matches(MC_BETA_NAME)) continue;
    
    bool duplicate_event = elist.is_duplicate(e.run(), e.subrun(), e.event());
    if(duplicate_event) continue;
    if(!e.is_preselected_reco()) continue;
    g.at("Gap_vs_r2_Preselected")->
      SetPoint(g.at("Gap_vs_r2_Preselected")->GetN(), e.r2_min(), e.gap_max());
    g.at("Beta_vs_r2_Preselected")->
      SetPoint(g.at("Beta_vs_r2_Preselected")->GetN(), e.r2_min(), e.beta());

    if (e.low_gap_max()){
      TGraph * G = g.at("Beta_vs_r2_Low_Gap");
      G->SetPoint(G->GetN(), e.r2_min(), e.beta());
    }

    // This will only work for MC graph, data graph will be empty.
    if (e.high_mc_fraction())
      g.at("Good_Reco")->SetPoint(g.at("Good_Reco")->GetN(), e.r2_min(), e.gap_max());

    if (e.is_linear())
      g.at("Signal")->SetPoint(g.at("Signal")->GetN(), e.r2_min(), e.gap_max());
  }

  delete input_file_name;
}

static double gdist(const double deltar2, const double deltalog10beta)
{
  return hypot(deltar2, deltalog10beta/4.);
}

// Remove similar points from the input graph, so that it can be drawn
// faster, but look the same.  This is *only* for plotting.  If you try
// to do something statistical with this graph later, it will be super wrong.
static void optimize_graph(TGraph * & g)
{
  TGraph * newg = (TGraph *)g->Clone("newg");

  // Is there a better way?
  while(newg->GetN()) newg->RemovePoint(newg->GetN()-1);

  for(int i = 0; i < g->GetN(); i++){
    const double r2 = g->GetX()[i];
    const double log10beta = log10(g->GetY()[i]);

    // This looks like a horrible O(N^2) algorithm, but since the point
    // of this is to minimize the number of points in newg, it is only ~O(N).
    bool pass = true;
    for(int j = 0; j < newg->GetN(); j++){
      const double nr2 = newg->GetX()[j];
      const double nlog10beta = log10(newg->GetY()[j]);

      const double dist = gdist(r2-nr2, log10beta - nlog10beta);
      const bool nearedge = r2 < 0.01;

      // Have to allow for closer points at the edge so it is visually filled
      // in where it needs to be.
      if((!nearedge && dist < 0.007) || (nearedge && dist < 0.001)){
        pass = false;
        break;
      }
    }
    if(pass) newg->SetPoint(newg->GetN(), g->GetX()[i], g->GetY()[i]);
  }

  printf("Reduced graph from %u to %u points\n", g->GetN(), newg->GetN());

  delete g;

  g = newg;
}


void draw_gap_r2(graphs_t const& g, std::string name)
{
  // For some reason, SetRangeUser does not work, so I am creating this
  // two-point grap to set the range.
  TGraph* range = new TGraph(2);
  range->SetPoint(0, 0, 0);
  range->SetPoint(1, 1, 1);

  TGraph* data = g.at("Data").at(name);
  TGraph* mc   = g.at("MC").at(name);
  TAxis* x     = range->GetXaxis();
  TAxis* y     = range->GetYaxis();

  x->SetTitle("r^{2}_{min}");
  x->CenterTitle();

  y->SetTitle("Time Gap Max");
  y->SetTitleOffset(1.2);
  y->CenterTitle();

  range->SetTitle("");
  range->SetMarkerColor(kWhite);
  
  data->SetMarkerStyle(2);
  data->SetMarkerSize(0.5);

  mc->SetMarkerStyle(4);
  mc->SetMarkerColor(kRed);
  mc->SetMarkerSize(0.5);
  
#if 0
  cans.new_canvas("Scatter Overlay - " + name, false, false, false, false, 2);
  range->Draw("AP");
  mc->Draw("SAME P");
  data->Draw("SAME P");
#endif

}

void save_data_graph(TGraph * g)
{
  std::ofstream out("nose_data.txt");

  if(!out.is_open()){
    fprintf(stderr, "Could not open nose_data.txt for writing\n");
    return; 
  }

  for(int i = 0; i < g->GetN(); i++)
    out << g->GetX()[i] << " " << g->GetY()[i] << "\n";
}

// The "nose plot"
void draw_beta_r2(graphs_t const& g, std::string name)
{
  TCanvas * c1 = new TCanvas;
  
  c1->SetCanvasSize(600, 400);
  c1->SetRightMargin(0.025);
  c1->SetTopMargin(0.03);
  c1->SetLeftMargin(0.14);
  c1->SetBottomMargin(0.14);

  c1->SetLogy();
  c1->SetLogz();
  c1->SetTickx();
  c1->SetTicky();
  
  // For some reason, SetRangeUser does not work, so I am creating this
  // two-point grap to set the range.
  TGraph* range = new TGraph(2);
  range->SetPoint(0, 0, 0);
  range->SetPoint(1, 1, 1);

  TGraph* data = dynamic_cast<TGraph *>(g.at("Data").at(name));
  if(data == NULL){
    fprintf(stderr, "data graph doesn't exist\n");
    _exit(1);
  }

  TGraph* mc   = g.at("MC").at(name);

  // TH2F *mc_heat = make_heat_map();
  const int n_bins = 19;
  double ybins[n_bins + 1];
  
  double low_bin = -5.0;
  double high_bin = -0.25;
  double delta_bin = (high_bin - low_bin) / n_bins;

  for (int bin = 0; bin != n_bins + 1; ++bin)
  {
    double exponent = low_bin + bin * delta_bin + delta_bin / 2;
    ybins[bin] = std::pow(10, exponent);
  }
  
  TH2F* mc_heat = new TH2F("mc_heat", "", 100, 0, 1.0001, n_bins, ybins);
  /*
    TGraph::GetN() usually returns the number of points in the TGraph object.
    However, the way that we implemented it above, TGraph does not know
    how many points it contains, so it simply returns TGraph::GetMaxSize().
    There are two ways to fix this:
      1. Fill the TH2F object above without using TGraph for MC.
      2. Use TGraph::Set() to tell it how many points it has once it 
         has been filled.
  */
  int n_graph_points = mc->GetN();
  for (int point = 0; point != n_graph_points; ++point)
  {
    double x, y;
    mc->GetPoint(point, x, y);
    mc_heat->Fill(x, y);
  }

  // only show the bins that have at least 10 entries
  mc_heat->SetMinimum(10);
  
  TAxis* x     = range->GetXaxis();
  TAxis* y     = range->GetYaxis();

  x->SetTitle("Correlation coefficient, r_{min}^{2}");
  x->CenterTitle();
  x->SetRangeUser(0, 1);

  x->SetTitleSize(textsize);
  x->SetLabelSize(textsize);
  y->SetTitleSize(textsize);
  y->SetLabelSize(textsize);

  y->SetTitle("Velocity (#beta)");
  y->SetTitleOffset(1.08);
  y->CenterTitle();
  y->SetRangeUser(1e-5, 1);

  range->SetTitle("");
  range->SetMarkerColor(kWhite);
  
  data->SetMarkerStyle(kFullCircle);
  data->SetMarkerColor(kGray+1);

  mc->SetMarkerStyle(kFullSquare);
  mc->SetMarkerColor(kBlack);
  mc->SetMarkerSize(0.9);

  optimize_graph(data);

  save_data_graph(data);

  std::string canvas_name_overlay = "Scatter_Overlay_" + name;
  range->Draw("AP");
  data->Draw("P");
  mc->Draw("p");

  if (name.find("Low_Gap") != std::string::npos)
  {
    TLine *line = new TLine(0.95, 1e-5, 0.95, 1e-2);
    line->SetLineColor(kGreen + 2);
    line->SetLineWidth(3);
    line->SetLineStyle(2);
    line->Draw();

    TLine *line2 = new TLine(0.95, 1e-2, 1, 1e-2);
    line2->SetLineColor(kGreen + 2);
    line2->SetLineWidth(3);
    line2->SetLineStyle(2);
    line2->Draw();
  }

  c1->SaveAs("scatter.pdf");
}

TH2F* make_heat_map()
{
  const int n_bins = 19;
  double ybins[n_bins + 1];
  
  double low_bin = -5.0;
  double high_bin = -0.25;
  double delta_bin = (high_bin - low_bin) / n_bins;

  for (int bin = 0; bin != n_bins + 1; ++bin)
  {
    double exponent = low_bin + bin * delta_bin + delta_bin / 2;
    ybins[bin] = std::pow(10, exponent);
  }
  
  return new TH2F("mc_heat", "", 40, 0, 1.0001, n_bins + 1, ybins);
}


int main()
{
  gStyle->SetOptStat("nemroui");

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

  graphs_t g;
  for (auto const& tree : trees)
  {
    std::string type = tree.first;
    (g[type])["Gap_vs_r2_Preselected"] = new TGraph;
    (g[type])["Beta_vs_r2_Preselected"] = new TGraph;
    (g[type])["Beta_vs_r2_Low_Gap"] = new TGraph;
    (g[type])["Signal"] = new TGraph;
    (g[type])["Good_Reco"] = new TGraph;
  }
  
  for (auto const& tree : trees)
  {
    std::string type = tree.first;
    fill_graphs(g.at(type), t.at(type), tree.second, type);
  }

  draw_beta_r2(g, "Beta_vs_r2_Low_Gap");
}
