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

enum whatplot_t { plotr2, plotfmax, plotNOTHING };

static const double textsize = 0.057;
static const int datacolor = kGray+1;
static const int mccolor = kRed+2;

typedef std::map<std::string, TGraph*> graph_t;
typedef std::map<std::string, graph_t> graphs_t;

using namespace MFRoot;

#define USECACHE

void fill_graphs(graph_t & g, tree_t const& t, TTree* tree,
                 const std::string & type, const whatplot_t plotwhat)
{
  std::cout << "Filling graphs for type: " << type << std::endl;

  std::string *input_file_name = new std::string("invalid");
  if (type != DATA_SAMPLE_NAME)
    tree->SetBranchAddress("input_file_name", &input_file_name);

  #ifdef USECACHE
  if(type == "Data"){
    if(plotwhat == plotr2){
      std::ifstream infilenose("nose_data.txt");
      if(infilenose.is_open()){
        printf("Reading nose data from text file and not filling anything else\n");

        double x, y;
        while(infilenose >> x >> y){
          TGraph * G = g.at("Beta_vs_r2_Low_Gap");
          G->SetPoint(G->GetN(), x, y);
        }
        return;
      }
    }

    if(plotwhat == plotfmax){
      std::ifstream infilemouth("mouth_data.txt");
      if(infilemouth.is_open()){
        printf("Reading mouth data from text file and not filling anything else\n");

        double x, y;
        while(infilemouth >> x >> y){
          TGraph * G = g.at("Beta_vs_Gap_Low_r2");
          G->SetPoint(G->GetN(), x, y);
        }
        return;
      }
    }
  }
  #endif

  Event_List elist(type);
  for(int64_t entry = 0; entry < tree->GetEntries(); ++entry){
    if(entry % 100000 == 0) printf("processing entry %ld\n", entry);

    tree->GetEntry(entry);

    Event_Info e(t, *input_file_name);

    if (!e.beta_value_matches(MC_BETA_NAME)) continue;
    if(elist.is_duplicate(e.run(), e.subrun(), e.event())) continue;
    if(!e.is_preselected_reco()) continue;

    // Limit the number of MC points to keep the plot uncluttered
    if ((type != "MC" || entry < 10000) && e.low_gap_max()){
      TGraph * G = g.at("Beta_vs_r2_Low_Gap");
      G->SetPoint(G->GetN(), e.r2_min(), e.beta());
    }

    if(e.high_r2_min()){
      TGraph * G = g.at("Beta_vs_Gap_Low_r2");
      G->SetPoint(G->GetN(), e.gap_max(), e.beta());
    }
  }

  delete input_file_name;
}

static double gdist(const double deltar2, const double deltalog10beta)
{
  // divide deltalog10beta by 5 because the range is from 10**-5
  // to 1, except the y axis is shorter, so divide by 4.
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

      // Distance in units of fraction of plot width
      const double dist = gdist(r2-nr2, log10beta - nlog10beta);
      const bool nearedge = r2 < 0.01;

      // Have to allow for closer points at the edge so it is visually filled
      // in where it needs to be.  Otherwise, it makes sense for the spacing
      // to be just somewhat smaller than the marker diameter.
      if((!nearedge && dist < 0.005) || (nearedge && dist < 0.003)){
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
}

void save_data_graph(TGraph * g, const char * const facepart)
{
  std::ofstream out(Form("%s_data.txt", facepart));

  if(!out.is_open()){
    fprintf(stderr, "Could not open %s_data.txt for writing\n", facepart);
    return; 
  }

  for(int i = 0; i < g->GetN(); i++)
    out << g->GetX()[i] << " " << g->GetY()[i] << "\n";
}

void draw_beta_fmax(graphs_t const& g, const std::string & name)
{
  TCanvas * c1 = new TCanvas;

  const double rat = 400/267.2;
  
  c1->SetCanvasSize(600, 400/rat);
  c1->SetRightMargin(0.025);
  c1->SetTopMargin(0.03 * rat);
  c1->SetLeftMargin(0.14);
  c1->SetBottomMargin(0.14 * rat);

  c1->SetLogy();
  c1->SetLogz();
  c1->SetTickx();
  c1->SetTicky();

  const double ymin = 1e-4;
  
  TH2D dum("dum", "", 1, 0, 1, 1, ymin, 1e-1);

  TGraph* data = dynamic_cast<TGraph *>(g.at("Data").at(name));
  if(data == NULL){
    fprintf(stderr, "data graph doesn't exist\n");
    _exit(1);
  }

  save_data_graph(data, "mouth");

  TGraph* mc = g.at("MC").at(name);

  TAxis* x = dum.GetXaxis();
  TAxis* y = dum.GetYaxis();

  x->SetTitle("Largest gap, f_{max}");
  y->SetTitle("Reconstructed speed (#beta)  ");

  x->SetTitleSize(textsize*rat);
  x->SetLabelSize(textsize*rat);
  y->SetTitleSize(textsize*rat);
  y->SetLabelSize(textsize*rat);
  x->CenterTitle();
  y->CenterTitle();

  y->SetTitleOffset(1.08/rat);
  x->SetTitleOffset(1.04);

  data->SetMarkerStyle(kFullCircle);
  data->SetMarkerColor(datacolor);

  mc->SetMarkerStyle(kFullSquare);
  mc->SetMarkerColor(mccolor);
  mc->SetMarkerSize(0.9);

  dum.Draw();
  data->Draw("P");
  mc->Draw("p");

  TLine *line = new TLine(0.2, ymin, 0.2, 1e-2);
  line->SetLineColor(kGreen + 2);
  line->SetLineWidth(3);
  line->SetLineStyle(2);
  line->Draw();

  TLine *line2 = new TLine(0.2, 1e-2, 0, 1e-2);
  line2->SetLineColor(kGreen + 2);
  line2->SetLineWidth(3);
  line2->SetLineStyle(2);
  line2->Draw();

  c1->SaveAs("scatterfmax.pdf");
}

void draw_beta_r2(graphs_t const& g, const std::string & name)
{
  TCanvas * c1 = new TCanvas;
  
  c1->SetCanvasSize(600, 400);
  c1->SetRightMargin(0.115);
  c1->SetTopMargin(0.03);
  c1->SetLeftMargin(0.14);
  c1->SetBottomMargin(0.14);

  c1->SetLogy();
  c1->SetLogz();
  c1->SetTickx();
  c1->SetTicky();
  
  const double ylow = 1e-5;
  const double yhigh = 1;

  const unsigned int nxbins = 40, nybins = 14;
  double ybins[nybins+1] = { 0 };
  for(unsigned int i = 0; i <= nybins; i++)
    ybins[i] = exp(log(ylow) + double(i)/nybins * (log(yhigh) - log(ylow)));
  
  TH2D heatmc("heatmc", "", nxbins, 0, 1, nybins, ybins);

  TGraph* data = dynamic_cast<TGraph *>(g.at("Data").at(name));
  if(data == NULL){
    fprintf(stderr, "data graph doesn't exist\n");
    _exit(1);
  }

  TGraph* mc = g.at("MC").at(name);

  TAxis* x = heatmc.GetXaxis();
  TAxis* y = heatmc.GetYaxis();
  TAxis* z = heatmc.GetZaxis();

  x->SetTitle("Correlation coefficient, r_{min}^{2}");
  y->SetTitle("Reconstructed speed (#beta)");
#if 0
  z->SetTitle("MC density");
#endif

  x->CenterTitle();
  y->CenterTitle();
  z->CenterTitle();
  x->SetTickSize(0.025);
  y->SetTickSize(0.025);
  x->SetTitleSize(textsize);
  x->SetLabelSize(textsize);
  y->SetTitleSize(textsize);
  y->SetLabelSize(textsize);
  z->SetTitleSize(textsize);
  z->SetLabelSize(textsize);

  z->SetRangeUser(9.999, 1000);

  heatmc.SetContour(8);

  y->SetTitleOffset(1.08);
  x->SetTitleOffset(1.04);
  z->SetLabelOffset(0.000);
  
  data->SetMarkerStyle(kFullCircle);
  data->SetMarkerColor(datacolor);

  mc->SetMarkerStyle(kFullSquare);
  mc->SetMarkerColor(mccolor);
  mc->SetMarkerSize(0.9);

  optimize_graph(data);

  save_data_graph(data, "nose");

  for(int i = 0; i < mc->GetN(); i++)
    heatmc.Fill(mc->GetX()[i], mc->GetY()[i]);

  heatmc.SavePrimitive(std::cout);

  heatmc.Draw("colz");
  data->Draw("P");
  mc->Draw("p");

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

  c1->RedrawAxis();

  c1->SaveAs("scatterr2.pdf");
}

int main(int argc, char ** argv)
{
  gStyle->SetOptStat(0);

  if(argc != 2){
    fprintf(stderr, "Tell me either r2 or fmax\n");
    exit(1);
  }

  const enum whatplot_t plotwhat = !strcmp(argv[1], "r2")?plotr2
                                  :!strcmp(argv[1], "fmax")?plotfmax
                                  :plotNOTHING;
  if(plotwhat == plotNOTHING){
    fprintf(stderr, "Eh, sonny?  Tell me either r2 or fmax\n");
    exit(1);
  }

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
  for (auto const& tree : trees) {
    std::string type = tree.first;
    (g[type])["Beta_vs_r2_Low_Gap"] = new TGraph;
    (g[type])["Beta_vs_Gap_Low_r2"] = new TGraph;
  }
  
  for (auto const& tree : trees) {
    const std::string type = tree.first;
    fill_graphs(g.at(type), t.at(type), tree.second, type, plotwhat);
  }

  if     (plotwhat == plotr2  ) draw_beta_r2(g,  "Beta_vs_r2_Low_Gap");
  else if(plotwhat == plotfmax) draw_beta_fmax(g,"Beta_vs_Gap_Low_r2");
}
