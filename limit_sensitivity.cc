/*
  Draw a plot of monopole limits and our sensitivity or limits.
*/

#include <TGraph.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TLegend.h>
#include <TTree.h>
#include <TStyle.h>
#include <TArrow.h>

#include <iostream>
#include <fstream>
#include <map>
#include <string>

//#define THIRTEENYEAR

static double textsize = 0.052;
static int font = 42;

double crossover(TGraph * a, TGraph * b)
{
  const int n = 1000;
  const double minb = -3, maxb = -1;
  for(int i = 0; i <= n; i++){
    const double beta = pow(10, minb + (maxb-minb)*double(i)/n);
    if(b->Eval(beta) < a->Eval(beta)) return beta;
  }
  fprintf(stderr, "Did not find crossover, oh no\n");
  return 0;
}

int main()
{
  gStyle->SetOptStat(0);

  std::map<std::string, TGraph*> nova;
  nova["halffast"] = new TGraph();
  nova["fullfast"] = new TGraph();
  nova["halfslow"] = new TGraph();
  nova["fullslow"] = new TGraph();
  nova["halfboth"] = new TGraph();
  nova["fullboth"] = new TGraph();

  {
    std::ifstream infiles("limitsensitivityfastdata.txt");
    std::ifstream infilef("limitsensitivityslowdata.txt");

    if(infiles.is_open()){
      double logbeta, half, full;
      while(infiles >> logbeta >> half >> full){
        nova["halffast"]->SetPoint(nova["halffast"]->GetN(), pow(10, logbeta), half);
        nova["fullfast"]->SetPoint(nova["fullfast"]->GetN(), pow(10, logbeta), full);
      }
    }

    if(infilef.is_open()){
      double logbeta, half, full;
      while(infilef >> logbeta >> half >> full){
        nova["halfslow"]->SetPoint(nova["halfslow"]->GetN(), pow(10, logbeta), half);
        nova["fullslow"]->SetPoint(nova["fullslow"]->GetN(), pow(10, logbeta), full);
      }
    }

    #if defined(DRAWFAST) && defined(DRAWSLOW)

    {
      const double crossbeta = crossover(nova["halfslow"], nova["halffast"]);
      TGraph * fill = nova["halfboth"];
      {
        TGraph * use = nova["halfslow"];
        for(int i = 0; i < use->GetN(); i++)
          if(use->GetX()[i] < crossbeta)
            fill->SetPoint(fill->GetN(), use->GetX()[i], use->GetY()[i]);
      }

      {
        TGraph * use = nova["halffast"];
        for(int i = 0; i < use->GetN(); i++)
          if(use->GetX()[i] > crossbeta)
            fill->SetPoint(fill->GetN(), use->GetX()[i], use->GetY()[i]);
      }
    }

    {
      const double crossbeta = crossover(nova["fullslow"], nova["fullfast"]);
      TGraph * fill = nova["fullboth"];
      {
        TGraph * use = nova["fullslow"];
        for(int i = 0; i < use->GetN(); i++)
          if(use->GetX()[i] < crossbeta)
            fill->SetPoint(fill->GetN(), use->GetX()[i], use->GetY()[i]);
      }
      {
        TGraph * use = nova["fullfast"];
        for(int i = 0; i < use->GetN(); i++)
          if(use->GetX()[i] > crossbeta)
            fill->SetPoint(fill->GetN(), use->GetX()[i], use->GetY()[i]);
      }
    }
    #endif
  }

  TCanvas *can = new TCanvas;

  can->SetCanvasSize(600, 400);
  can->SetRightMargin(0.010);
  can->SetTopMargin(0.027);
  can->SetLeftMargin(0.115);
  can->SetBottomMargin(0.125);

  can->SetFrameLineWidth(2);

  can->SetLogy();

  can->SetLogx();

  const double xmin = pow(10, -4.4);
  const double ymax = 1e-11;

  TH2D dum("dum", "", 1, xmin, 1, 1, 1e-19, ymax);

  dum.Draw();

  const double alpha = 0.88;

  nova.at("halffast")->SetLineWidth(2);
  nova.at("halffast")->SetFillStyle(1001);
  nova.at("halffast")->SetLineColor(kBlack);
  nova.at("halffast")->SetFillColorAlpha(kYellow, alpha);

  nova.at("fullfast")->SetLineWidth(2);
  nova.at("fullfast")->SetFillStyle(1001);
  nova.at("fullfast")->SetLineColor(kBlack);
  nova.at("fullfast")->SetFillColorAlpha(kYellow-10, alpha);

  nova.at("halfslow")->SetLineWidth(2);
  nova.at("halfslow")->SetFillStyle(1001);
  nova.at("halfslow")->SetLineColor(kBlack);
  nova.at("halfslow")->SetFillColorAlpha(kYellow, alpha);

  nova.at("fullslow")->SetLineWidth(2);
  nova.at("fullslow")->SetFillStyle(1001);
  nova.at("fullslow")->SetLineColor(kBlack);
  nova.at("fullslow")->SetFillColorAlpha(kYellow-10, alpha);

  nova.at("halfboth")->SetLineWidth(2);
  nova.at("halfboth")->SetFillStyle(1001);
  nova.at("halfboth")->SetLineColor(kBlack);
  nova.at("halfboth")->SetFillColorAlpha(kYellow, alpha);

  nova.at("fullboth")->SetLineWidth(2);
  nova.at("fullboth")->SetFillStyle(1001);
  nova.at("fullboth")->SetLineColor(kBlack);
  nova.at("fullboth")->SetFillColorAlpha(kYellow-10, alpha);

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
  y->SetRangeUser(1e-16, ymax);

  x->SetTickSize(0.015); // Smaller than default (0.03)
  y->SetTickSize(0.015);

  TGraph slimlight;

  slimlight.SetPoint(slimlight.GetN(), 1, ymax);
  slimlight.SetPoint(slimlight.GetN(), 0.000089, ymax);
  slimlight.SetPoint(slimlight.GetN(), 0.00008988841109054591, 3.183e-14);
  slimlight.SetPoint(slimlight.GetN(), 0.00009143471406700329, 2.904e-14);
  slimlight.SetPoint(slimlight.GetN(), 0.00009380418666398144, 2.678e-14);
  slimlight.SetPoint(slimlight.GetN(), 0.00009623506263980887, 2.614e-14);
  slimlight.SetPoint(slimlight.GetN(), 0.00009915081680951927, 2.722e-14);
  slimlight.SetPoint(slimlight.GetN(), 0.00010042731504499836, 2.850e-14);
  slimlight.SetPoint(slimlight.GetN(), 0.000102154913233549,   3.057e-14);
  slimlight.SetPoint(slimlight.GetN(), 0.00010302982519024635, 3.360e-14);
  slimlight.SetPoint(slimlight.GetN(), 0.00010480219298964032, 3.683e-14);
  slimlight.SetPoint(slimlight.GetN(), 0.00010660504989847922, 4.522e-14);
  slimlight.SetPoint(slimlight.GetN(), 0.00010797751623277094, 5.037e-14);
  slimlight.SetPoint(slimlight.GetN(), 0.00010890229622637306, 5.971e-14);
  slimlight.SetPoint(slimlight.GetN(), 0.00010890229622637306, 5.971e-14);
  slimlight.SetPoint(slimlight.GetN(), 0.000109, ymax);

  slimlight.SetPoint(slimlight.GetN(), 0.009, ymax);
  slimlight.SetPoint(slimlight.GetN(), 0.009034047151333588, 2.4114e-14);
  slimlight.SetPoint(slimlight.GetN(), 0.009424638605877472, 1.7950e-14);
  slimlight.SetPoint(slimlight.GetN(), 0.009749236106412613, 1.3506e-14);
  slimlight.SetPoint(slimlight.GetN(), 0.010344413737342377, 8.9340e-15);
  slimlight.SetPoint(slimlight.GetN(), 0.011163339112402325, 5.9414e-15);
  slimlight.SetPoint(slimlight.GetN(), 0.012047095495477358, 4.4465e-15);
  slimlight.SetPoint(slimlight.GetN(), 0.013678214555560873, 3.1538e-15);
  slimlight.SetPoint(slimlight.GetN(), 0.015929637110667732, 2.5038e-15);
  slimlight.SetPoint(slimlight.GetN(), 0.017782794100389247, 2.2249e-15);
  slimlight.SetPoint(slimlight.GetN(), 0.019353731242483577, 2.0309e-15);
  slimlight.SetPoint(slimlight.GetN(), 0.02142310187038004,  1.8241e-15);
  slimlight.SetPoint(slimlight.GetN(), 0.022924196246146077, 1.7288e-15);
  slimlight.SetPoint(slimlight.GetN(), 0.026472448838535192, 1.6037e-15);
  slimlight.SetPoint(slimlight.GetN(), 0.03299000284067294,  1.4638e-15);
  slimlight.SetPoint(slimlight.GetN(), 0.03809626095526314,  1.4174e-15);
  slimlight.SetPoint(slimlight.GetN(), 0.043992875835166305, 1.3725e-15);
  slimlight.SetPoint(slimlight.GetN(), 0.057194325303422125, 1.3506e-15);
  slimlight.SetPoint(slimlight.GetN(), 0.06549022320294305,  1.3362e-15);
  slimlight.SetPoint(slimlight.GetN(), 0.07956742329247525,  1.3219e-15);
  slimlight.SetPoint(slimlight.GetN(), 0.09667053400910426,  1.3219e-15);
  slimlight.SetPoint(slimlight.GetN(), 0.14269598863094754,  1.3078e-15);
  slimlight.SetPoint(slimlight.GetN(), 0.17782794100389238,  1.3008e-15);
  slimlight.SetPoint(slimlight.GetN(), 0.23315624847200295,  1.2938e-15);
  slimlight.SetPoint(slimlight.GetN(), 0.387467512045613,    1.2938e-15);
  slimlight.SetPoint(slimlight.GetN(), 0.596671499785315,    1.2869e-15);
  slimlight.SetPoint(slimlight.GetN(), 1,                    1.2869e-15);

  slimlight.SetLineWidth(1);

  slimlight.SetFillStyle(1001);
  slimlight.SetFillColorAlpha(kGreen, alpha);

  TGraph slimheavy; // irrelevant
  slimheavy.SetPoint(slimheavy.GetN(), 0.05, 0.65e-15);
  slimheavy.SetPoint(slimheavy.GetN(), 0.00, 0.65e-15);

  slimheavy.SetLineWidth(2);

  TGraph berkeley;
  berkeley.SetPoint(berkeley.GetN(), 1, ymax);
  berkeley.SetPoint(berkeley.GetN(), 1, 1e-12);
  berkeley.SetPoint(berkeley.GetN(), 0.007, 1e-12);
  berkeley.SetPoint(berkeley.GetN(), 0.007, ymax);

  berkeley.SetPoint(berkeley.GetN(), 0.00012095437388542935, ymax);
  berkeley.SetPoint(berkeley.GetN(), 0.00012095437388542935, 9.34236263460415e-13);
  berkeley.SetPoint(berkeley.GetN(), 0.00011601205430941776, 4.5017913416809247e-13);
  berkeley.SetPoint(berkeley.GetN(), 0.00010897494036169146, 2.3031958162714385e-13);
  berkeley.SetPoint(berkeley.GetN(), 0.00010057962089638911, 1.3442545617806353e-13);
  berkeley.SetPoint(berkeley.GetN(), 0.00009175129438744992, 8.634450182087942e-14);
  berkeley.SetPoint(berkeley.GetN(), 0.00007883822934105696, 6.17793522326172e-14);
  berkeley.SetPoint(berkeley.GetN(), 0.0000669887256861911, 5.616139785470251e-14);
  berkeley.SetPoint(berkeley.GetN(), 0.000057307694675768425, 6.332872527747738e-14);
  berkeley.SetPoint(berkeley.GetN(), 0.000051155301825344104, 7.952497699458978e-14);
  berkeley.SetPoint(berkeley.GetN(), 0.00004707071661256225, 1.0729119231853001e-13);
  berkeley.SetPoint(berkeley.GetN(), 0.00004385550639408299, 1.6315378463696875e-13);
  berkeley.SetPoint(berkeley.GetN(), 0.000041365362146698854, 2.6656819649214686e-13);
  berkeley.SetPoint(berkeley.GetN(), 0.00003948583726170778, 4.2521511443329243e-13);
  berkeley.SetPoint(berkeley.GetN(), 0.000038365545318902525, 6.163107515893746e-13);
  berkeley.SetPoint(berkeley.GetN(), 0.000037505300118674406, 9.149160153048394e-13);
  berkeley.SetPoint(berkeley.GetN(), 0.000036191937950676105, 2.4418099105342384e-12);
  berkeley.SetPoint(berkeley.GetN(), 0.000036191937950676105, ymax);


  berkeley.SetLineWidth(2);
  berkeley.SetLineStyle(7);
  berkeley.SetFillStyle(1001);
  berkeley.SetFillColorAlpha(kGray, alpha);


  TGraph cabrera;
  cabrera.SetPoint(cabrera.GetN(), xmin, 7.2e-13 * 2);
  cabrera.SetPoint(cabrera.GetN(), 1   , 7.2e-13 * 2);
  cabrera.SetPoint(cabrera.GetN(), 1   , ymax);
  cabrera.SetPoint(cabrera.GetN(), xmin, ymax);

  cabrera.SetLineWidth(2);
  cabrera.SetFillColorAlpha(kViolet, alpha);

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
  macroheavy.SetPoint(macroheavy.GetN(), (1), ymax);
  macroheavy.SetPoint(macroheavy.GetN(), (xmin), ymax);

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
  macroheavy.SetFillColorAlpha(kCyan-10, alpha);


  macrolight.SetLineWidth(2);
  macrolight.SetFillStyle(1001);
  macrolight.SetLineColor(kBlack);
  macrolight.SetFillColorAlpha(kCyan, alpha);

  TGraph icecube;
  icecube.SetPoint(icecube.GetN(), (0.995), 1e-11); // just for drawing
  icecube.SetPoint(icecube.GetN(), (0.8), 1e-11); // just for drawing
  {
    std::ifstream icecubein("icecube.txt");
    double beta, lim;
    while(icecubein >> beta >> lim)
      icecube.SetPoint(icecube.GetN(), beta, lim);
  }

  icecube.SetFillStyle(1001);
  icecube.SetFillColorAlpha(kRed, alpha);
  icecube.SetLineColor(kBlack);
  icecube.SetLineWidth(2);

  TGraph antares;
  antares.SetPoint(antares.GetN(), 0.995, 1e-11); // just for drawing
  antares.SetPoint(antares.GetN(), 0.5945, 1e-11); // just for drawing
  {
    std::ifstream antaresin("antares.txt");
    double beta, lim;
    while(antaresin >> beta >> lim)
      antares.SetPoint(antares.GetN(), beta, lim);
  }

  antares.SetFillStyle(1001);
  antares.SetFillColorAlpha(kOrange, alpha);
  antares.SetLineColor(kBlack);
  antares.SetLineWidth(2);


  macroheavy.Draw("lf"); // 1e16
  macrolight.Draw("lf"); // 1e10

#if defined(DRAWFAST) && defined(DRAWSLOW)
  nova.at("fullboth")->Draw("lf");
#elif defined(DRAWFAST)
  nova.at("fullfast")->Draw("lf"); // 1e16 (same as MACRO, although we say 2e15)
#elif defined(DRAWSLOW)
  nova.at("fullslow")->Draw("lf");
#endif

  antares.Draw("lf"); // 1e10(?)
  icecube.Draw("lf"); // 1e6

#if defined(DRAWFAST) && defined(DRAWSLOW)
  nova.at("halfboth")->Draw("lf"); // 5e8
#elif defined(DRAWFAST)
  nova.at("halffast")->Draw("lf"); // 5e8
#elif defined(DRAWSLOW)
  nova.at("halfslow")->Draw("lf");
#endif

  cabrera  .Draw("lf"); // Same as NOvA (?)

  berkeley .Draw("lf"); // ~1e5, but somewhat lower altitude than SLIM

  slimlight.Draw("lf"); // 1e5


  {  // Cabrera
    const double x = 0.23,
                 y = 0.865;
    TLegend *l = new TLegend(x,      y,
                             x+0.15, y+0.12);
    l->SetTextSize(textsize);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextFont(font);
    l->AddEntry((TH1D*)NULL, "Cabrera, surface", "");
    l->Draw();
  }
  {
    const double icex = 0.72,
                 icey = 0.13;
    #ifdef THIRTEENYEAR
    TLegend *l = new TLegend(icex,      icey+textsize/2,
                             icex+0.15, icey+2.5*textsize);
    #else
    TLegend *l = new TLegend(icex,      icey,
                             icex+0.15, icey+3*textsize);
    #endif
    l->SetTextSize(textsize);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextFont(font);
    l->AddEntry((TH1D*)NULL, "IceCube,", "");
    #ifndef THIRTEENYEAR
    l->AddEntry((TH1D*)NULL, "under ice", "");
    #endif
    l->AddEntry((TH1D*)NULL, "> 10^{8}#kern[-0.3]{ }GeV", "");
    l->Draw();

    TArrow * a = new TArrow(0.36, 1.0e-18,
                            0.75, 1.0e-18, 0.011, "|>");
    a->SetLineWidth(2);
    a->Draw();
  }
  {
    const double antx = 0.67,
                 anty = 0.29;
    #ifdef THIRTEENYEAR
    TLegend *l = new TLegend(antx,      anty+textsize/2,
                             antx+0.15, anty+2.5*textsize);
    #else
    TLegend *l = new TLegend(antx,      anty,
                             antx+0.15, anty+3*textsize);
    #endif
    l->SetTextSize(textsize);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextFont(font);
    l->AddEntry((TH1D*)NULL, "ANTARES,", "");
    #ifndef THIRTEENYEAR
    l->AddEntry((TH1D*)NULL, "underwater", "");
    #endif
    l->AddEntry((TH1D*)NULL, "> 10^{10}#kern[-0.3]{ }GeV", "");
    l->Draw();

    TArrow * a = new TArrow(0.30, 3.0e-17,
                            0.59, 3.0e-17, 0.011, "|>");
    a->SetLineWidth(2);
    a->Draw();
  }
  {
    const double macrox = 0.12,
                 macroy = 0.30;
    TLegend *l = new TLegend(macrox,      macroy,
                             macrox+0.15, macroy+0.12);
    l->SetTextSize(textsize);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextFont(font);
    l->SetTextAlign(22);
    l->AddEntry((TH1D*)NULL, "MACRO", "");
    l->AddEntry(&macrolight, "> 10^{10}#kern[-0.3]{ }GeV", "");
    l->Draw();

    TArrow * a = new TArrow(1.2e-4, 6.0e-17,
                            1.2e-4, 7.0e-16, 0.011, "|>");
    a->SetLineWidth(2);
    a->Draw();
  }

  {
    const double macrox = 0.22,
                 macroy = 0.16;
    TLegend *l = new TLegend(macrox,      macroy,
                             macrox+0.15, macroy+0.12
                             );
    l->SetTextSize(textsize);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextFont(font);
    l->SetTextAlign(22);
    l->AddEntry((TH1D*)NULL, "MACRO", "");
    l->AddEntry((TH1D*)NULL, "> 10^{16}#kern[-0.3]{ }GeV", "");
    l->Draw();

    TArrow * a = new TArrow(6.00e-4, 4.0e-18,
                            2.00e-4, 2.0e-16, 0.011, "|>");
    a->SetLineWidth(2);
    a->Draw();
  }

  {
    const double slimx = 0.70,
                 slimy = 0.67;
    TLegend *l = new TLegend(slimx,      slimy,
                             slimx+0.15, slimy+textsize*3);
    l->SetTextSize(textsize);
    l->SetBorderSize(0);
    l->SetTextAlign(22);
    l->SetFillStyle(0);
    l->SetTextFont(font);
    l->AddEntry((TH1D*)NULL, "SLIM,", "");
    l->AddEntry((TH1D*)NULL, "mountaintop", "");
    l->AddEntry(&slimlight, "> 10^{5}#kern[-0.3]{ }GeV", "");
    l->Draw();
  }
  { // Berkeley 1983 with Price 1984 re-analysis, as per Groom
    const double x = 0.10,
                 y = 0.59;
    TLegend *l = new TLegend(x,      y,
                             x+0.15, y+textsize);
    l->SetTextSize(textsize);
    l->SetBorderSize(0);
    l->SetTextAlign(22);
    l->SetFillStyle(0);
    l->SetTextFont(font);
    l->AddEntry((TH1D*)NULL, "Berkeley", "");
    l->Draw();

    TArrow * a = new TArrow(0.9e-4, 1.0e-14,
                            0.7e-4, 4.5e-14, 0.011, "|>");
    a->SetLineWidth(2);
    a->Draw();
  }

  #if defined(DRAWFAST) && defined(DRAWSLOW)
  {
    const double thisworkx = 0.26;
    const double thisworky = 0.54;
    TLegend *l = new TLegend(thisworkx,      thisworky,
                             thisworkx+0.15, thisworky +0.06*5
    );
    l->SetTextSize(textsize);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextFont(font);
    l->SetTextAlign(12);
    #ifdef THIRTEENYEAR
    l->AddEntry((TH1D*)NULL, "", "");
    l->AddEntry((TH1D*)NULL, "Sensitivity", "");
    l->AddEntry((TH1D*)NULL, "11 live-year", "");
    #else
    l->AddEntry((TH1D*)NULL, "NOvA", "");
    l->AddEntry((TH1D*)NULL, "7.5 (slow), 7.4 (fast)", "");
    l->AddEntry((TH1D*)NULL, "live-year,", "");
    #endif
    l->AddEntry((TH1D*)NULL, ">10^{9}#kern[-0.3]{ }GeV", "");
    l->AddEntry((TH1D*)NULL, "  >10^{7}#kern[-0.3]{ }GeV for #beta > 0.1", "");
    l->Draw();
  }
  {
    const double novaheavyx = 0.46;
    const double novaheavyy = 0.21;
    TLegend *l = new TLegend(novaheavyx,      novaheavyy,
                             novaheavyx+0.15, novaheavyy+textsize*3);
    l->SetTextSize(textsize);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextFont(font);
    l->SetTextAlign(22);
    l->SetTextColor(kBlack);
    #ifdef THIRTEENYEAR
    l->AddEntry((TH1D*)NULL, "Sensitivity", "");
    #else
    l->AddEntry((TH1D*)NULL, "NOvA", "");
    #endif
    l->AddEntry((TH1D*)NULL, "> 10^{15}#kern[-0.3]{ }GeV", "");
    l->AddEntry((TH1D*)NULL, "> 10^{13}#kern[-0.3]{ }GeV for #beta > 0.1", "");
    l->Draw();

    TArrow * a = new TArrow(1.0e-2, 2.5e-17,
    #ifdef THIRTEENYEAR
                            1.0e-2, 0.7e-16, 0.011, "|>");
    #else
                            1.0e-2, 1.1e-16, 0.011, "|>");
    #endif
    a->SetLineColor(kBlack);
    a->SetFillColor(kBlack);
    a->SetLineWidth(2);
    a->Draw();
  }
  #elif defined(DRAWFAST)
  {
    const double thisworkx = 0.390;
    const double thisworky = 0.515;
    TLegend *l = new TLegend(thisworkx,      thisworky,
                             thisworkx+0.15, thisworky +0.06*5*0.8
    );
    l->SetTextSize(textsize*0.8);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextFont(font);
    l->SetTextAlign(12);
    l->AddEntry((TH1D*)NULL, "NOvA", "");
    l->AddEntry((TH1D*)NULL, "7.4 live-year", "");
    l->AddEntry((TH1D*)NULL, "  fast-only", "");
    l->AddEntry
      ((TH1D*)NULL,
      "        >10^{9}#kern[-0.3]{ }GeV",
    "");
    l->AddEntry((TH1D*)NULL, "                  >10^{7}#kern[-0.3]{ }GeV for #beta > 0.1", "");
    l->Draw();
  }
  {
    const double novaheavyx = 0.46;
    const double novaheavyy = 0.21;
    TLegend *l = new TLegend(novaheavyx,      novaheavyy,
                             novaheavyx+0.15, novaheavyy+textsize*3);
    l->SetTextSize(textsize);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextFont(font);
    l->SetTextAlign(22);
    l->SetTextColor(kBlack);
    l->AddEntry((TH1D*)NULL, "NOvA", "");
    l->AddEntry((TH1D*)NULL, "> 10^{15}#kern[-0.3]{ }GeV", "");
    l->AddEntry((TH1D*)NULL, "> 10^{13}#kern[-0.3]{ }GeV for #beta > 0.1", "");
    l->Draw();

    TArrow * a = new TArrow(1.0e-2, 2.5e-17,
    #ifdef THIRTEENYEAR
                            1.0e-2, 0.7e-16, 0.011, "|>");
    #else
                            1.0e-2, 1.1e-16, 0.011, "|>");
    #endif
    a->SetLineColor(kBlack);
    a->SetFillColor(kBlack);
    a->SetLineWidth(2);
    a->Draw();
  }
  #elif defined(DRAWSLOW)
  {
    const double thisworkx = 0.3;
    const double thisworky = 0.55;
    TLegend *l = new TLegend(thisworkx,      thisworky,
                             thisworkx+0.15, thisworky +textsize*4
    );
    l->SetTextSize(textsize);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextFont(font);
    l->SetTextAlign(12);
    l->AddEntry((TH1D*)NULL, "NOvA", "");
    l->AddEntry((TH1D*)NULL, "7.5 live-year", "");
    l->AddEntry((TH1D*)NULL, "slow-only", "");
    l->AddEntry((TH1D*)NULL, ">10^{9}#kern[-0.3]{ }GeV", "");
    l->Draw();
  }
  {
    const double novaheavyx = 0.42,
                 novaheavyy = 0.26;
    TLegend *l = new TLegend(novaheavyx,      novaheavyy,
                             novaheavyx+0.15, novaheavyy+textsize*2);
    l->SetTextSize(textsize);
    l->SetBorderSize(0);
    l->SetFillStyle(0);
    l->SetTextFont(font);
    l->SetTextAlign(22);
    l->SetTextColor(kBlack);
    l->AddEntry((TH1D*)NULL, "NOvA", "");
    l->AddEntry((TH1D*)NULL, "> 10^{16}#kern[-0.3]{ }GeV", "");
    l->Draw();

    TArrow * a = new TArrow(4.0e-3, 2.0e-17,
                            4.0e-3, 1.1e-16, 0.011, "|>");
    a->SetLineColor(kBlack);
    a->SetFillColor(kBlack);
    a->SetLineWidth(2);
    a->Draw();
  }
  #endif

  can->RedrawAxis();

  std::string speed;
  #ifdef DRAWFAST
    #ifdef DRAWSLOW
      speed = "slowfast";
    #else
      speed = "fast";
    #endif
  #elif DRAWSLOW
    speed = "slow";
  #endif


  can->SaveAs(Form("limit_sensitivity_%s_plot.pdf", speed.c_str()));
}
