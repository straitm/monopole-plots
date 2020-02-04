#ifndef MONO_MFROOT_INCLUDE_GUARD
#define MONO_MFROOT_INCLUDE_GUARD

/*
  The MFRoot namespace contains several convenience functions for Root.
  Their main advantage is that they perform error checking and output
  nicely formatted error output in case something is invalid 
  (e.g. a non-existing histogram name).
*/

#include <TCanvas.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2.h>

#include <iostream>
#include <map>


namespace MFRoot
{
  typedef std::map<std::string, double> tree_t;
  typedef std::map<std::string, TH1*> hist_t;

  typedef std::map<std::string, tree_t> trees_t;
  typedef std::map<std::string, hist_t> hists_t;

  class Canvas_Manager;
  class Tree_Info;

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

  double calculate_integral(TH1* hist, double x_min, double x_max);
  void copy(hist_t & h, std::string from, std::string to);
  TH1* divide_with_poisson_errors(TH1* numerator, TH1* denominator,
				  std::string const& result_hist_name);
  void fill(tree_t const& t, hist_t const& h, std::string const name);


  TH1* get(TFile* file, std::string hist_name);
  TH2* get2D(TFile* file, std::string hist_name);
  TH1* integral_plot(TH1 *hist);
  TFile* open(std::string file_name);
  void resize_stats_box(TH1* hist, double const x1, double const x2,
			double const y1, double const y2);
  void wilson_interval(double const& k, double const& n,
		       double & center, double & error);
  TH1* wilson_interval(TH1* numerator, TH1* denominator,
		       std::string const& result_hist_name);
};


class MFRoot::Canvas_Manager
{
public:
  Canvas_Manager();

  // format = 1 is rectangular canvas
  // format = 2 is square canvas
  void new_canvas(std::string canvas_name, 
		  bool logx = false, bool logy = false,
		  bool gridx = false, bool gridy = false,
		  int format = 1);

  void add_canvas(TCanvas *canvas, std::string canvas_name);
  
  TCanvas* get(std::string canvas_name);
  
  void save(std::string canvas_name, 
	    std::string path = "", std::string extension = "pdf");
  void save_all(std::string path);
  void save_all(std::string path, std::string extension);

private:
  std::map<std::string, TCanvas*> canvases_;
  unsigned n_canvas_;
};


class MFRoot::Tree_Info
{
public:
  Tree_Info();
  Tree_Info(tree_t const& t);

  double get(std::string branch_name) const;

  
private:
  tree_t t_;  
};

  
namespace nova
{
  void Preliminary();
  void Simulation();
};


#endif
