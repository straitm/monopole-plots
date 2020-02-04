#include "Event_Info.hh"

#include <algorithm>
#include <TMath.h>

#include "Constants.hh"
#include "MFRoot.hh"

using MFRoot::tree_t;

Event_Info::Event_Info(tree_t const& t) :
  t_(t),
  beta_long_name_("invalid"),
  beta_short_name_("invalid")
{
}

Event_Info::Event_Info(tree_t const& t, const std::string & input_file_name) :
  t_(t),
  beta_long_name_("invalid"),
  beta_short_name_("invalid")
{
  if (is_mc()) extract_beta_value(input_file_name);
}

double Event_Info::area_projected() const
{
  double result = 0.0;

  for (auto const& area_normal : AREA_NORMALS)
  {
    TVector3 A = area_normal.second;
    TVector3 v = mc_velocity_unit_vector();

    // The monopole is going into the detector and therefore has a negative
    // dot product with the outward pointing surface normal.
    if (A.Dot(v) < 0)
      result += std::fabs(A.Dot(v));
  }

  return result;
}

double Event_Info::area_total() const
{
  double result = 0.0;

  for (auto const& area_normal : AREA_NORMALS) {
    TVector3 A = area_normal.second;
    result += A.Mag();
  }

  return result;
}

double Event_Info::beta() const
{
  if (!is_preselected_reco()) {
    std::cerr << "\n\n\t\tThis event contains no tracks that pass the "
	      << "reco preselection." << std::endl;
    _exit(1);
  }
  
  return get("beta", first_preselected_reco_track());
}



double Event_Info::beta_from_file_name() const
{
  if (beta_long_name_ == "invalid") {
    std::cerr << "\n\n\t\tThe beta_long_name_ value has not been assigned."
	      << std::endl;
    _exit(1);
  }

  double result = std::stod(beta_long_name_);

  if (mc_beta() > 0) {
    double beta_difference = std::abs(mc_beta() - result) / mc_beta();

    if (beta_difference > 0.01) {
      std::cerr << "\n\n\t\tThe beta values differ by more than 1%:"
		<< "beta_values: MC, from name, relative difference = "
		<< mc_beta() << ", " << result << ", " << beta_difference
		<< std::endl;
      _exit(1);
    }
  }

  return result;
}

std::string Event_Info::beta_name() const
{
  if (beta_short_name_ == "invalid")
  {
    std::cerr << "\n\n\t\tThe beta_short_name_ value has not been assigned."
	      << std::endl;
    _exit(1);
  }

  return beta_short_name_;
}

bool Event_Info::beta_value_matches(const std::string & beta_value) const
{
  return !is_mc() || beta_name() == beta_value;
}

int Event_Info::event() const
{
  return get("event_number");
}

int Event_Info::first_preselected_reco_track() const
{
  for (int n_track = 1; n_track <= n_tracks_min(); ++n_track)
    if (is_preselected_reco_track(n_track))
      return n_track;

  return 0;
}

void Event_Info::extract_beta_value(const std::string & name)
{
  // sample file name:
  // ...monopole_merge_gain100_n1000_beta2.51188643e-03_birksmodC_64.root

  std::size_t it = name.find("_beta");
  if (it == std::string::npos)
    it = name.find("_Beta");
  
  std::size_t start = it + 5;

  std::size_t end = it + 19;

  if (it != std::string::npos &&
      end != std::string::npos)
  {
    beta_long_name_ = name.substr(start, end - start);
    beta_short_name_ = name.substr(start, 3) + name.substr(end - 4, 4);
  } else {
    std::cerr << "\n\n\t\tInput file name branch value not formatted properly:"
	      << "\n\t\t  " << name << "\n" << std::endl;
    _exit(1);
  }
}

double Event_Info::get(std::string branch_name) const
{
  return MFRoot::get(t_, branch_name);
}

double Event_Info::get(std::string parameter, unsigned n_track) const
{
  std::string branch_name =
    "reco_track_" + std::to_string(n_track) + "_" + parameter;

  return get(branch_name);
}

double Event_Info::gap_max() const
{
  int n_track = first_preselected_reco_track();

  double dt = get("dt", n_track);
  if (dt == 0){
    fprintf(stderr, "\n\n\tThis event has zero time duration!\n\n");
    _exit(1);
  }
  
  double gap_xz = get("max_time_gap_xz", n_track) / dt;
  double gap_yz = get("max_time_gap_yz", n_track) / dt;

  return std::max(gap_xz, gap_yz);
}

double Event_Info::gen_beta() const
{
  double result = get("gen_beta");

  // factor out the incorrect mass used for the initial beta (= p/m)
  // calculation and use the correct monopole mass
  result *= GENERATOR_INCORRECT_MASS / MONOPOLE_MASS;

  double gen_mc_diff = (result - beta_from_file_name()) / beta_from_file_name();
  if (gen_mc_diff > 1e-2)
    std::cout << "gen_beta       = " << result
	      << "\nbeta_from_file = " << beta_from_file_name()
	      << "\ngen_mc_diff    = " << gen_mc_diff
	      << std::endl;
  
  return result;
}

bool Event_Info::high_mc_fraction() const
{
  return mc_fraction() > MC_FRACTION_CUT;
}

bool Event_Info::high_r2_min() const
{
  return r2_min() > R2_MIN_CUT;
}

bool Event_Info::hit_detector() const
{
  return get("mc_monopole_hit_detector") == 1;
}

bool Event_Info::is_data() const
{
  return !is_mc();
}

bool Event_Info::is_linear() const
{
  return is_preselected_reco() && high_r2_min() && low_gap_max();
}

bool Event_Info::is_mc() const
{
  return get("is_mc");
}

bool Event_Info::is_preselected_mc() const
{
  if (get("mc_n_hits_x") < 20) return false;
  if (get("mc_n_hits_y") < 20) return false;
  if (get("mc_dplane_x") < 10) return false;
  if (get("mc_dplane_y") < 10) return false;
  if (get("mc_length") < 1000) return false;

  return true;
}

bool Event_Info::is_preselected_reco() const
{
  if (n_tracks() <= 0) return false;

  for (int n_track = 1; n_track <= n_tracks_min(); ++n_track)
    if (is_preselected_reco_track(n_track))
      return true;

  return false;
}

bool Event_Info::is_preselected_reco_track(unsigned n_track) const
{
  if (get("n_hits_x", n_track) < 20) return false;
  if (get("n_hits_y", n_track) < 20) return false;
  if (get("dplane_x", n_track) < 10) return false;
  if (get("dplane_y", n_track) < 10) return false;
  if (get("length", n_track) < 1000) return false;

  return true;
}


bool Event_Info::is_signal() const
{
  return hit_detector() && is_triggered() && is_linear() && is_slow();
}

bool Event_Info::is_slow() const
{
  return beta() <= SLOW_BETA_CUT;
}

bool Event_Info::is_triggered() const
{
  int n_triggers = get("ddt_trigger_decisions");

  if (n_triggers > 1)
    std::cout << "More than one trigger decision!"
	      << "\n\tddt_trigger_decisions = " << n_triggers
	      << std::endl;

  return n_triggers >= 1;
}


bool Event_Info::low_gap_max() const
{
  return gap_max() <= GAP_MAX_CUT;
}

double Event_Info::mc_beta() const
{
  return get("mc_beta");
}

double Event_Info::mc_fraction() const
{
  double result = 0.0;
  for (unsigned n = 1; n <= MAX_RECO_N_TRACKS; ++n)
    result = std::max(result, mc_fraction(n));

  return result;
}

double Event_Info::mc_fraction(unsigned n_track) const
{
  if (!hit_detector()) return 0.0;
  
  double total_mc_hits = get("mc_n_hits_all");
  if (total_mc_hits <= 0) return 0.0;

  double mc_hits_on_track = get("n_mc_hits", n_track);
  if (mc_hits_on_track <= 0) return 0.0;
  
  return mc_hits_on_track / total_mc_hits;
}

double Event_Info::mc_theta(std::string view) const
{
  return theta("mc", view);
}

double Event_Info::mc_velocity() const
{
  return get("mc_velocity");
}

TVector3 Event_Info::mc_velocity_vector() const
{
  double vx = get("mc_dxdt");
  double vy = get("mc_dydt");
  double vz = get("mc_dzdt");

  return TVector3(vx, vy, vz);
}

TVector3 Event_Info::mc_velocity_unit_vector() const
{
  TVector3 result = mc_velocity_vector();

  // normalize vector:
  result *= (1.0 / result.Mag());

  return result;
}

int Event_Info::n_tracks() const
{
  return get("reco_n_tracks");
}

int Event_Info::n_tracks_min() const
{
  return std::min(n_tracks(), MAX_RECO_N_TRACKS);
}

double Event_Info::r2_min() const
{
  int n_track = first_preselected_reco_track();

  double r2_xt = get("r2_xt", n_track);
  double r2_yt = get("r2_yt", n_track);

  return std::min(r2_xt, r2_yt);
}

double Event_Info::reco_theta(std::string view, unsigned n_track) const
{
  return theta("reco_track_" + std::to_string(n_track), view);
}

int Event_Info::run() const
{
  return get("run_number");
}

int Event_Info::subrun() const
{
  return get("subrun_number");
}

double Event_Info::theta(int n_track, std::string view) const
{
  std::string prefix = "reco_track_" + std::to_string(n_track);

  return theta(prefix, view);
}

double Event_Info::theta(std::string prefix, std::string view) const
{
  std::string branch_name;
  if (view == "xz") branch_name = prefix + "_dxdz";
  else if (view == "yz") branch_name = prefix + "_dydz";
  else{
    fprintf(stderr, "\n\n\t\tInvalid View Value!\n");
    _exit(1);
  }

  double dvdz = get(branch_name);

  return theta(dvdz);
}

double Event_Info::theta(double dvdz) const
{
  return std::atan(dvdz) * 180 / TMath::Pi();
}
