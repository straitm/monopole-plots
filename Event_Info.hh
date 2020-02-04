#ifndef MONO_ROOT_EVENT_INFO_HH
#define MONO_ROOT_EVENT_INFO_HH

/*
  This class will serve as a central object to determine derived
  event quantities, so they are calculated in the same way in all
  root macros.
 */

#include "MFRoot.hh"

#include <TVector3.h>

class Event_Info
{
public:
  Event_Info(MFRoot::tree_t const& t);
  Event_Info(MFRoot::tree_t const& t, const std::string & input_file_name);

  bool beta_value_matches(const std::string & beta_value) const;
  bool high_mc_fraction() const;
  bool high_r2_min() const;
  bool hit_detector() const;
  bool is_data() const;
  bool is_linear() const;
  bool is_mc() const;
  bool is_preselected_mc() const;
  bool is_preselected_reco() const;
  bool is_preselected_reco_track(unsigned n_track) const;
  bool is_signal() const;
  bool is_slow() const;
  bool is_triggered() const;
  bool low_gap_max() const;
  
  double      area_total() const;
  double      area_projected() const;
  double      beta() const;
  double      beta_from_file_name() const;
  std::string beta_name() const;
  int         event() const;
  int         first_preselected_reco_track() const;
  double      gap_max() const;
  double      gen_beta() const;
  double      mc_beta() const;
  double      mc_fraction() const;
  double      mc_fraction(unsigned n_track) const;
  double      mc_theta(std::string view) const;
  double      mc_velocity() const;
  TVector3    mc_velocity_vector() const;
  TVector3    mc_velocity_unit_vector() const;
  int         n_tracks() const;
  int         n_tracks_min() const;
  double      r2_min() const;
  double      reco_theta(std::string view, unsigned n_track) const;
  int         run() const;
  int         subrun() const;
  double      theta(int n_track, std::string view) const;
  double      theta(std::string prefix, std::string view) const;
  double      theta(double dvdz) const;

  double get(std::string branch_name) const;
  double get(std::string parameter, unsigned n_track) const;  

  
private:
  void extract_beta_value(const std::string & name);
  
  MFRoot::tree_t t_;

  std::string beta_long_name_;
  std::string beta_short_name_;
};

#endif
