#ifndef MONO_CONSTANTS_ROOT_HH
#define MONO_CONSTANTS_ROOT_HH

#include <TVector3.h>

#include <map>
#include <vector>

namespace Constants
{
  const std::string DATA_CHECK_FILE =
    "~/data/mono/190220_DataCheck_Tree_Full.root";
  // "~/data/mono/190215_DataCheck_Tree.root";
  // "~/data/mono/170802_DataCheck_Tree.root";
  
  const std::string DATA_RECO_FILE =
    // 10% training sample:
    // "~/data/mono/170807_Monopole_Tree_DDslowmono.root";
    // 100% data sample:
    "data/180913_Monopole_Tree_DDslowmono.root";

  const std::string DATA_QUALITY_FILE =
    "/Users/mfrank/data/mono/170802_Good_Runs.csv";
  
  const std::string MC_RECO_FILE_LOW_ENERGY =
    "~/data/mono/180629_Monopole_Tree_MC_LowE.root";

  const std::string MC_RECO_FILE_HIGH_ENERGY =
    "~/data/mono/180629_Monopole_Tree_MC_HighE.root";

  const std::string MC_RECO_FILE =
    "data/191119_Monpole_Tree_MC_0.9dEdx.root";

  const std::string DATA_TRIG_FILE =
    "~/data/mono/150528_Monopole_Hist_Data_Gap30_Sparse10_withRate.root";

  const std::string MC_TRIG_FILE =
    "~/data/mono/150623_Monopole_Hist_MC.root";

  
  const std::string DATA_SAMPLE_NAME = "Data";

  const std::string MC_BETA_NAME      = "1.0e-03";
  const std::string MC_BETA_NICE_NAME = "10^{-3}";
  // const std::string MC_BETA_NAME      = "3.1e-04";
  // const std::string MC_BETA_NICE_NAME = "3#times10^{-4}";
  // const std::string MC_BETA_NAME      = "1.0e-04";
  // const std::string MC_BETA_NICE_NAME = "10^{-4}";
  // const std::string MC_BETA_NAME      = "5.0e-03";
  // const std::string MC_BETA_NICE_NAME = "5#times10^{-3}";
  // const std::string MC_BETA_NAME      = "5.0e-04";
  // const std::string MC_BETA_NICE_NAME = "5#times10^{-4}";


  
  const double DEFAULT_VALUE = -9999;
  
  const int MAX_RECO_N_TRACKS = 3.0;

  const double MC_FRACTION_CUT = 0.5;
  const double R2_MIN_CUT      = 0.95;
  const double GAP_MAX_CUT     = 0.2;

  const double SAMPLES_PER_MICROSECOND   =  2.0; // one every 500 ns
  const double TDC_TICKS_PER_MICROSECOND = 64.0;

  const double N_SUBRUNS_PER_RUN = 64.0;

  const double SLOW_BETA_CUT = 1e-2;

  const double GENERATOR_INCORRECT_MASS = -2.14748e9; // GeV / c^2
  const double MONOPOLE_MASS            =  1e16;      // GeV / c^2

  const double SPEED_OF_LIGHT  = 30.0; // cm / ns

  const double LIVE_TIME               = 8.213e6; // s
  const double LIVE_TIME_DDT_CORRECTED = 8.202e6; // s
  
  // Area Normal Vectors (units are m^2)
  // * from DocDB 13,623 (total = 4,173.2 m^2)
  // std::map<std::string, TVector3> AREA_NORMALS =
  //   {
  //     {"front" , {   0.0,    0.0, -240.0} } ,
  //     {"back"  , {   0.0,    0.0,  240.0} } ,
  //     {"east"  , {-923.3,    0.0,    0.0} } ,
  //     {"west"  , { 923.3,    0.0,    0.0} } ,
  //     {"top"   , {   0.0,  923.3,    0.0} } ,
  //     {"bottom", {   0.0, -923.3,    0.0} }
  //   };
  // * from wiki, consistent with my tech note (total = 4,082.48 m^2)
  // https://cdcvs.fnal.gov/redmine/projects/novaart/wiki/Detector_Specifications) 
  const double FD_DX = 7.65 + 7.58; // m
  const double FD_DY = 7.65 + 7.49; // m
  const double FD_DZ = 59.62;       // m

  const double A_BACK = FD_DX * FD_DY; // 230.582 m^2
  const double A_WEST = FD_DY * FD_DZ; // 902.647 m^2
  const double A_TOP  = FD_DX * FD_DZ; // 908.013 m^2

#if 0
  std::map<std::string, TVector3> AREA_NORMALS =
    {
      {"front" , {          0.0,          0.0, -1.0 * A_BACK } } ,
      {"back"  , {          0.0,          0.0,        A_BACK } } , 
      {"east"  , {-1.0 * A_WEST,          0.0,            0.0} } ,
      {"west"  , {       A_WEST,          0.0,            0.0} } ,
      {"top"   , {          0.0,        A_TOP,            0.0} } ,
      {"bottom", {          0.0, -1.0 * A_TOP,            0.0} }
    };
#endif
  
  
  const unsigned int nbranch = 88;
  const char * BRANCH_NAMES[nbranch] =
    {
      "ddt_trigger_decisions",
      "event_number",
      "gen_dxdz",
      "gen_dydz",
      "gen_beta",
      "is_mc",
      "mc_adc_mean",
      "mc_beta",
      "mc_dEdx",
      "mc_dEdx_Birks",
      "mc_dplane_x",
      "mc_dplane_y",
      "mc_dplane_all",
      "mc_dx",
      "mc_dxdt",
      "mc_dxdz",
      "mc_dy",
      "mc_dydt",
      "mc_dydz",
      "mc_dz",
      "mc_dzdt",
      "mc_length",
      "mc_monopole_hit_detector",
      "mc_n_hits_all",
      "mc_n_hits_all_without_noise",
      "mc_n_hits_x",
      "mc_n_hits_y",
      "mc_n_hits_ge_100_ADC",
      "mc_n_hits_ge_150_ADC",
      "mc_velocity",
      "reco_n_tracks",
      "reco_track_1_beta",
      "reco_track_1_dplane_x",
      "reco_track_1_dplane_y",
      "reco_track_1_dt",
      "reco_track_1_dx",
      "reco_track_1_dxdz",
      "reco_track_1_dy",
      "reco_track_1_dydz",
      "reco_track_1_dz",
      "reco_track_1_length",
      "reco_track_1_max_time_gap_xz",
      "reco_track_1_max_time_gap_yz",
      "reco_track_1_n_hits",
      "reco_track_1_n_hits_x",
      "reco_track_1_n_hits_y",
      "reco_track_1_n_mc_hits",
      "reco_track_1_r2_xt",
      "reco_track_1_r2_yt",
      "reco_track_2_beta",
      "reco_track_2_dplane_x",
      "reco_track_2_dplane_y",
      "reco_track_2_dt",
      "reco_track_2_dx",
      "reco_track_2_dxdz",
      "reco_track_2_dy",
      "reco_track_2_dydz",
      "reco_track_2_dz",
      "reco_track_2_length",
      "reco_track_2_max_time_gap_xz",
      "reco_track_2_max_time_gap_yz",
      "reco_track_2_n_hits",
      "reco_track_2_n_hits_x",
      "reco_track_2_n_hits_y",
      "reco_track_2_n_mc_hits",
      "reco_track_2_r2_xt",
      "reco_track_2_r2_yt",
      "reco_track_3_beta",
      "reco_track_3_dplane_x",
      "reco_track_3_dplane_y",
      "reco_track_3_dt",
      "reco_track_3_dx",
      "reco_track_3_dxdz",
      "reco_track_3_dy",
      "reco_track_3_dydz",
      "reco_track_3_dz",
      "reco_track_3_length",
      "reco_track_3_max_time_gap_xz",
      "reco_track_3_max_time_gap_yz",
      "reco_track_3_n_hits",
      "reco_track_3_n_hits_x",
      "reco_track_3_n_hits_y",
      "reco_track_3_n_mc_hits",
      "reco_track_3_r2_xt",
      "reco_track_3_r2_yt",
      "run_number",
      "slice_sum_adc_max",
      "subrun_number"
    };

  const unsigned int ndcbranch = 19;
  const char * DATA_CHECK_BRANCH_NAMES[ndcbranch] =
    {
      "run_duration_in_seconds",
      "run_farm70_in_configuration",
      "run_n_active_channels",
      "run_n_active_dcms",
      "run_n_active_diblocks",
      "run_n_buffer_nodes",
      "run_n_subruns",
      "run_n_triggers",      
      "run_number",
      "run_partition",
      "subrun_number",
      "subrun_duration_in_seconds",
      "subrun_n_events",
      "subrun_n_events_empty",
      "subrun_n_events_slowmono",
      "subrun_n_hits",
      "subrun_prescale_max",
      "subrun_prescale_mean",
      "subrun_prescale_total",
    };
}

#endif
