#ifndef DATA_CHECK_RUN_INF_HH
#define DATA_CHECK_RUN_INF_HH

/*
  A class to hold the run-/subrun-level information from the data check trees.
 */

#include "MFRoot.hh"
#include "Subrun_Info.hh"

#include <map>

typedef std::map<int, Subrun_Info> subruns_t;

class Run_Info
{
public:
  Run_Info();
  Run_Info(MFRoot::tree_t const& t);

  void add(Subrun_Info const& info);
  void print() const;
  int  run() const;

  subruns_t subruns() const;
  
  double get(std::string branch_name) const;

  
private:
  MFRoot::tree_t t_;  
  subruns_t subruns_;
};


#endif
