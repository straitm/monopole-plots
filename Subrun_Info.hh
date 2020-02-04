#ifndef DATA_CHECK_SUBRUN_INFO
#define DATA_CHECK_SUBRUN_INFO

/*
  A class to hold the subrun-level information from the data check trees.
 */

#include "MFRoot.hh"


class Subrun_Info : public MFRoot::Tree_Info
{
public:
  // Subrun_Info(MFRoot::tree_t const& t) : t_(t) {};
  using MFRoot::Tree_Info::Tree_Info;
  
  double duration() const;
  int    run() const;
  int    subrun() const;

// private:
//   MFRoot::tree_t t_;
};


#endif
