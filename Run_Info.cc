#include "Run_Info.hh"

#include "MFRoot.hh"

using MFRoot::tree_t;



#ifndef __CINT__

#include <boost/format.hpp>



Run_Info::Run_Info()
{}



Run_Info::Run_Info(tree_t const& t) : t_(t)
{}



void Run_Info::add(Subrun_Info const& info)
{
  subruns_.emplace(info.subrun(), info);
}



double Run_Info::get(std::string branch_name) const
{
  return MFRoot::get(t_, branch_name);
}



void Run_Info::print() const
{
  double sum_subrun_duration = 0;
  for (auto const& subrun : subruns_)
    sum_subrun_duration += subrun.second.duration();

  boost::format form("%5u %7u %14u %18.1f %20.1f");
  std::cout <<
    form %
    run() %
    get("run_n_subruns") %
    subruns_.size() %
    get("run_duration_in_seconds") %
    sum_subrun_duration
	    << std::endl;
}



int Run_Info::run() const
{
  return get("run_number");
}



subruns_t Run_Info::subruns() const
{
  return subruns_;
}



#endif
