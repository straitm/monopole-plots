#include "Subrun_Info.hh"



double Subrun_Info::duration() const
{
  return get("subrun_duration_in_seconds");
}



int Subrun_Info::run() const
{
  return get("run_number");
}



int Subrun_Info::subrun() const
{
  return get("subrun_number");
}
