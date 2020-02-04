#include "Event_List.hh"

#include "Constants.hh"

#include <iostream>



Event_List::Event_List() : ignore_events_(false), n_duplicate_events_(0)
{
}



Event_List::Event_List(std::string type) : n_duplicate_events_(0)
{
  // MC events are not numbered uniquely because all of the files use the same
  // run number.  We will therefore ignore any non-data events.
  
  ignore_events_ = true;
  if (type == DATA_SAMPLE_NAME)
    ignore_events_ = false;
}



bool Event_List::is_duplicate(int run, int subrun, int event)
{
  if (ignore_events_)
    return false;

  bool found_duplicate = false;
  
  auto r_it = runs_.find(run);
  if (r_it != runs_.end())
  {
    auto & subruns = r_it->second;
    auto s_it = subruns.find(subrun);
    if (s_it != subruns.end())
    {
      auto & events = s_it->second;
      auto e_it = events.find(event);
      if (e_it != events.end())
	found_duplicate = true;
    }
  }

  if (!found_duplicate)
    (runs_[run])[subrun].insert(event);
  else
    ++n_duplicate_events_;

  return found_duplicate;
}



int Event_List::n_duplicate_events() const
{
  return n_duplicate_events_;
}
