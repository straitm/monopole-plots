#ifndef GUARD_EVENT_LIST_MONO_HH
#define GUARD_EVENT_LIST_MONO_HH

/*
  This class maintains a list of events (run/subrun/event) 
  and reports duplicate events.
 */

#include <map>
#include <set>

class Event_List
{
public:
  Event_List();
  Event_List(std::string type);

  bool is_duplicate(int run, int subrun, int event);
  int  n_duplicate_events() const;
  
private:
  bool ignore_events_;
  int  n_duplicate_events_;
  
  std::map<int, std::map<int, std::set<int> > > runs_;
};

#endif
