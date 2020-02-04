#ifndef GUARD_DATA_QUALITY_MONO
#define GUARD_DATA_QUALITY_MONO

/*
 Class to parse the data quality information from the database
 passed in through a CSV file (generated from the DB information).
*/

#include <map>

class Data_Quality
{
public:
  Data_Quality() = default;
  Data_Quality(std::string line);

  std::string comment() const;
  bool is_good() const;
  bool is_good_with_14_diblocks() const;
  int flag() const;
  int n_good_diblocks() const;
  int run() const;
  int subrun() const;
  
private:
  std::string comment_;
  int flag_;
  int n_good_diblocks_;
  int run_;
  int subrun_;
};

typedef std::map<int, std::map<int, Data_Quality> > dq_t;

dq_t parse_data_quality_csv_file();


#endif
